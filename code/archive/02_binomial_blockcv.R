# remotes::install_github("inlabru-org/fmesher", ref = "stable")
library(sf)
library(sp)
library(sdmTMB)
library(lubridate)
library(dplyr)
library(future)
library(INLA)
library(inlabru)
library(viridis)
plan(multisession)
haul <- readRDS("data/haul_cleaned.rds")
catch <- readRDS("data/catch_cleaned.rds")
haul$trawl_id <- as.numeric(haul$trawl_id)
# initial loop over the cutoff values
df <- expand.grid(
  "cutoff" = round(exp(seq(log(12), log(175), length.out=8))),
  "folds" = c(10),
  "bin_width" = seq(10,150,by=20),
  "species" = unique(catch$common_name),
  "n" = NA,
  "blocks"=NA,
  "present_dens_ll" = NA,
  "present_sd_ll" = NA,
  "present_converged" = NA,
  "range" = NA,
  "sigma_O" = NA,
  "sigma_E" = NA
  # "total_dens_ll" = NA,
  # "total_sd_ll" = NA,  
  # "total_converged" = NA,
  #"selection" = c("random", "systematic")
)
# df <- dplyr::filter(df, !(range==25 & selection == "systematic"),
#                     !(range==125 & selection == "systematic"))
df <- dplyr::arrange(df, -cutoff, species, bin_width)

set.seed(2021)

for (i in 1:nrow(df)) {
  
  catch_sub <- dplyr::filter(catch, common_name == df$species[i])
  
  # Join catch and haul data
  joined_dat <- haul |>
    dplyr::filter(trawl_id %in% catch$trawl_id) |>
    left_join(catch_sub[,c("date","cpue_kg_km2","depth_m","trawl_id")], 
              by = "trawl_id") |>
    dplyr::filter(!is.na(cpue_kg_km2)) |>
    dplyr::mutate(present = ifelse(cpue_kg_km2 > 0, 1, 0))
  joined_dat$date <- lubridate::ymd(as.POSIXct.Date(joined_dat$date))
  joined_dat$year <- year(joined_dat$date)
  joined_dat$fyear <- as.factor(joined_dat$year)
  joined_dat$yday <- yday(joined_dat$date)
  # assign folds with blockCV -- first convert to sf
  #pa_data <- sf::st_as_sf(joined_dat, coords = c("X", "Y"), crs = 32610)
  # est_range <- blockCV::cv_spatial_autocor(x = pa_data, 
  #                             column = "present")
  # sb <- try(blockCV::cv_spatial(
  #   x = pa_data,
  #   column = "present",
  #   size = df$range[i], # cv_spatial() uses meters
  #   k = df$folds[i], # k has to be smaller than number of blocks
  #   selection = as.character(df$selection[i])
  # ), silent=TRUE)
  # 
  #if(class(sb)!="try-error") {
  
  breaks <- seq(min(joined_dat$Y), max(joined_dat$Y), by = df$bin_width[i])
  bins <- cut(joined_dat$Y, breaks = breaks, include.lowest = TRUE, labels=FALSE)
  bins[which(is.na(bins))] <- max(bins,na.rm=T) + 1
  bins_to_folds <- data.frame(bin = 1:max(bins), fold = rep(1:10, 1000)[1:max(bins)])
  
  joined_dat$bin <- bins
  joined_dat <- dplyr::left_join(joined_dat, bins_to_folds)
  
    coordinates(joined_dat) <- c("X", "Y")
    
    # create boundary, same for all meshes
    boundary <- inla.nonconvex.hull(coordinates(joined_dat),
        convex = -0.05 # Negative values are interpreted as fractions of the approximate initial set
    )
    
    # create mesh
    inla_mesh <- inla.mesh.2d(
      loc = coordinates(joined_dat),
      boundary = boundary,
      offset = c(-0.05, -0.05),
      max.n = 5000,
      max.n.strict = 5000,
      cutoff = df$cutoff[i],
      max.edge = c(100, 500)
    )
    joined_dat <- as.data.frame(joined_dat)
    mesh <- make_mesh(joined_dat, c("X", "Y"), mesh = inla_mesh)
    df$n[i] <- mesh$mesh$n
    
    # fit sdmTMB model
    # use delta-model to look at separate likelihoods for pres-abs and pos
    # don't include depth as predictor, as it's confounded with spatial field
    #joined_dat$zday <- scale(joined_dat$yday)
    fit_present <- try(sdmTMB_cv(present ~ -1 + as.factor(year),
                     data = joined_dat,
                     mesh = mesh,
                     time = "year",
                     spatial = "on",
                     spatiotemporal = "iid",
                     family = binomial(),
                     k_folds = df$folds[i],
                     fold_ids = joined_dat$fold), silent=TRUE)
    

    # fit_total <- sdmTMB_cv(cpue_kg_km2 ~ -1 + as.factor(year) + zday + I(zday^2),
    #                          data = as.data.frame(joined_dat),
    #                          mesh = mesh,
    #                          time = "year",
    #                          spatial = "on",
    #                          spatiotemporal = "iid",
    #                          family = delta_gamma(),
    #                          k_folds = df$folds[i],
    #                          fold_ids = "fold")
    # 
    if(class(fit_present) != "try-error") {
      df$present_dens_ll[i] <- mean(fit_present$fold_loglik) # total log density
      df$present_sd_ll[i] <- sd(fit_present$fold_loglik)
      df$present_converged[i] <- fit_present$converged
      
      tidy_pars <- lapply(fit_present$models, tidy, effects = "ran_pars")
      ranges <- unlist(lapply(lapply(tidy_pars, getElement, 2), getElement, 1))
      df$range[i] <- mean(ranges)
      
      sigma_O <- unlist(lapply(lapply(tidy_pars, getElement, 2), getElement, 2))
      df$sigma_O[i] <- mean(sigma_O)
      
      sigma_E <- unlist(lapply(lapply(tidy_pars, getElement, 2), getElement, 3))
      df$sigma_E[i] <- mean(sigma_E)
      
      tidy_ran <- purrr::map_dfr(tidy_pars, .f = ~ .x, .id = "id")
      tidy_pars <- lapply(fit_present$models, tidy)
      tidy_fixef <- purrr::map_dfr(tidy_pars, .f = ~ .x, .id = "id")
      tidy_all <- rbind(tidy_ran, tidy_fixef)
      tidy_all$i <- i
      
      if(i == 1) {
        all_pars_to_save <- tidy_all
      } else {
        all_pars_to_save <- rbind(all_pars_to_save, tidy_all)
      }
    }
    # df$total_dens_ll[i] <- fit_total$sum_loglik # total log density
    # df$total_sd_ll[i] <- sd(fit_total$fold_loglik)
    # df$total_converged[i] <- fit_total$converged
  #}
  saveRDS(df, "output/02_binomial_blockCV.rds")
  saveRDS(all_pars_to_save, "output/all_pars_to_save_binomial_blockCV.rds")
  print(i)
}

library(ggplot2)
d <- readRDS("output/02_binomial_blockCV.rds")

# filter out the number of species that don't converge enough
d <- dplyr::filter(d, present_converged == TRUE) |>
  dplyr::group_by(species) |>
  dplyr::mutate(nobs = n()) |>
  dplyr::filter(nobs >= 20) |>
  dplyr::select(-nobs)
d$species <- as.factor(as.character(d$species))

# Bring in the random 
d_random <- readRDS("output/02_binomial_randomCV.rds")
d_random <- dplyr::filter(d_random, present_converged == TRUE) |>
  dplyr::filter(species %in% d$species) |>
  dplyr::mutate(bin_width = NA)

# Larger bins result in more widely spaced test regions, beyond the estimated
# range. These regions are no longer correlated with the training data and predictions
# become more uncertain / not as good
dplyr::filter(d) |>
  ggplot(aes(n, present_dens_ll, group = bin_width, col = bin_width)) + 
  geom_line() + 
  facet_wrap(~species, scale = "free_y")


d |>
ggplot(aes(n, range, group = bin_width, col = bin_width)) + 
  geom_line() + 
  facet_wrap(~species, scale = "free_y") + 
  xlab("Mesh vertices") + ylab("Estimated spatial range (km)") + 
  scale_color_viridis(option="magma", begin = 0.2, end = 0.8, name = "Strip width (km)") + 
  theme_bw() + 
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size=5),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) + 
  geom_point(data = d_random, aes(n, range), col="black", alpha=0.5)
ggsave("figures/groundfish_range_v_n.png", height = 5, width = 7)


library(scales)
d |> # divide by 12670 to get average ll per obs
  ggplot(aes(n, present_dens_ll / 12670, group = bin_width, col = bin_width)) + 
  geom_line() + 
  facet_wrap(~species, scale = "free_y") + 
  xlab("Mesh vertices") + ylab("Predicted log likelihood") + 
  scale_color_viridis(option="magma", begin = 0.2, end = 0.8, name = "Strip width (km)") + 
  theme_bw() + 
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size=5),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=8)) + 
  geom_point(data = d_random, aes(n, present_dens_ll / 12670), col="black", alpha=0.5)# + 
  #scale_y_continuous(labels = function(x) format(x, scientific=TRUE))
ggsave("figures/groundfish_loglik_v_n.png", height = 5, width = 7)





