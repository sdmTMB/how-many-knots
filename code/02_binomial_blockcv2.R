# remotes::install_github("inlabru-org/fmesher", ref = "stable")
library(sf)
# library(sp)
library(sdmTMB)
library(lubridate)
library(dplyr)
library(future)
library(viridis)
# plan(multisession)
haul <- readRDS("data/haul_cleaned.rds")
catch <- readRDS("data/catch_cleaned.rds")
haul$trawl_id <- as.numeric(haul$trawl_id)

run_cv <- function(cutoff, bin_width, species, bins_to_folds = bin_to_folds, do_full_fit = FALSE, parallel = TRUE) {
  catch_sub <- dplyr::filter(catch, common_name == species)

  # Join catch and haul data
  joined_dat <- haul |>
    dplyr::filter(trawl_id %in% catch$trawl_id) |>
    left_join(catch_sub[, c("date", "cpue_kg_km2", "depth_m", "trawl_id")],
      by = "trawl_id"
    ) |>
    dplyr::filter(!is.na(cpue_kg_km2)) |>
    dplyr::mutate(present = ifelse(cpue_kg_km2 > 0, 1, 0))
  joined_dat$date <- lubridate::ymd(as.POSIXct.Date(joined_dat$date))
  joined_dat$year <- year(joined_dat$date)
  joined_dat$fyear <- as.factor(joined_dat$year)
  joined_dat$yday <- yday(joined_dat$date)

  breaks <- seq(min(joined_dat$Y), max(joined_dat$Y), by = bin_width)
  bins <- cut(joined_dat$Y, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  bins[which(is.na(bins))] <- max(bins, na.rm = TRUE) + 1
  bins_to_folds <- data.frame(bin = 1:max(bins), fold = rep(1:10, 1000)[1:max(bins)])

  joined_dat$bin <- bins
  joined_dat <- dplyr::left_join(joined_dat, bins_to_folds)

  sp::coordinates(joined_dat) <- c("X", "Y")

  # create boundary, same for all meshes
  boundary <- fmesher::fm_nonconvex_hull_inla(sp::coordinates(joined_dat),
    convex = -0.05 # Negative values are interpreted as fractions of the approximate initial set
  )

  # create mesh
  inla_mesh <- fmesher::fm_mesh_2d_inla(
    loc = sp::coordinates(joined_dat),
    boundary = boundary,
    offset = c(-0.05, -0.05),
    max.n = 5000,
    max.n.strict = 5000,
    cutoff = cutoff,
    max.edge = c(100, 500)
  )
  joined_dat <- as.data.frame(joined_dat)
  mesh <- make_mesh(joined_dat, c("X", "Y"), mesh = inla_mesh)

  if (!do_full_fit) {
    cat("Cross validation model fit\n")
    fit_cv <- try(sdmTMB_cv(
      present ~ -1 + as.factor(year) + poly(log_depth_scaled, 2),
      data = joined_dat,
      mesh = mesh,
      time = "year",
      spatial = "on",
      spatiotemporal = "iid",
      anisotropy = TRUE,
      share_range = TRUE,
      family = binomial(),
      parallel = parallel,
      k_folds = length(unique(joined_dat$fold)),
      fold_ids = joined_dat$fold
    ), silent = TRUE)
  }

  if (do_full_fit) {
    cat("Full model fit\n")
    fit_full <- try(sdmTMB(
      present ~ -1 + as.factor(year) + poly(log_depth_scaled, 2),
      data = joined_dat,
      mesh = mesh,
      time = "year",
      spatial = "on",
      spatiotemporal = "iid",
      anisotropy = TRUE,
      share_range = TRUE,
      family = binomial()
    ), silent = TRUE)
    if (class(fit_full) != "try-error") {
      tidy_pars <- tidy(fit_full, "ran_pars")
      tidy_pars <- filter(tidy_pars, term != "range")
      aniso <- sdmTMB:::print_anisotropy(fit_full, return_dat = TRUE)
      tidy_pars <- bind_rows(
        tidy_pars,
        tibble(term = "range_a", estimate = as.numeric(aniso$sp$a)),
        tibble(term = "range_b", estimate = as.numeric(aniso$sp$b)),
        tibble(term = "angle", estimate = as.numeric(aniso$sp$degree))
      )
      tidy_pars$cutoff <- cutoff
      tidy_pars$bin_width <- bin_width
      tidy_pars$species <- species
      tidy_pars$n <- mesh$mesh$n
      return(tidy_pars)
    } else {
      return(
        data.frame(cutoff = cutoff, bin_width = bin_width, species = species, converged = FALSE, n = mesh$mesh$n)
      )
    }
  }

  if (class(fit_cv) != "try-error") {
    ll <- data.frame(cutoff = cutoff, bin_width = bin_width, species = species)
    ll$n <- mesh$mesh$n
    ll$converged <- fit_cv$converged
    ll$present_dens_ll <- mean(fit_cv$fold_loglik) # total log density
    ll$present_sd_ll <- sd(fit_cv$fold_loglik)
    ll$present_converged <- fit_cv$converged
    ll$cAIC <- as.numeric(cAIC(fit_full, what = "cAIC"))
    edf <- cAIC(fit_full, what = "EDF")
    ll$edf_epsilon <- edf[["epsilon_st"]]
    ll$edf_omega <- edf[["omega_s"]]
    return(ll)
  } else {
    return(
      data.frame(cutoff = cutoff, bin_width = bin_width, species = species, converged = FALSE, n = mesh$mesh$n)
    )
  }
}

if (FALSE) {
  df <- expand.grid(
    cutoff = 40,
    bin_width = 40,
    species = "sablefish"
  )
  tictoc::tic()
  plan(multisession)
  out <- purrr::pmap(df[1, ], run_cv, parallel = TRUE)
  tictoc::toc()
  out2 <- purrr::pmap(df[1, ], run_cv, do_full_fit = TRUE, parallel = FALSE)
}

df <- expand.grid(
  cutoff = round(exp(seq(log(12), log(175), length.out = 20))),
  bin_width = seq(10, 150, by = 20),
  species = c("sablefish", "arrowtooth", "lingcod", "widow rockfish")
)
nrow(df)
cores <- 40
nrow(df) * 560 / 60 / 60 / cores

set.seed(1234)
plan(multicore, workers = 40L)
out <- furrr::future_pmap(df, run_cv, parallel = FALSE)
out2 <- furrr::future_pmap(df, run_cv, do_full_fit = TRUE, parallel = FALSE)
plan(sequential)

save(out, "output/02-binomial-blockCV-2025-03-28.rds")
save(out2, "output/02-binomial-blockCV-pars-2025-03-28.rds")

## library(ggplot2)
## d <- readRDS("output/02_binomial_blockCV.rds")
##
## # filter out the number of species that don't converge enough
## d <- dplyr::filter(d, present_converged == TRUE) |>
##   dplyr::group_by(species) |>
##   dplyr::mutate(nobs = n()) |>
##   dplyr::filter(nobs >= 10) |>
##   dplyr::select(-nobs)
## d$species <- as.factor(as.character(d$species))
##
## # Bring in the random
## d_random <- readRDS("output/02_binomial_randomCV.rds")
## d_random <- dplyr::filter(d_random, present_converged == TRUE) |>
##   dplyr::filter(species %in% d$species) |>
##   dplyr::mutate(bin_width = NA)
##
## # Larger bins result in more widely spaced test regions, beyond the estimated
## # range. These regions are no longer correlated with the training data and predictions
## # become more uncertain / not as good
## dplyr::filter(d) |>
##   ggplot(aes(n, present_dens_ll, group = bin_width, col = bin_width)) +
##   geom_line() +
##   facet_wrap(~species, scale = "free_y")
##
##
## d |>
##   ggplot(aes(n, range, group = bin_width, col = bin_width)) +
##   geom_line() +
##   facet_wrap(~species, scale = "free_y") +
##   xlab("Mesh vertices") +
##   ylab("Estimated spatial range (km)") +
##   scale_color_viridis(option = "magma", begin = 0.2, end = 0.8, name = "Strip width (km)") +
##   theme_bw() +
##   theme(
##     strip.background = element_rect(fill = "white"),
##     strip.text = element_text(size = 5),
##     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
##   ) +
##   geom_point(data = d_random, aes(n, range), col = "black", alpha = 0.5)
## ggsave("figures/groundfish_range_v_n.png", height = 5, width = 7)
##
##
## library(scales)
## d |> # divide by 12670 to get average ll per obs
##   ggplot(aes(n, present_dens_ll / 12670, group = bin_width, col = bin_width)) +
##   geom_line() +
##   facet_wrap(~species, scale = "free_y") +
##   xlab("Mesh vertices") +
##   ylab("Predicted log likelihood") +
##   scale_color_viridis(option = "magma", begin = 0.2, end = 0.8, name = "Strip width (km)") +
##   theme_bw() +
##   theme(
##     strip.background = element_rect(fill = "white"),
##     strip.text = element_text(size = 5),
##     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
##     axis.text.y = element_text(size = 8)
##   ) +
##   geom_point(data = d_random, aes(n, present_dens_ll / 12670), col = "black", alpha = 0.5) # +
## # scale_y_continuous(labels = function(x) format(x, scientific=TRUE))
## ggsave("figures/groundfish_loglik_v_n.png", height = 5, width = 7)
