library(INLA)
library(dplyr)
library(sdmTMB)
is_rstudio <- !is.na(Sys.getenv("RSTUDIO", unset = NA))
is_unix <- .Platform$OS.type == "unix"
library(future)
if (is_rstudio || !is_unix) plan(multisession) else plan(multicore)

haul <- readRDS("data/haul_cleaned.rds")
catch <- readRDS("data/catch_cleaned.rds")

# initial loop over the cutoff values
df <- expand.grid(
  "cutoff" = seq(10, 120, by = 5),
  "blocks" = c(20),
  "species" = unique(catch$common_name),
  "n" = NA,
  "dens_ll" = NA
)

set.seed(2021)

for (i in 1:nrow(df)) {
  catch_sub <- dplyr::filter(catch, common_name == df$species[i])

  # Join catch and haul data
  haul_new <- haul %>%
    left_join(catch_sub, by = "trawl_id") %>%
    dplyr::select(trawl_id, X, Y,
      latitude = latitude_dd.x,
      longitude = longitude_dd.x,
      year = year,
      log_depth_scaled,
      log_depth_scaled2,
      cpue_kg_km2
    ) %>%
    dplyr::filter(!is.na(year))

  # convert coordinates to km
  haul_new$X <- haul_new$X / 1000
  haul_new$Y <- haul_new$Y / 1000
  # create occurrence field
  #haul_new$present <- ifelse(haul_new$cpue_kg_km2 > 0, 1, 0)
  coordinates(haul_new) <- c("X", "Y")

  # create boundary, same for all meshes
  boundary <- inla.nonconvex.hull(coordinates(haul_new),
    convex = -0.05
  )
  
  #n_folds <- 10
  n_blocks <- df$blocks[i]
  
  haul_new$fold <- NULL
  haul_new$block <- 1
  for (jj in 1:n_blocks) {
    haul_new$block[which(haul_new$latitude < quantile(haul_new$latitude, 1 - jj * (1 / n_blocks)))] <- jj + 1
  }
  # now assign folds
  #block_fold <- data.frame(
  #  "block" = 1:n_blocks,
  #  "fold" = rep(1:n_folds, n_blocks / n_folds)
  #)
  #haul_new <- dplyr::left_join(as.data.frame(haul_new), block_fold, by = "block")
  # return to SpatialPointsDataFrame
  #coordinates(haul_new) <- c("X", "Y")

  # randomly hold out two blocks for the test case -- blocks 5 and 15 are test
  haul_new$fold <- ifelse(haul_new$block %in% c(5,15), 2, 1)
  
  # create mesh
  mesh <- inla.mesh.2d(
    loc = coordinates(haul_new),
    boundary = boundary,
    offset = c(-0.05, -0.05),
    max.n = 5000,
    max.n.strict = 5000,
    cutoff = df$cutoff[i],
    max.edge = c(100, 500)
  )
  df$n[i] <- mesh$n

  haul_df <- as.data.frame(haul_new)
  mesh_sdmTMB <- sdmTMB::make_mesh(data = haul_df, xy_cols = c("X", "Y"), mesh = mesh)
  
  fit <- sdmTMB::sdmTMB_cv(
    formula = cpue_kg_km2 ~ log_depth_scaled + log_depth_scaled2 + as.factor(year),
    data = haul_df,
    time="year",
    mesh = mesh_sdmTMB,
    parallel = TRUE,
    spatial = "on",
    spatiotemporal = "iid",
    fold_ids = haul_df$fold,
    family = tweedie(link = "log"),
    priors = sdmTMBpriors(
      matern_s = pc_matern(
        range_gt = 5, range_prob = 0.05,
        sigma_lt = 20, sigma_prob = 0.05
      )
    )
  )
  
  df$dens_ll[i] <- fit$fold_loglik[2]
  saveRDS(df, "output/09_TMB_tweedie-index.rds")
}

plan(sequential)

library(sdmTMB)
library(dplyr)
library(ggplot2)
library(viridis)

df = readRDS("output/09_TMB_tweedie-index.rds")
# filter out POP because of low occurrence
df = dplyr::filter(df, n<1000, species!="Pacific ocean perch")
g = ggplot(df, 
       aes(cutoff, dens_ll)) + 
  geom_line() +
  facet_wrap(~species, scale="free_y") + 
  ylab("Predicive log density") + 
  xlab("Cutoff distance (km)")

# find best cutoff / species
best = dplyr::group_by(df, species) %>% 
  dplyr::filter(dens_ll == max(dens_ll))
highres = best
highres$cutoff <- 20
index_models <- rbind(best, highres)

grid = readRDS("data/wc_grid.rds")
#grid = dplyr::rename(grid, lon = X, lat = Y)
grid = dplyr::mutate(grid,
                     X = X*10, # scale to units of km
                     Y = Y*10,
                     depth_scaled = as.numeric(scale(-depth)),
                     depth_scaled2 = depth_scaled^2)

grid$cell = seq(1,nrow(grid))
pred_grid = expand.grid(cell = grid$cell, year = 2003:2019)
pred_grid = dplyr::left_join(pred_grid, grid)
pred_grid$year = pred_grid$year
#pred_grid$time = as.numeric(pred_grid$year) - floor(mean(unique(as.numeric(pred_grid$year))))

index = list() # list for indices
#for(i in 1:nrow(index_models)){
for(i in c(1,2,4,5)) {  
  catch_sub <- dplyr::filter(catch, common_name == index_models$species[i])
  
  # Join catch and haul data
  haul_new <- haul %>%
    left_join(catch_sub, by = "trawl_id") %>%
    dplyr::select(trawl_id, X, Y,
                  latitude = latitude_dd.x,
                  longitude = longitude_dd.x,
                  year = year,
                  log_depth_scaled,
                  log_depth_scaled2,
                  cpue_kg_km2
    ) %>%
    dplyr::filter(!is.na(year))
  
  # convert coordinates to km
  haul_new$X <- haul_new$X / 1000
  haul_new$Y <- haul_new$Y / 1000
  # create occurrence field
  #haul_new$present <- ifelse(haul_new$cpue_kg_km2 > 0, 1, 0)
  coordinates(haul_new) <- c("X", "Y")
  
  # create boundary, same for all meshes
  boundary <- inla.nonconvex.hull(coordinates(haul_new),
                                  convex = -0.05
  )
  
  mesh <- inla.mesh.2d(
    loc = coordinates(haul_new),
    boundary = boundary,
    offset = c(-0.05, -0.05),
    max.n = 5000,
    max.n.strict = 5000,
    cutoff = index_models$cutoff[i],
    max.edge = c(100, 500)
  )
  #df$n[i] <- mesh$n
  
  haul_df <- as.data.frame(haul_new)
  mesh_sdmTMB <- sdmTMB::make_mesh(data = haul_df, 
                                   xy_cols = c("X", "Y"), mesh = mesh)
  
  fit <- sdmTMB(
    formula = cpue_kg_km2 ~ log_depth_scaled + log_depth_scaled2 + as.factor(year),
    data = haul_df,
    time="year",
    mesh = mesh_sdmTMB,
    spatial = "on",
    spatiotemporal = "iid",
    family = tweedie(link = "log"),
    priors = sdmTMBpriors(
      matern_s = pc_matern(
        range_gt = 5, range_prob = 0.05,
        sigma_lt = 20, sigma_prob = 0.05
      )
    )
  )
  
  predictions <- predict(fit, newdata = pred_grid, sims = 500)
  #mean_null <- apply(null_predictions, 1, mean)
  #sd_null <- apply(null_predictions, 1, sd)
  #null_predictions_summ[[i]] <- cbind(mean_null, sd_null)
  index[[i]] <- get_index_sims(predictions)

}



