library(INLA)
library(dplyr)
library(sdmTMB)
is_rstudio <- !is.na(Sys.getenv("RSTUDIO", unset = NA))
is_unix <- .Platform$OS.type == "unix"
library(future)
if (is_rstudio || !is_unix) plan(multisession) else plan(multicore)

haul <- readRDS("data/haul_all_cleaned.rds")
catch <- readRDS("data/catch_all_cleaned.rds")
catch$year <- as.numeric(substr(catch$date,1,4))

# initial loop over the cutoff values
df <- expand.grid(
  "cutoff" = seq(3, 120, by = 12),
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
    formula = cpue_kg_km2 ~ log_depth_scaled + log_depth_scaled2,
    data = haul_df,
    time="year",
    mesh = mesh_sdmTMB,
    parallel = TRUE,
    spatial = "on",
    spatiotemporal = "off",
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
