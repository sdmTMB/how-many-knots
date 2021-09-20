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
  "cutoff" = seq(3, 120, by = 6),
  "blocks" = c(10),
  "species" = c(
    "Dover sole", "petrale sole", "darkblotched rockfish",
    "lingcod"
  ),
  "n" = NA,
  "dens_ll" = NA
)

set.seed(2021)

for (i in 1:nrow(df)) {
  catch_sub <- dplyr::filter(catch, common_name == df$species[i])

  # Join catch and haul data
  haul_new <- haul %>%
    left_join(catch_sub, by = "trawl_id") %>%
    select(trawl_id, X, Y,
      latitude = latitude_dd.x,
      longitude = longitude_dd.x,
      year = year,
      log_depth_scaled,
      log_depth_scaled2,
      cpue_kg_km2
    )
  # Set NA CPUEs to 0
  haul_new$cpue_kg_km2[which(is.na(haul_new$cpue_kg_km2))] <- 0

  # convert coordinates to km
  haul_new$X <- haul_new$X / 1000
  haul_new$Y <- haul_new$Y / 1000
  # create occurrence field
  haul_new$present <- ifelse(haul_new$cpue_kg_km2 > 0, 1, 0)
  coordinates(haul_new) <- c("X", "Y")

  # create boundary, same for all meshes
  boundary <- inla.nonconvex.hull(coordinates(haul_new),
    convex = -0.05
  )

  n_folds <- 10
  n_blocks <- df$blocks[i]
  # first assign blocks
  haul_new$fold <- NULL
  haul_new$block <- 1
  for (jj in 1:n_blocks) {
    haul_new$block[which(haul_new$latitude < quantile(haul_new$latitude, 1 - jj * (1 / n_blocks)))] <- jj + 1
  }
  # now assign folds
  block_fold <- data.frame(
    "block" = 1:n_blocks,
    "fold" = rep(1:n_folds, n_blocks / n_folds)
  )
  haul_new <- dplyr::left_join(as.data.frame(haul_new), block_fold, by = "block")
  # return to SpatialPointsDataFrame
  coordinates(haul_new) <- c("X", "Y")

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
  # mesh_sdmTMB <- sdmTMB::make_mesh(data = haul_df, xy_cols = c("X", "Y"), cutoff = df$cutoff[i])
  # mesh_sdmTMB$mesh$n
  #
  # mesh_sdmTMB <- sdmTMB::make_mesh(data = haul_df, xy_cols = c("X", "Y"), cutoff = 10)
  # m0 <- sdmTMB(present ~ 1, data = haul_df, spde = mesh_sdmTMB, family = binomial(link = "logit"), silent = FALSE)
  #
  # mesh_sdmTMB <- sdmTMB::make_mesh(data = haul_df, xy_cols = c("X", "Y"), cutoff = 10)
  # m1 <- sdmTMB(present ~ 1, data = haul_df, spde = mesh_sdmTMB, family = binomial(link = "logit"), silent = FALSE, priors = sdmTMBpriors(
  #   matern_s = pc_matern(
  #     range_gt = 5, range_prob = 0.05,
  #     sigma_lt = 20, sigma_prob = 0.05
  #   )))
  # m1
  #
  # mesh_sdmTMB <- sdmTMB::make_mesh(data = haul_df, xy_cols = c("X", "Y"), cutoff = 100)
  # m2 <- sdmTMB(present ~ 1, data = haul_df, spde = mesh_sdmTMB, family = binomial(link = "logit"), silent = FALSE, priors = sdmTMBpriors(
  #   matern_s = pc_matern(
  #     range_gt = 5, range_prob = 0.05,
  #     sigma_lt = 20, sigma_prob = 0.05
  #   )))
  # m1
  # m2
  #
  fit <- sdmTMB::sdmTMB_cv(
    formula = present ~ 1,
    data = haul_df,
    spde = mesh_sdmTMB,
    parallel = TRUE,
    fold_ids = haul_df$fold,
    family = binomial(link = "logit"),
    priors = sdmTMBpriors(
      matern_s = pc_matern(
        range_gt = 5, range_prob = 0.05,
        sigma_lt = 20, sigma_prob = 0.05
      )
    )
  )
  df$dens_ll[i] <- fit$sum_loglik
}
saveRDS(df, "output/08_binom_dens_4species_TMB.rds")

plan(sequential)
