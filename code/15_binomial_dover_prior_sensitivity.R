# remotes::install_github("inlabru-org/fmesher", ref = "stable")
# library(sf)
library(inlabru)
library(INLA)
library(dplyr)
# library(fmesher)
# https://www.maths.ed.ac.uk/~flindgre/2018/07/22/spatially-varying-mesh-quality/

haul <- readRDS("data/haul_cleaned.rds")
catch <- readRDS("data/catch_cleaned.rds")

# species are dover sole, sablefish, petrale sole, lingcod, darkblotched rockfish
# g = group_by(catch, common_name) %>% dplyr::summarise(p = sum(present)) %>% dplyr::arrange(-p)

# initial loop over the cutoff values
prior_sigma <- data.frame(
  "prior.sigma_cutoff" = c(1, 5, 10, rep(5, 3)),
  "prior.sigma_thresh" = c(rep(0.05, 3), 0.01, 0.1, 0.4),
  "prior.range_cutoff" = 20,
  "prior.range_thresh" = 0.05,
  "sensitivity" = c(rep("sigma_cutoff", 3), rep("sigma_p", 3))
)
prior_range <- data.frame(
  "prior.range_cutoff" = c(5, 20, 40, rep(20, 3)),
  "prior.range_thresh" = c(rep(0.05, 3), 0.01, 0.1, 0.4),
  "prior.sigma_cutoff" = 5,
  "prior.sigma_thresh" = 0.05,
  "sensitivity" = c(rep("range_cutoff", 3), rep("range_p", 3))
)
priors <- rbind(prior_sigma, prior_range)
priors$prior_id <- seq(1, nrow(priors))

df <- expand.grid(
  "cutoff" = seq(3, 120, by = 6),
  "blocks" = c(20),
  "species" = c(
    "Dover sole"
  ),
  "n" = NA,
  "dens_ll" = NA,
  "dens_ll_train" = NA,
  prior_id = seq(1, max(priors$prior_id))
)
df <- dplyr::left_join(df, priors)


set.seed(2021)

for (i in 1:nrow(df)) {
  catch_sub <- dplyr::filter(catch, common_name == df$species[i],
                             year==2018)

  # Join catch and haul data
  haul_new <- haul %>%
    dplyr::filter(trawl_id %in% catch_sub$trawl_id) %>% 
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
  haul_new <- dplyr::left_join(as.data.frame(haul_new), block_fold)
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

  # use PC prior for matern model
  matern <-
    inla.spde2.pcmatern(mesh,
      prior.sigma = c(df$prior.sigma_cutoff[i], df$prior.sigma_thresh[i]),
      prior.range = c(df$prior.range_cutoff[i], df$prior.range_thresh[i])
    )
  # components is equivalent to formula
  # components <- present ~ Intercept + field(
  #   main = coordinates,
  #   model = matern
  # )
  components <- present ~ Intercept(1) + log_depth_scaled + log_depth_scaled2 + field(
    main = coordinates,
    model = matern
  )

  # do cross validation here
  haul_new$pred <- NA
  haul_new$predtrain <- NA
  for (k in 1:max(haul_new$fold)) {
    fit_train <- try(bru(components,
      haul_new[which(haul_new$fold != k), ],
      family = "binomial"
    ), silent = TRUE)
    test_indx <- which(haul_new$fold == k)
    # if model didn't have problems
    if (class(fit_train)[1] == "bru") {
      if (fit_train$ok == TRUE) {
        pred_test <- predict(fit_train,
          data = haul_new[test_indx, , drop = FALSE],
          formula = ~ Intercept + log_depth_scaled + log_depth_scaled2 + field
          #formula = ~ Intercept + field
        )
        haul_new$pred[test_indx] <- pred_test$mean

        if (k == 1) {
          pred_train <- predict(fit_train,
            data = haul_new[which(haul_new$fold != k), ],
            formula = ~ Intercept + log_depth_scaled + log_depth_scaled2 + field
            #formula = ~ Intercept + field
          )
          haul_new$predtrain[which(haul_new$fold != k)] <- pred_train$mean
        }
      } else {
        # model had issues
        haul_new$pred[test_indx] <- NA

        if (k == 1) {
          haul_new$predtrain[which(haul_new$fold != k)] <- NA
        }
      }
    }
  }
  # calculate total log density
  df$dens_ll[i] <-
    sum(dbinom(
      haul_new$present,
      size = 1,
      prob = plogis(haul_new$pred),
      log = TRUE
    ))

  # calculate total log density - training data
  df$dens_ll_train[i] <-
    sum(dbinom(
      haul_new$present,
      size = 1,
      prob = plogis(haul_new$predtrain),
      log = TRUE
    ), na.rm = T)
  saveRDS(df, "output/15_binom_dens_4species.rds")
}
