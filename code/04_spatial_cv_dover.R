# remotes::install_github("inlabru-org/fmesher", ref = "stable")
library(sf)
library(inlabru)
library(INLA)
library(dplyr)
# library(fmesher)
# https://www.maths.ed.ac.uk/~flindgre/2018/07/22/spatially-varying-mesh-quality/
library(blockCV)
library(future)
plan(multisession)

dover <- readRDS("data/doversole_cleaned.rds")
# convert coordinates to km
dover$X <- dover$X / 1000
dover$Y <- dover$Y / 1000
# create occurrence field
dover$present <- ifelse(dover$cpue_kg_km2 > 0, 1, 0)

set.seed(2021)

coordinates(dover) <- c("X", "Y")

# initial loop over the cutoff values
df <- expand.grid(
  "cutoff" = seq(5, 135, by = 10),
  "range" = seq(25, 150, by = 25),
  "folds" = c(10, 40),
  "n" = NA,
  "dens_ll" = NA
)
# create boundary, same for all meshes
boundary <- inla.nonconvex.hull(coordinates(dover),
  convex = -0.05
)

for (i in 1:nrow(df)) {
  # n_blocks <- df$blocks[i]
  # # first assign blocks
  # dover$fold <- NULL
  # dover$block <- 1
  # for (jj in 1:n_blocks) {
  #   dover$block[which(dover$latitude < quantile(dover$latitude, 1 - jj * (1 / n_blocks)))] = jj + 1
  # }
  # # now assign folds
  # block_fold <- data.frame("block" = 1:n_blocks,
  #                          "fold" = rep(1:n_folds, n_blocks / n_folds))
  # dover <- dplyr::left_join(as.data.frame(dover), block_fold)
  # # return to SpatialPointsDataFrame
  # coordinates(dover) <- c("X", "Y")
  #
  pa_data <- sf::st_as_sf(dover, coords = c("longitude", "latitude"))
  sb <- spatialBlock(
    speciesData = pa_data,
    species = "present",
    theRange = df$range[i],
    k = df$folds[i],
    selection = "systematic",
    showBlocks = FALSE
  )
  dover$fold <- sb$foldID
  df$blocks[i] <- nrow(as.data.frame(sb$blocks))

  # create mesh
  mesh <- inla.mesh.2d(
    loc = coordinates(dover),
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
      prior.sigma = c(10, 0.01),
      prior.range = c(1, 0.01)
    )
  # components is equivalent to formula
  components <- present ~ Intercept + field(
    main = coordinates,
    model = matern
  )

  # do cross validation here
  dover$pred <- NA
  for (k in 1:max(dover$fold)) {
    fit_train <- try(bru(components,
      dover[which(dover$fold != k), ],
      family = "binomial"
    ), silent = TRUE)
    test_indx <- which(dover$fold == k)
    # if model didn't have problems
    if (class(fit_train)[1] == "bru") {
      if (fit_train$ok) {
        
        pred_test <- predict(fit_train, 
          data = dover[test_indx, , drop = FALSE], 
          formula = ~ Intercept + field
        )
        dover$pred[test_indx] <- pred_test$mean
      } else {
        # model had issues
        dover$pred[test_indx] <- NA
      }
    }
  }
  # calculate total log density
  df$dens_ll[i] <-
    sum(dbinom(
      dover$present,
      size = 1,
      prob = plogis(dover$pred),
      log = TRUE
    ))
  saveRDS(df, "output/04_binomial_cv_df.rds")
}
