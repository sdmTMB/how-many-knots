# remotes::install_github("inlabru-org/fmesher", ref = "stable")
library(sf)
library(inlabru)
library(INLA)
library(dplyr)
library(ggplot2)
# library(fmesher)
# https://www.maths.ed.ac.uk/~flindgre/2018/07/22/spatially-varying-mesh-quality/
library(future)
plan(multisession)

haul <- readRDS("data/june_bottom_temp.rds")
haul[, c("X", "Y")] <- haul[, c("X", "Y")] / 1000
# species are dover sole, sablefish, petrale sole, lingcod, darkblotched rockfish
# g = group_by(catch, common_name) %>% dplyr::summarise(p = sum(present)) %>% dplyr::arrange(-p)

# initial loop over the cutoff values
df <- expand.grid(
  "cutoff" = seq(20, 160, by = 2.5),
  "blocks" = 10,
  "n" = NA,
  "dens_ll_train" = NA,
  "dens_ll_test" = NA,
  "tau" = NA
)

set.seed(2021)
meshes <- list()

for (i in 1:nrow(df)) {

  # Set NA CPUEs to 0
  # haul_new$cpue_kg_km2[which(is.na(haul_new$cpue_kg_km2))] <- 0
  haul_new <- dplyr::filter(
    haul,
    !is.na(temperature_at_gear_c_der),
  )
  # create occurrence field
  # haul_new$present <- ifelse(haul_new$cpue_kg_km2 > 0, 1, 0)
  coordinates(haul_new) <- c("X", "Y")

  # create boundary, same for all meshes
  boundary <- inla.nonconvex.hull(coordinates(haul_new),
    convex = -0.05
  )
  haul_new$fold <- 2 # training data
  haul_new$fold[seq(1, 175, by = 7)] <- 1 # test

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
  meshes[[i]] <- mesh

  df$n[i] <- mesh$n

  # use PC prior for matern model
  matern <-
    inla.spde2.pcmatern(mesh,
      prior.sigma = c(5, 0.05),
      prior.range = c(20, 0.05)
    )
  # components is equivalent to formula
  components <- temperature_at_gear_c_der ~ Intercept(1) + log_depth_scaled + log_depth_scaled2 +
    field(
      main = coordinates,
      model = matern
    )

  # do cross validation here
  haul_new$pred <- NA
  haul_new$se <- NA
  test_ll <- 0
  train_ll <- 0

  for (k in 1:1) {
    fit_train <- try(bru(components,
      haul_new[which(haul_new$fold != k), ],
      family = "gaussian"
    ), silent = TRUE)
    test_indx <- which(haul_new$fold == k)
    haul_new$pred[-test_indx] <- fit_train$summary.fitted.values$mean[1:150]
    haul_new$se[-test_indx] <- fit_train$summary.fitted.values$sd[1:150]
    # get precision parameter and calculate density for training data
    tau <- fit_train$summary.hyperpar$mean[1]
    train_ll[k] <- sum(dnorm(haul_new$temperature_at_gear_c_der[-test_indx],
      mean = haul_new$pred[-test_indx],
      sd = sqrt(1 / tau),
      log = TRUE
    ))

    # if model didn't have problems
    if (class(fit_train)[1] == "bru") {
      if (fit_train$ok) {
        pred_test <- predict(fit_train,
          data = haul_new[test_indx, , drop = FALSE],
          formula = ~ Intercept + log_depth_scaled + log_depth_scaled2 + field
        )
        haul_new$pred[test_indx] <- pred_test$mean
        haul_new$se[test_indx] <- pred_test$sd
      } else {
        # model had issues
        haul_new$pred[test_indx] <- NA
      }
    }

    # get precision parameter for the normal observations
    tau <- fit_train$summary.hyperpar$mean[1]
    # calculate density for test data
    test_ll[k] <- sum(dnorm(haul_new$temperature_at_gear_c_der[test_indx],
      mean = haul_new$pred[test_indx],
      sd = sqrt(1 / tau),
      log = TRUE
    ))
  }

  # save the data frame
  haul_new$n <- df$n[i]
  haul_new$cutoff <- df$cutoff[i]
  if (i == 1) {
    all_df <- haul_new
  } else {
    all_df <- rbind(all_df, haul_new)
  }

  # calculate total log density
  df$tau[i] <- tau
  df$dens_ll_test[i] <- sum(test_ll)
  df$dens_ll_train[i] <- sum(train_ll)
  saveRDS(df, file = "output/temp_model_df.rds")
  saveRDS(all_df, file = "output/temp_model_all_est.rds")
}

saveRDS(meshes, file = "output/temp_meshes.rds")
