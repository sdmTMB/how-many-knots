# remotes::install_github("inlabru-org/fmesher", ref = "stable")
library(sf)
library(inlabru)
library(INLA)
library(dplyr)
library(ggplot2)
library(future)
plan(multisession)

set.seed(2021)

haul <- readRDS("data/haul_cleaned.rds")
haul$date_formatted <- as.Date(as.numeric(haul$date_formatted), origin = "1970-01-01")
haul$year <- lubridate::year(haul$date_formatted)
haul$month <- lubridate::month(haul$date_formatted)
haul$yday <- lubridate::yday(haul$date_formatted)
haul$zday <- scale(haul$yday)
haul_new <- dplyr::filter(haul, year==2018,
                      !is.na(temperature_at_gear_c_der))
haul_new$fold_id <- rep(1:10, length.out = nrow(haul_new))
haul_new$pred_train <- NA

ns <- dplyr::group_by(haul_new, fold_id) |>
  dplyr::summarise(n = n(), n_other = nrow(haul_new)-n) |>
  dplyr::ungroup() |>
  dplyr::summarise(mean_n = mean(n), mean_other = mean(n_other))

#haul$fold <- sample(c(1,2), size=nrow(haul), replace=T, prob = c(0.1,0.9))
# initial loop over the cutoff values
df <- expand.grid(
  "cutoff" = exp(seq(log(8), log(500), length.out=200)),
  "n" = NA,
  "dens_ll_train" = NA,
  "dens_ll_test" = NA,
  "rmse_test" = NA,
  "rmse_train" = NA,
  "converged" = NA
)

for (i in 1:nrow(df)) {

  # This is all to use a custom mesh
  # haul_new$present <- ifelse(haul_new$cpue_kg_km2 > 0, 1, 0)
  #coordinates(haul_new) <- c("X", "Y")
  # create boundary, same for all meshes
  #boundary <- inla.nonconvex.hull(coordinates(haul_new),
  #  convex = -0.05
  #)
  # mesh <- inla.mesh.2d(
  #   loc = coordinates(haul_new),
  #   boundary = boundary,
  #   offset = c(-0.05, -0.05),
  #   max.n = 5000,
  #   max.n.strict = 5000,
  #   cutoff = df$cutoff[i],
  #   max.edge = c(100, 500)
  # )
  #haul_new <- as.data.frame(cbind(haul_new, coordinates(haul_new)))
  
  #haul_new$fold <- 2 # training data
  #haul_new$fold[seq(1, nrow(haul_new), by = 7)] <- 1 # test

  mesh_cv <- make_mesh(data = haul_new, 
                    xy_cols = c("X","Y"), 
                    cutoff = df$cutoff[i])
  df$n[i] <- mesh_cv$mesh$n

  # do cross validation here
  
  fit_cv <- sdmTMB::sdmTMB_cv(temperature_at_gear_c_der ~ 1 + zday + I(zday^2),
                            spatial="on",
                            fold_ids = haul_new$fold_id,
                            k_folds = max(haul_new$fold_id),
                            mesh = mesh_cv,
                            data = haul_new)
  df$converged[i] <- fit_cv$converged
  
  # also need in sample predictions
  ll_train <- NA
  rmse_train <- NA
  for(ii in 1:max(haul_new$fold_id)) {
    # predict to whole dataset, test and train
    pred <- predict(fit_cv$models[[ii]]) |>
      dplyr::filter(fold_id != ii)
    haul_train <- dplyr::filter(haul_new, fold_id != ii)
    # get phi
    phi <- as.numeric(tidy(fit_cv$models[[ii]],"ran_pars")[2,2])
    ll_train[ii] <- sum(dnorm(haul_train$temperature_at_gear_c_der, pred$est, phi, log=TRUE))
    rmse_train[ii] <- sqrt(mean( (haul_train$temperature_at_gear_c_der - pred$est)^2 ))
  }

  df$dens_ll_train[i] <- mean(ll_train)
  
  # calculate total log density
  df$dens_ll_test[i] <- mean(fit_cv$fold_loglik)
 
  rmse_test <- dplyr::group_by(fit_cv$data, fold_id) |>
    dplyr::summarise(rmse = sqrt(mean((cv_predicted - temperature_at_gear_c_der)^2)) )
  df$rmse_test[i] <- mean(rmse_test$rmse)
  df$rmse_train[i] <- mean(rmse_train)
  saveRDS(df, file = "output/temp_model_df.rds")
  
  # Save parameters
  fixef <- lapply(fit_cv$models, tidy)
  fixef <- purrr::map_dfr(fixef, ~ .x, .id = "index")
  fixef <- dplyr::group_by(fixef, term) |>
    dplyr::summarise(mean_estimate = mean(estimate), mean_se = mean(std.error))
  
  ranef <- lapply(fit_cv$models, tidy, "ran_pars")
  ranef <- purrr::map_dfr(ranef, ~ .x, .id = "index")
  ranef <- dplyr::group_by(ranef, term) |>
    dplyr::summarise(mean_estimate = mean(estimate), mean_se = mean(std.error))
  pars <- rbind(fixef, ranef)
  pars$cutoff <- df$cutoff[i]
  
  if(i==1) {
   all_pars <- pars
  } else {
   all_pars <- rbind(all_pars, pars)
  }
  saveRDS(all_pars, file = "output/temp_model_all_est.rds")
}

# dplyr::filter(df, converged==TRUE, n < nrow(haul_new)) |>
#   ggplot(aes(n, rmse_train)) + geom_point()
# dplyr::filter(df, converged==TRUE, n < nrow(haul_new)) |>
#   ggplot(aes(n, rmse_test)) + geom_line()

dplyr::filter(df, converged==TRUE, n < nrow(haul_new)) |>
  ggplot(aes(n, dens_ll_train)) + geom_point()
dplyr::filter(df, converged==TRUE, n < nrow(haul_new)) |>
  ggplot(aes(n, dens_ll_test)) + geom_point() + geom_line()

