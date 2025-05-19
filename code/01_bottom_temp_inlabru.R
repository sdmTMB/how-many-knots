library(sf)
library(inlabru)
library(INLA)
library(dplyr)
library(ggplot2)
library(future)
library(sdmTMB)
plan(multisession)
library(sp)

set.seed(2021)

haul <- readRDS("data/haul_cleaned.rds")
haul$date_formatted <- as.Date(as.numeric(haul$date_formatted), origin = "1970-01-01")
haul$year <- lubridate::year(haul$date_formatted)
haul$month <- lubridate::month(haul$date_formatted)
haul$yday <- lubridate::yday(haul$date_formatted)
haul$zday <- as.numeric(scale(haul$yday))
haul$log_depth_scaled <- as.numeric(haul$log_depth_scaled)
haul$log_depth_scaled2 <- as.numeric(haul$log_depth_scaled2)
haul_new <- dplyr::filter(haul, year == 2018, !is.na(temperature_at_gear_c_der))
haul_new <- dplyr::filter(haul_new, !is.na(log_depth_scaled))
haul_new$fold_id <- rep(1:10, length.out = nrow(haul_new))
haul_new$pred_train <- NA
haul <- haul_new

ns <- dplyr::group_by(haul_new, fold_id) |>
  dplyr::summarise(n = n(), n_other = nrow(haul_new) - n) |>
  dplyr::ungroup() |>
  dplyr::summarise(mean_n = mean(n), mean_other = mean(n_other))

# haul$fold <- sample(c(1,2), size=nrow(haul), replace=T, prob = c(0.1,0.9))
# initial loop over the cutoff values

sp::coordinates(haul) <- c("X", "Y")
boundary <- INLA::inla.nonconvex.hull(sp::coordinates(haul),
  convex = -0.05
)

run_cv <- function(cutoff) {
  df <- data.frame(cutoff = cutoff)
  mesh <- INLA::inla.mesh.2d(
    loc = sp::coordinates(haul),
    boundary = boundary,
    offset = c(-0.05, -0.05),
    max.n = 5000,
    max.n.strict = 5000,
    cutoff = cutoff,
    max.edge = c(100, 500)
  )
  df$n <- mesh$n

  # Create SPDE model with PC priors
  matern <- INLA::inla.spde2.pcmatern(
    mesh,
    prior.sigma = c(10, 0.05),
    prior.range = c(500, 0.05)
  )

  components <- temperature_at_gear_c_der ~ log_depth_scaled + I(log_depth_scaled^2) +
    field(main = coordinates, model = matern)


  df$ll_test <- 0
  df$ll_train <- 0

  for (ii in 1:10) {
    print(ii)
    fit_train <- tryCatch(
      {
        bru(
          components,
          data = haul[haul$fold_id != ii, , drop = FALSE],
          family = "gaussian"
        )
      },
      error = function(e) {
        message("    -> bru() failed: ", e$message)
        return(NULL)
      }
    )

    if (is.null(fit_train)) {
      fold_failed <- TRUE
      break
    }

    pred_train <- tryCatch(
      {
        predict(fit_train,
          newdata = haul[haul$fold_id != ii, , drop = FALSE],
          ~ log_depth_scaled + I(log_depth_scaled^2) + field
        )
      },
      error = function(e) {
        message("    -> pred_train failed: ", e$message)
        return(NULL)
      }
    )

    pred_test <- tryCatch(
      {
        predict(fit_train,
          newdata = haul[haul$fold_id == ii, , drop = FALSE],
          ~ log_depth_scaled + I(log_depth_scaled^2) + field
        )
      },
      error = function(e) {
        message("    -> pred_test failed: ", e$message)
        return(NULL)
      }
    )

    tau <- tryCatch(
      {
        fit_train$summary.hyperpar$mean[1]
      },
      error = function(e) {
        message("    -> failed to extract tau: ", e$message)
        return(NA)
      }
    )

    # If any part failed, mark the fold as failed
    if (is.null(pred_train) || is.null(pred_test) || is.na(tau)) {
      fold_failed <- TRUE
      break
    }

    # Otherwise accumulate log-likelihoods
    df$ll_train <- df$ll_train +
      sum(dnorm(
        haul$temperature_at_gear_c_der[haul$fold_id != ii],
        mean = pred_train$mean,
        sd = sqrt(1 / tau),
        log = TRUE
      ))

    df$ll_test <- df$ll_test +
      sum(dnorm(
        haul$temperature_at_gear_c_der[haul$fold_id == ii],
        mean = pred_test$mean,
        sd = sqrt(1 / tau),
        log = TRUE
      ))
  }
  df$tau <- tau
  df
}

# out <- run_cv(cutoff = 25)

torun <- data.frame(cutoff = exp(seq(log(8), log(500), length.out = 24)))
nrow(torun)
torun$cutoff

plan(multisession, workers = 8)
out <- furrr::future_pmap(torun, run_cv)
plan(sequential)

out <- out |> dplyr::bind_rows()

out |>
  ggplot(aes(n, ll_test)) +
  geom_point() +
  xlab("Mesh vertices") +
  ylab("EDF") +
  ggtitle("Temperature example")
# dplyr::filter(df, converged==TRUE, n < nrow(haul_new)) |>
#   ggplot(aes(n, rmse_train)) + geom_point()
# dplyr::filter(df, converged==TRUE, n < nrow(haul_new)) |>
#   ggplot(aes(n, rmse_test)) + geom_line()

# dplyr::filter(df, converged == TRUE, n < nrow(haul_new)) |>
#   ggplot(aes(n, dens_ll_train)) +
#   geom_point()
# dplyr::filter(df, converged == TRUE, n < nrow(haul_new)) |>
#   ggplot(aes(n, dens_ll_test)) +
#   geom_point() +
#   geom_line()
