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

run_cv <- function(cutoff, folds = 1:10, run_inla = TRUE) {
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

  components <- temperature_at_gear_c_der ~ log_depth_scaled + log_depth_scaled2 + Intercept(1) +
    field(main = geometry, model = matern)

  df$ll_test <- 0
  df$ll_train <- 0
  df$ll_test_sdmTMB <- 0
  df$ll_train_sdmTMB <- 0
  df$ll_test_sdmTMB_noprior <- 0
  df$ll_train_sdmTMB_noprior <- 0

  for (ii in folds) {
    print(ii)
    if (run_inla) {
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
      # summary(fit_train)
      if (is.null(fit_train)) {
        fold_failed <- TRUE
        break
      }
    }

    this_dat <- haul[haul$fold_id != ii, , drop = FALSE]
    co <- sp::coordinates(this_dat)
    this_dat <- as.data.frame(this_dat)
    this_dat$X <- co[, 1]
    this_dat$Y <- co[, 2]
    this_mesh <- make_mesh(this_dat, c("X", "Y"), mesh = mesh)
    fit_train_sdmTMB <- tryCatch(sdmTMB(
      temperature_at_gear_c_der ~ log_depth_scaled + I(log_depth_scaled^2),
      data = this_dat,
      mesh = this_mesh,
      priors = sdmTMBpriors(pc_matern(range_gt = 500, sigma_lt = 10))
    ), error = function(e) {
      return(NULL)
    })

    fit_train_sdmTMB_noprior <- tryCatch(sdmTMB(
      temperature_at_gear_c_der ~ log_depth_scaled + I(log_depth_scaled^2),
      data = this_dat,
      mesh = this_mesh,
    ), error = function(e) {
      return(NULL)
    })

    if (is.null(fit_train_sdmTMB_noprior)) {
      break
    }

    if (run_inla) {
      pred_train <- tryCatch(
        {
          predict(fit_train,
            newdata = haul[haul$fold_id != ii, , drop = FALSE],
            ~ log_depth_scaled + log_depth_scaled2 + field + Intercept
          )
        },
        error = function(e) {
          message("    -> pred_train failed: ", e$message)
          return(NULL)
        }
      )
    }

    pred_train_sdmTMB <- predict(fit_train_sdmTMB, newdata = NULL)
    pred_train_sdmTMB_noprior <- predict(fit_train_sdmTMB_noprior, newdata = NULL)

    if (run_inla) {
      pred_test <- tryCatch(
        {
          predict(fit_train,
            newdata = haul[haul$fold_id == ii, , drop = FALSE],
            ~ log_depth_scaled + log_depth_scaled2 + field + Intercept
          )
        },
        error = function(e) {
          message("    -> pred_test failed: ", e$message)
          return(NULL)
        }
      )
    }

    nd <- haul[haul$fold_id == ii, , drop = FALSE]
    co <- sp::coordinates(nd)
    nd <- as.data.frame(nd)
    nd$X <- co[, 1]
    nd$Y <- co[, 2]
    pred_test_sdmTMB <- predict(fit_train_sdmTMB, newdata = nd)
    pred_test_sdmTMB_noprior <- predict(fit_train_sdmTMB_noprior, newdata = nd)

    if (FALSE) {
      plot(pred_test_sdmTMB$est, pred_test$mean)
      abline(0, 1)
      plot(pred_train_sdmTMB$est, pred_train$mean)
      abline(0, 1)
    }

    if (run_inla) {
      tau <- tryCatch(
        {
          fit_train$summary.hyperpar$mean[1]
        },
        error = function(e) {
          message("    -> failed to extract tau: ", e$message)
          return(NA)
        }
      )
    }

    phi <- exp(get_pars(fit_train_sdmTMB)$ln_phi)
    phi_noprior <- exp(get_pars(fit_train_sdmTMB_noprior)$ln_phi)

    if (run_inla) {
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

    df$ll_train_sdmTMB <- df$ll_train_sdmTMB +
      sum(dnorm(
        this_dat$temperature_at_gear_c_der,
        mean = pred_train_sdmTMB$est,
        sd = phi,
        log = TRUE
      ))

    df$ll_test_sdmTMB <- df$ll_test_sdmTMB +
      sum(dnorm(
        nd$temperature_at_gear_c_der,
        mean = pred_test_sdmTMB$est,
        sd = phi,
        log = TRUE
      ))

    df$ll_train_sdmTMB_noprior <- df$ll_train_sdmTMB_noprior +
      sum(dnorm(
        this_dat$temperature_at_gear_c_der,
        mean = pred_train_sdmTMB_noprior$est,
        sd = phi,
        log = TRUE
      ))

    df$ll_test_sdmTMB_noprior <- df$ll_test_sdmTMB_noprior +
      sum(dnorm(
        nd$temperature_at_gear_c_der,
        mean = pred_test_sdmTMB_noprior$est,
        sd = phi,
        log = TRUE
      ))
  }

  if (run_inla) {
    df$tau <- tau
  }
  df$phi <- phi
  df$phi_noprior <- phi_noprior
  df
}

# out <- run_cv(cutoff = 100, folds = 1:2)

# torun <- data.frame(cutoff = exp(seq(log(8), log(500), length.out = 30)))
torun <- data.frame(cutoff = exp(seq(log(8), log(500), length.out = 30)))
nrow(torun)
torun$cutoff

plan(multisession, workers = 10)
out <- furrr::future_pmap(torun, run_cv)
plan(sequential)

out2 <- out |> dplyr::bind_rows()
head(out2)

out3 <- tidyr::pivot_longer(select(out2, cutoff, n, ll_test, ll_train), cols = ll_test:ll_train)

g1 <- tidyr::pivot_longer(select(out2, cutoff, n, ll_test_sdmTMB, ll_train_sdmTMB), cols = ll_test_sdmTMB:ll_train_sdmTMB) |>
  filter(!(name == "ll_test_sdmTMB" & value < -500)) |>
  filter(!(name == "ll_train_sdmTMB" & value < -4000)) |>
  filter(value != 0) |>
  filter(value < -250) |>
  ggplot(aes(n, value, colour = name)) +
  facet_wrap(~name, scales = "free_y") +
  geom_point() +
  xlab("Mesh vertices") +
  ggtitle("sdmTMB: Temperature example")

g3 <- tidyr::pivot_longer(select(out2, cutoff, n, ll_test_sdmTMB_noprior, ll_train_sdmTMB_noprior), cols = ll_test_sdmTMB_noprior:ll_train_sdmTMB_noprior) |>
  filter(!(name == "ll_test_sdmTMB_noprior" & value < -500)) |>
  filter(!(name == "ll_train_sdmTMB" & value < -4000)) |>
  filter(value != 0) |>
  filter(value < -250) |>
  ggplot(aes(n, value, colour = name)) +
  facet_wrap(~name, scales = "free_y") +
  geom_point() +
  xlab("Mesh vertices") +
  ggtitle("sdmTMB no PC prior: Temperature example")

g2 <- out3 |>
  filter(value != 0) |>
  filter(value < -250) |>
  filter(value > -4000) |>
  filter(!(name == "ll_test" & value < -500)) |>
  # filter(!(name == "ll_test" & value != 0)) |>
  # filter(!(name == "ll_test" & value < -800)) |>
  ggplot(aes(n, value, colour = name)) +
  facet_wrap(~name, scales = "free_y") +
  geom_point() +
  xlab("Mesh vertices") +
  ggtitle("inlabru: Temperature example")
library(patchwork)
theme_set(theme_light())
g1 / g3/ g2

out2 |>
  ggplot(aes(n, ll_train)) +
  geom_point() +
  xlab("Mesh vertices") +
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
