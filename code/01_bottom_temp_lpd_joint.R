# Simple version of 01_bottom_temp_inlabru_mgcv.R:
# out-of-sample LPD across 6 mesh complexities with sdmTMB_cv()
library(dplyr)
library(ggplot2)
library(sdmTMB)
library(future)
library(mgcv)

haul <- readRDS("data/haul_cleaned.rds")
haul$date_formatted <- as.Date(as.numeric(haul$date_formatted), origin = "1970-01-01")
haul$year <- lubridate::year(haul$date_formatted)
haul$yday <- lubridate::yday(haul$date_formatted)
haul$zday <- as.numeric(scale(haul$yday))
haul$log_depth_scaled <- as.numeric(haul$log_depth_scaled)
haul <- haul |>
  filter(year == 2018, !is.na(temperature_at_gear_c_der)) |>
  filter(!is.na(log_depth_scaled), !is.na(zday))

set.seed(123)
haul$fold_id <- sample(1:10, size = nrow(haul), replace = TRUE)

boundary <- INLA::inla.nonconvex.hull(as.matrix(haul[, c("X", "Y")]), convex = -0.05)

make_inla_mesh <- function(cutoff) {
  INLA::inla.mesh.2d(
    loc = as.matrix(haul[, c("X", "Y")]),
    boundary = boundary,
    offset = c(-0.05, -0.05),
    max.n = 5000,
    max.n.strict = 5000,
    cutoff = cutoff,
    max.edge = c(100, 500)
  )
}

run_cv <- function(cutoff) {
  mesh <- make_inla_mesh(cutoff)
  sdmTMB_mesh <- make_mesh(haul, c("X", "Y"), mesh = mesh)
  fit <- sdmTMB_cv(
    temperature_at_gear_c_der ~ zday + I(zday^2),
    data = haul,
    mesh = sdmTMB_mesh,
    fold_ids = haul$fold_id,
    predictive = "joint",
    nsim = 500L,
    parallel = FALSE
  )
  data.frame(
    model = "sdmTMB",
    cutoff = cutoff,
    n = mesh$n,
    lpd_total = fit$sum_loglik,
    lpd = fit$sum_loglik / nrow(haul)
  )
}

# mgcv analogue of sdmTMB_cv(predictive = "random"): the spline coefficients
# (the penalized 'random effects') are drawn from their Bayesian posterior
# N(beta_hat, Vp) with the scale fixed at its estimate. Each held-out
# observation is scored as log(mean_s p(y_i | beta_s)) over nsim draws.
# As in sdmTMB, k = round(mesh vertices / 3) is the basis dimension.
run_cv_mgcv <- function(cutoff, nsim = 200L) {
  n_mesh <- make_inla_mesh(cutoff)$n
  k <- round(n_mesh / 3)
  ll_fold <- vapply(sort(unique(haul$fold_id)), function(ii) {
    train <- haul[haul$fold_id != ii, ]
    test <- haul[haul$fold_id == ii, ]
    fit <- tryCatch(
      mgcv::gam(temperature_at_gear_c_der ~ zday + I(zday^2) + s(X, Y, k = k), data = train),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NA_real_)
    Xp <- predict(fit, newdata = test, type = "lpmatrix")
    beta <- mgcv::rmvn(nsim, coef(fit), fit$Vp) # nsim x p
    mu <- Xp %*% t(beta) # n_test x nsim
    ll <- dnorm(test$temperature_at_gear_c_der, mean = mu, sd = sqrt(fit$sig2), log = TRUE)
    m <- apply(ll, 1, max)
    sum(m + log(rowSums(exp(ll - m))) - log(nsim))
  }, numeric(1))
  data.frame(
    model = "mgcv",
    cutoff = cutoff,
    n = k, # basis dimension
    lpd_total = sum(ll_fold),
    lpd = sum(ll_fold) / nrow(haul)
  )
}

cutoffs <- exp(seq(log(8), log(500), length.out = 6))

plan(multisession, workers = 6L)
opts <- furrr::furrr_options(seed = TRUE)
out <- bind_rows(
  furrr::future_map(cutoffs, run_cv, .options = opts),
  furrr::future_map(cutoffs, run_cv_mgcv, .options = opts)
)
plan(sequential)
print(out)

ggplot(out, aes(n, lpd)) +
  facet_wrap(~model, scales = "free_x") +
  geom_line() +
  geom_point() +
  xlab("Mesh vertices (sdmTMB)\nor basis dimension k (mgcv)") +
  ylab("Out-of-sample average\nlog predictive density")
