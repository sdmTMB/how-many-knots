library(geoR)
library(ggplot2)
library(sdmTMB)
library(purrr)
library(future)
library(dplyr)
library(tidyr)
plan(multisession)

kappa <- sqrt(8) / range

create_geoR_range <- function(range) {
  range / sqrt(8)
}
create_geoR_range(0.1)

set.seed(1)
sim1 <- grf(2000, cov.pars = c(1, create_geoR_range(0.2)), cov.model = "matern", kappa = 1)
dat <- data.frame(x = sim1$coords[, 1], y = sim1$coords[, 2], z = sim1$data, obs = sim1$data + rnorm(2000, 0, 1))

ggplot(dat, aes(x, y, colour = z)) +
  geom_point() +
  scale_colour_viridis_c()

ggplot(dat, aes(x, y, colour = obs)) +
  geom_point() +
  scale_colour_viridis_c()

mesh <- make_mesh(dat, c("x", "y"), cutoff = 0.01)
mesh$mesh$n
plot(mesh)

fit <- sdmTMB(obs ~ 1, data = dat, mesh = mesh)
summary(fit)
sanity(fit)
tidy(fit, "ran_pars")

cutoffs <- seq(0.005, 0.15, length.out = 20)

set.seed(1)
dat$fold_id <- sample(1:10, nrow(dat), replace = TRUE)

ret <- furrr::future_map_dfr(cutoffs, \(x) {
  print(x)
  mesh <- make_mesh(dat, c("x", "y"), cutoff = x)
  fit_cv <- sdmTMB_cv(
    obs ~ 1,
    data = dat, mesh = mesh,
    k_folds = 10, parallel = FALSE,
    fold_ids = "fold_id"
  )
  tidy_out <- map_dfr(fit_cv$models, \(m) tidy(m, "ran_pars"))
  grab_coef <- function(.term) {
    tidy_out |>
      filter(term == .term) |>
      pull(estimate) |>
      mean()
  }
  phi <- grab_coef("phi")
  sigma_O <- grab_coef("sigma_O")
  range <- grab_coef("range")
  rmse_true <- sqrt(mean((fit_cv$data$z - fit_cv$data$cv_predicted)^2))
  data.frame(
    leftout_loglik = fit_cv$sum_loglik,
    cutoff = x,
    n = mesh$mesh$n,
    phi = phi,
    sigma_O = sigma_O,
    range = range,
    rmse_true = rmse_true
  )
})

true <- data.frame(
  name = c("phi", "sigma_O", "range"),
  value = c(1, 1, 0.2)
)

pivot_longer(ret, cols = -n) |>
  ggplot(aes(n, value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~name, scales = "free_y") +
  geom_hline(data = true, mapping = aes(yintercept = value)) +
  ggsidekick::theme_sleek()
