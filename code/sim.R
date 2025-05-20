library(geoR)
library(ggplot2)
library(sdmTMB)
library(purrr)
library(future)
library(dplyr)
library(tidyr)
plan(multisession)

create_geoR_range <- function(range) {
  range / sqrt(8)
}
create_geoR_range(0.1)

true <- data.frame(
  name = c("phi", "sigma_O", "range"),
  value = c(2, 1, 0.1)
)

N <- 2000
.range <- true$value[true$name == "range"]
.phi <- true$value[true$name == "phi"]
.sigma_O <- true$value[true$name == "sigma_O"]

set.seed(1)
sim1 <- grf(N, cov.pars = c(.sigma_O, create_geoR_range(.range)), cov.model = "matern", kappa = 1)

dat <- data.frame(x = sim1$coords[, 1], y = sim1$coords[, 2], z = sim1$data, obs = sim1$data + rnorm(N, 0, .phi))
set.seed(1)
dat$fold_id <- sample(1:10, nrow(dat), replace = TRUE)

if (FALSE) {
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
}

# cutoffs <- seq(0.005, 0.15, length.out = 20)
cutoffs <- exp(seq(log(0.005), log(0.15), length.out = 20))


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

title <- paste0(
  "N = ", N,
  ", range = ", .range,
  ", obs. SD = ", .phi,
  ", sigma_O = ", .sigma_O
)
filename <- paste0("figures/", gsub(" ", "", title), ".pdf")

pivot_longer(ret, cols = -n) |>
  ggplot(aes(n, value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~name, scales = "free_y") +
  geom_hline(data = true, mapping = aes(yintercept = value), lty = 2) +
  ggsidekick::theme_sleek() +
  ggtitle(title, subtitle = "10-fold cross validation")

ggsave(filename, width = 8, height = 5)
