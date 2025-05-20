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

sim_and_fit <- function(N, .phi, .sigma_O, .range, .seed = 1) {
  true <- data.frame(
    name = c("phi", "sigma_O", "range"),
    value = c(.phi, .sigma_O, .range)
  )

  set.seed(.seed)
  sim1 <- grf(N, cov.pars = c(.sigma_O, create_geoR_range(.range)), cov.model = "matern", kappa = 1)
  dat <- data.frame(x = sim1$coords[, 1], y = sim1$coords[, 2], z = sim1$data, obs = sim1$data + rnorm(N, 0, .phi))
  dat$fold_id <- sample(1:10, nrow(dat), replace = TRUE)

  cutoffs <- exp(seq(log(0.005), log(0.15), length.out = 10))

  ret <- tryCatch(
    {
      purrr::map_dfr(cutoffs, \(x) {
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
          mesh_n = mesh$mesh$n,
          phi_hat = phi,
          sigma_O_hat = sigma_O,
          range_hat = range,
          rmse_true = rmse_true
        )
      })
    },
    error = function(e) {
      data.frame(
        leftout_loglik = NA,
        cutoff = NA,
        mesh_n = NA,
        phi_hat = NA,
        sigma_O_hat = NA,
        range_hat = NA,
        rmse_true = NA
      )
    }
  )

  ret$range <- .range
  ret$phi <- .phi
  ret$sigma_O <- .sigma_O
  ret$N <- N
  ret$seed <- .seed
  ret
}

# sim_and_fit(N = 200, .phi = 1, .sigma_O = 1, .range = 0.3)

torun <- expand.grid(
  .phi = seq(0.05, 3, length.out = 4),
  .sigma_O = seq(0.1, 2, length.out = 4),
  .range = seq(0.05, 0.6, length.out = 4),
  N = c(1000)
)
nrow(torun)
plan(multicore)
out <- furrr::future_pmap_dfr(torun, sim_and_fit)

saveRDS(out, file = "output/sim-grid-output.rds")

out <- readRDS("output/sim-grid-output.rds")

glimpse(out)

pivot_longer(out, cols = c(leftout_loglik, phi_hat, sigma_O_hat, range_hat, rmse_true)) |>
  filter(name == "leftout_loglik") |>
  ggplot(aes(mesh_n, value)) +
  facet_grid(
    paste("phi:", round(phi, 2))
    ~ paste("sigma_O:", round(sigma_O, 2)),
    scales = "free"
  ) +
  geom_line(aes(colour = factor(round(range, 2)))) +
  ggsidekick::theme_sleek() +
  scale_colour_viridis_d(end = 0.9, option = "C")
ggsave("figures/sim-grid-ll.pdf", width = 9, height = 6)

pivot_longer(out, cols = c(leftout_loglik, phi_hat, sigma_O_hat, range_hat, rmse_true)) |>
  filter(name == "rmse_true") |>
  ggplot(aes(mesh_n, value)) +
  facet_grid(
    paste("phi:", round(phi, 2))
    ~ paste("sigma_O:", round(sigma_O, 2)),
    scales = "free"
  ) +
  geom_line(aes(colour = factor(round(range, 2)))) +
  ggsidekick::theme_sleek() +
  scale_colour_viridis_d(end = 0.9, option = "C")
ggsave("figures/sim-grid-rmse.pdf", width = 9, height = 6)

make_panels <- function(.term) {
  x <- pivot_longer(out, cols = c(leftout_loglik, phi_hat, sigma_O_hat, range_hat, rmse_true)) |>
    filter(name == .term)

  if (.term == "sigma_O_hat") {
    x <- filter(x, value < 20)
  }

  g <- x |> 
    ggplot(aes(mesh_n, value)) +
    facet_grid(
      paste("phi:", round(phi, 2))
      ~ paste("sigma_O:", round(sigma_O, 2)),
      scales = "free"
    ) +
    ylab(.term) +
    xlab("Mesh knots") +
    ggtitle(.term) +
    labs(colour = "Range") +
    geom_line(aes(colour = factor(round(range, 2)))) +
    ggsidekick::theme_sleek() +
    scale_colour_viridis_d(end = 0.9, option = "C")

  # if (.term == "phi_hat") {
  #   g <- g + geom_hline
  # }

  ggsave(paste0("figures/sim-grid-", .term, ".pdf"),
    width = 9, height = 6
  )
}
make_panels("rmse_true")
make_panels("phi_hat")
make_panels("sigma_O_hat")
make_panels("range_hat")
make_panels("leftout_loglik")


# title <- paste0(
#   "N = ", N,
#   ", range = ", .range,
#   ", obs. SD = ", .phi,
#   ", sigma_O = ", .sigma_O
# )
# filename <- paste0("figures/", gsub(" ", "", title), ".pdf")
#
# pivot_longer(ret, cols = -n) |>
#   ggplot(aes(n, value)) +
#   geom_point() +
#   geom_line() +
#   facet_wrap(~name, scales = "free_y") +
#   geom_hline(data = true, mapping = aes(yintercept = value), lty = 2) +
#   ggsidekick::theme_sleek() +
#   ggtitle(title, subtitle = "10-fold cross validation")
