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

f <- "output/sim-grid-output.rds"
if (!file.exists(f)) {
  out <- furrr::future_pmap_dfr(torun, sim_and_fit)
  saveRDS(out, file = f)
} else {
  out <- readRDS(f)
}

f <- "output/sim-grid-output-small.rds"
if (!file.exists(f)) {
  torun <- expand.grid(
    .phi = c(0.05, 1.5, 3),
    .sigma_O = c(0.1, 1, 2),
    .range = c(0.05, 0.3, 0.6),
    N = c(1000)
  )
  set.seed(1)
  out_small <- furrr::future_pmap_dfr(torun, sim_and_fit, .options = furrr::furrr_options(seed = TRUE))
  saveRDS(out_small, file = f)
} else {
  out_small <- readRDS(f)
}

make_panels <- function(.term, data = out) {
  x <- pivot_longer(data, cols = c(leftout_loglik, phi_hat, sigma_O_hat, range_hat, rmse_true)) |>
    filter(name == .term) |>
    mutate(sigma_O_clean = paste0("Spatial SD: ", round(sigma_O, 2))) |>
    mutate(phi_clean = paste0("Observation\nSD: ", round(phi, 2)))

  if (.term == "sigma_O_hat") {
    x <- filter(x, value < 10)
    x <- filter(x, value > 0.01)
  }
  if (.term == "range_hat") {
    x <- filter(x, value < 20)
  }

  true <- select(x, -value, -seed, -name, -cutoff) |> distinct()

  g <- x |>
    ggplot(aes(mesh_n, value)) +
    facet_grid(phi_clean ~ sigma_O_clean, scales = "free") +
    ylab(.term) +
    xlab("Mesh knots") +
    # ggtitle(.term) +
    labs(colour = "Range") +
    geom_line(aes(colour = factor(round(range, 2)))) +
    ggsidekick::theme_sleek() +
    scale_x_continuous(breaks = seq(0, 1500, 500)) +
    scale_colour_viridis_d(end = 0.9, option = "C")

  if (.term == "range_hat") {
    g <- g +
      geom_hline(data = true, mapping = aes(yintercept = range, colour = factor(round(range, 2))), lty = 2) +
      ylab("Spatial range")
  }
  if (.term == "sigma_O_hat") {
    g <- g +
      geom_hline(data = true, mapping = aes(yintercept = sigma_O), lty = 2) +
      ylab("Spatial SD (sigma_O)")
  }
  if (.term == "phi_hat") {
    g <- g +
      geom_hline(data = true, mapping = aes(yintercept = phi), lty = 2) +
      ylab("Observation SD (phi)")
  }
  if (.term == "rmse_true") {
    g <- g + ylab("Predictive RMSE from truth")
  }
  if (.term == "leftout_loglik") {
    g <- g + ylab("Log predictive density")
  }

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

# version for main text?
x <- pivot_longer(out_small, cols = c(leftout_loglik, phi_hat, sigma_O_hat, range_hat, rmse_true)) |>
  filter(name == "leftout_loglik")
# true <- select(x, -value, -seed, -name, -cutoff) |> distinct()
x |>
  mutate(sigma_O_clean = paste0("Spatial SD: ", round(sigma_O, 2))) |>
  mutate(phi_clean = paste0("Observation\nSD: ", round(phi, 2))) |>
  ggplot(aes(mesh_n, value)) +
  facet_grid(phi_clean ~ sigma_O_clean, scales = "free") +
  ylab("Log predictive density") +
  xlab("Mesh vertices") +
  labs(colour = "Range") +
  geom_line(aes(colour = factor(round(range, 2)))) +
  ggsidekick::theme_sleek() +
  scale_colour_viridis_d(end = 0.9, option = "C") +
  tagger::tag_facets(tag_prefix = "(", position = list(x = 0.08, y = 0.89), tag = "panel") +
  theme(tagger.panel.tag.text = element_text(colour = "grey30"))
ggsave("figures/sim-grid-small-lpd.pdf", width = 6, height = 4)
ggsave("figures/sim-grid-small-lpd.png", width = 6, height = 4)
