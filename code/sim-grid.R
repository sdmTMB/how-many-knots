# questions:
# - what ranges should we evaluate here if we want to ballpark match the NWFSC survey?
# - survey X range = ~700 km
# - survey Y range = ~1800 km
# - survey mean nearest neighbour distance = 7.146
# - survey mean points per year = 645
# - grid is 12,000 cells 2 nm x 1.5 nm, 1 nm = 1.852 km

# - for real data: arrowtooth range = ~ 500 km in long direction and ~ 200 in short direction
# - sablefish is much smaller: maybe 100 in long direction and 50 in short direction? very rough

# # so at max, range is this fraction of Y
# 500/1800

# # so at min, range is this fraction of X
# 50/700

# 100/1800
# 300/700

# so 0.05 to 0.3

# in terms of sampling density, mean NN dist:
# 7 / 700
# 7 / 1800

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
  sim1 <- grf(N, cov.pars = c(.sigma_O^2, create_geoR_range(.range)), cov.model = "matern", kappa = 1)
  dat <- data.frame(x = sim1$coords[, 1], y = sim1$coords[, 2], z = sim1$data, obs = exp(sim1$data + rnorm(N, 0, .phi)))
  dat$fold_id <- sample(1:10, nrow(dat), replace = TRUE)

#  browser()
#  mean(spatstat.geom::nndist(X = dat$x,Y = dat$y))

  cutoffs <- exp(seq(log(0.005), log(0.15), length.out = 10))

  ret <- tryCatch(
    {
      purrr::map_dfr(cutoffs, \(x) {
        mesh <- make_mesh(dat, c("x", "y"), cutoff = x)
        fit_cv <- sdmTMB_cv(
          obs ~ 1,
          data = dat, mesh = mesh,
          k_folds = 10, parallel = FALSE,
          fold_ids = "fold_id", family = lognormal(link = "log")
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
        rmse_true <- sqrt(mean((exp(fit_cv$data$z) - fit_cv$data$cv_predicted)^2))
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
  .phi = seq(0.05, 2, length.out = 4),
  .sigma_O = seq(0.1, 2, length.out = 4),
  .range = seq(0.03, 0.3, length.out = 4),
  N = c(1000)
)
nrow(torun)
plan(multisession, workers = 5L)

# f <- "output/sim-grid-output.rds"
# if (!file.exists(f)) {
#   out <- furrr::future_pmap_dfr(torun, sim_and_fit)
# #  out <- purrr::pmap_dfr(torun, sim_and_fit)
#   saveRDS(out, file = f)
# } else {
#   out <- readRDS(f)
# }

f <- "output/sim-grid-output-small-lognormal.rds"
if (!file.exists(f)) {
  torun <- expand.grid(
    .phi = c(0.5, 1, 2),
    .sigma_O = c(0.5, 1, 2),
    .range = c(0.05, 0.1, 0.2),
    N = c(1000)
  )
  set.seed(1)
  out_small <- furrr::future_pmap_dfr(torun, sim_and_fit, .options = furrr::furrr_options(seed = TRUE))
  saveRDS(out_small, file = f)
} else {
  out_small <- readRDS(f)
}

make_panels <- function(.term, data = out_small) {
  x <- pivot_longer(data, cols = c(leftout_loglik, phi_hat, sigma_O_hat, range_hat, rmse_true)) |>
    filter(name == .term) |>
    mutate(sigma_O_clean = paste0("Spatial SD: ", round(sigma_O, 2))) |>
    mutate(phi_clean = paste0("Observation\nlog-SD: ", round(phi, 2)))

  if (.term == "sigma_O_hat") {
    x <- filter(x, value < 10)
    x <- filter(x, value > 0.01)
  }
  if (.term == "range_hat") {
    x <- filter(x, value < 20)
  }

  if (.term %in% c("leftout_loglik", "rmse_true")) {
    x <- x |>
      group_by(phi, sigma_O, range, seed) |>
      mutate(value = if (.term == "leftout_loglik") {
        value - max(value, na.rm = TRUE)
      } else {
        value - min(value, na.rm = TRUE)
      }) |>
      ungroup()
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
      ylab("Observation log-SD (phi)")
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
    width = 6, height = 4
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
ggplot(
  data = x |>
    group_by(phi, sigma_O, range, seed) |>
    # mutate(value = (value - min(value, na.rm = TRUE)) /
    #   (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))) |>
    mutate(value = value - max(value, na.rm = TRUE)) |>
    ungroup() |>
    mutate(sigma_O_clean = paste0("Spatial SD: ", round(sigma_O, 2))) |>
    mutate(phi_clean = paste0("Observation\nlog-SD: ", round(phi, 2))),
  aes(mesh_n, value)
) +
  facet_grid(phi_clean ~ sigma_O_clean, scales = "free") +
  ylab("Log predictive density") +
  xlab("Mesh vertices") +
  labs(colour = "Range") +
  geom_line(aes(colour = factor(round(range, 2)))) +
  ggsidekick::theme_sleek() +
  scale_colour_viridis_d(end = 0.9, option = "C") +
  tagger::tag_facets(tag_prefix = "(", position = "bl", tag = "panel") +
  theme(tagger.panel.tag.text = element_text(colour = "grey30"))
ggsave("figures/sim-grid-small-lpd.pdf", width = 6, height = 4)
ggsave("figures/sim-grid-small-lpd.png", width = 6, height = 4)

out$scenario <- paste0("phi = ", round(out$phi, 2), ", sigma_O = ", round(out$sigma_O, 2), ", range = ", round(out$range, 2))

make_scenario_plot <- function(scen) {
  x1 <- filter(out, scenario == scen) |>
    pivot_longer(cols = c(ends_with("hat"), rmse_true, leftout_loglik))

  trues <- filter(out, scenario == scen) |>
    pivot_longer(cols = c(phi, sigma_O, range)) |>
    rename(true_value = value) |>
    mutate(name = paste0(name, "_hat")) |>
    select(mesh_n, name, true_value)

  x1 <- x1 |> left_join(trues) |>
    mutate(name = factor(name, levels = c("leftout_loglik", "rmse_true", "phi_hat", "range_hat", "sigma_O_hat")))

  x1 |>
    ggplot(aes(mesh_n, value)) + geom_line() +
    geom_line(aes(y = true_value), lty = 2) +
    facet_wrap(~name, scales = "free_y", nrow = 5) +
    ggtitle(scen) +
    ggsidekick::theme_sleek() + ylab("Value") + xlab("Mesh vertices")
}

# unique(out$scenario)
# make_scenario_plot("phi = 1.35, sigma_O = 2, range = 0.23")
# make_scenario_plot("phi = 1.03, sigma_O = 2, range = 0.6")

# # an extreme top right scenario:
# make_scenario_plot("phi = 0.05, sigma_O = 2, range = 0.23")

# # more obs. error
# # now phi asymptotes
# make_scenario_plot("phi = 1.03, sigma_O = 2, range = 0.23")

# make_scenario_plot("phi = 0.05, sigma_O = 0.73, range = 0.6")

# make_scenario_plot("phi = 2.02, sigma_O = 0.1, range = 0.6")

filter(out, phi_hat > phi)
