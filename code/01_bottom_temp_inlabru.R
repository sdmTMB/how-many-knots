library(sf)
library(inlabru)
library(INLA)
library(dplyr)
library(ggplot2)
library(future)
library(sdmTMB)
plan(multisession)
library(sp)
library(patchwork)
library(mgcv)

haul <- readRDS("data/haul_cleaned.rds")
haul$date_formatted <- as.Date(as.numeric(haul$date_formatted), origin = "1970-01-01")
haul$year <- lubridate::year(haul$date_formatted)
haul$month <- lubridate::month(haul$date_formatted)
haul$yday <- lubridate::yday(haul$date_formatted)
haul$zday <- as.numeric(scale(haul$yday))
haul$zday2 <- haul$zday^2
haul$log_depth_scaled <- as.numeric(haul$log_depth_scaled)
haul$log_depth_scaled2 <- as.numeric(haul$log_depth_scaled2)
haul_new <- dplyr::filter(haul, year == 2018, !is.na(temperature_at_gear_c_der))
haul_new <- dplyr::filter(haul_new, !is.na(log_depth_scaled))
haul_new <- dplyr::filter(haul_new, !is.na(zday))
set.seed(123)
# haul_new$fold_id <- rep(1:10, length.out = nrow(haul_new))
haul_new$fold_id <- sample(1:10, size = nrow(haul_new), replace = TRUE)
haul_new$pred_train <- NA
haul <- haul_new

ns <- dplyr::group_by(haul_new, fold_id) |>
  dplyr::summarise(n = n(), n_other = nrow(haul_new) - n) |>
  dplyr::ungroup() |>
  dplyr::summarise(mean_n = mean(n), mean_other = mean(n_other))

sp::coordinates(haul) <- c("X", "Y")
boundary <- INLA::inla.nonconvex.hull(sp::coordinates(haul),
  convex = -0.05
)

run_cv <- function(cutoff, folds = 1:10, run_inla = TRUE, range_gt = 100, sigma_lt = 10, run_mgcv = TRUE, run_noprior = TRUE) {
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

  if (!is.na(range_gt)) {
    matern <- INLA::inla.spde2.pcmatern(
      mesh,
      prior.sigma = c(sigma_lt, 0.05),
      prior.range = c(range_gt, 0.05)
    )

    components <- temperature_at_gear_c_der ~ zday + zday2 + Intercept(1) +
      field(main = geometry, model = matern)
  }

  df$ll_test <- 0
  df$ll_train <- 0
  df$ll_test_sdmTMB <- 0
  df$ll_train_sdmTMB <- 0
  df$ll_test_sdmTMB_noprior <- 0
  df$ll_train_sdmTMB_noprior <- 0
  df$ll_test_mgcv <- 0
  df$ll_train_mgcv <- 0

  df$ss_sdmTMB <- df$ss_sdmTMB_noprior <- df$ss_mgcv <- df$ss_inla <- 0

  df$range_gt <- range_gt
  df$sigma_lt <- sigma_lt

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
    if (!is.na(range_gt)) {
      fit_train_sdmTMB <- tryCatch(sdmTMB(
        temperature_at_gear_c_der ~ zday + I(zday^2),
        data = this_dat,
        mesh = this_mesh,
        priors = sdmTMBpriors(pc_matern(range_gt = range_gt, sigma_lt = sigma_lt))
      ), error = function(e) {
        return(NULL)
      })
      if (is.null(fit_train_sdmTMB)) {
        break
      }
    }

    if (run_noprior) {
      fit_train_sdmTMB_noprior <- tryCatch(sdmTMB(
        temperature_at_gear_c_der ~ zday + I(zday^2),
        data = this_dat,
        mesh = this_mesh,
      ), error = function(e) {
        return(NULL)
      })
      if (is.null(fit_train_sdmTMB_noprior)) {
        break
      }
    }

    if (run_mgcv) {
      mgcv_k <- round(df$n / 3)
      if (mgcv_k > 350) run_mgcv <- FALSE
      df$mgcv_k <- mgcv_k
      fit_train_mgcv <- tryCatch(mgcv::gam(
        temperature_at_gear_c_der ~ zday + I(zday^2) + s(X, Y, k = mgcv_k),
        data = this_dat,
        mesh = this_mesh,
      ), error = function(e) {
        return(NULL)
      })
      if (is.null(fit_train_mgcv)) {
        run_mgcv <- FALSE # failed!
      }
    }

    if (run_mgcv) {
      df$mgcv_edf <- summary(fit_train_mgcv)$s.table[, "edf"]
    }

    if (!is.na(range_gt)) {
      df$sdmTMB_edf <- cAIC(fit_train_sdmTMB, "EDF")[[1]]
      df$sdmTMB_cAIC <- cAIC(fit_train_sdmTMB, "cAIC")[[1]]
    }
    if (run_noprior) {
      df$sdmTMB_edf_noprior <- cAIC(fit_train_sdmTMB_noprior, "EDF")[[1]]
      df$sdmTMB_cAIC_noprior <- cAIC(fit_train_sdmTMB_noprior, "cAIC")[[1]]
    }

    if (run_inla) {
      pred_train <- tryCatch(
        {
          predict(fit_train,
            newdata = haul[haul$fold_id != ii, , drop = FALSE],
            ~ zday + zday2 + field + Intercept
          )
        },
        error = function(e) {
          message("    -> pred_train failed: ", e$message)
          return(NULL)
        }
      )
    }

    if (!is.na(range_gt)) {
      pred_train_sdmTMB <- predict(fit_train_sdmTMB, newdata = NULL)
    }
    if (run_noprior) {
      pred_train_sdmTMB_noprior <- predict(fit_train_sdmTMB_noprior, newdata = NULL)
    }
    if (run_mgcv) {
      pred_train_mgcv <- predict(fit_train_mgcv)
    }

    if (run_inla) {
      pred_test <- tryCatch(
        {
          predict(fit_train,
            newdata = haul[haul$fold_id == ii, , drop = FALSE],
            ~ zday + zday2 + field + Intercept
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
    if (!is.na(range_gt)) {
      pred_test_sdmTMB <- predict(fit_train_sdmTMB, newdata = nd)
    }
    if (run_noprior) {
      pred_test_sdmTMB_noprior <- predict(fit_train_sdmTMB_noprior, newdata = nd)
    }

    if (run_mgcv) {
      pred_test_mgcv <- predict(fit_train_mgcv, newdata = nd)
    }

    if (FALSE) {
      plot(pred_test_sdmTMB$est, pred_test$mean)
      abline(0, 1)
      plot(pred_train_sdmTMB$est, pred_train$mean)
      abline(0, 1)
      plot(pred_train_sdmTMB$est, pred_train_mgcv)
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

    if (!is.na(range_gt)) {
      phi <- exp(get_pars(fit_train_sdmTMB)$ln_phi)
    }
    if (run_noprior) {
      phi_noprior <- exp(get_pars(fit_train_sdmTMB_noprior)$ln_phi)
    }
    if (run_mgcv) {
      phi_mgcv <- summary(fit_train_mgcv)$scale
    }

    if (run_inla) {
      # if any part failed, mark the fold as failed
      if (is.null(pred_train) || is.null(pred_test) || is.na(tau)) {
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

      df$ss_inla <- df$ss_inla +
        sum((haul$temperature_at_gear_c_der[haul$fold_id == ii] - pred_test$mean)^2)
    }

    if (!is.na(range_gt)) {
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

      df$ss_sdmTMB <- df$ss_sdmTMB +
        sum((nd$temperature_at_gear_c_der - pred_test_sdmTMB$est)^2)
    }

    if (run_noprior) {
      df$ll_train_sdmTMB_noprior <- df$ll_train_sdmTMB_noprior +
        sum(dnorm(
          this_dat$temperature_at_gear_c_der,
          mean = pred_train_sdmTMB_noprior$est,
          sd = phi_noprior,
          log = TRUE
        ))

      df$ll_test_sdmTMB_noprior <- df$ll_test_sdmTMB_noprior +
        sum(dnorm(
          nd$temperature_at_gear_c_der,
          mean = pred_test_sdmTMB_noprior$est,
          sd = phi_noprior,
          log = TRUE
        ))

      df$ss_sdmTMB_noprior <- df$ss_sdmTMB_noprior +
        sum((nd$temperature_at_gear_c_der - pred_test_sdmTMB_noprior$est)^2)
      df$sigma_sdmTMB_noprior <- phi_noprior
    }

    if (run_mgcv) {
      df$ll_train_mgcv <- df$ll_train_mgcv +
        sum(dnorm(
          this_dat$temperature_at_gear_c_der,
          mean = pred_train_mgcv,
          sd = phi_mgcv,
          log = TRUE
        ))

      df$ll_test_mgcv <- df$ll_test_mgcv +
        sum(dnorm(
          nd$temperature_at_gear_c_der,
          mean = pred_test_mgcv,
          sd = phi_mgcv,
          log = TRUE
        ))

      df$ss_mgcv <- df$ss_mgcv +
        sum((nd$temperature_at_gear_c_der - pred_test_mgcv)^2)
    }
  } # end CV folds loop

  df
}

out <- run_cv(cutoff = 100, folds = 1:2, run_noprior = FALSE, run_inla = FALSE, run_mgcv = FALSE)

out <- run_cv(cutoff = 100, folds = 1:2)

# torun <- data.frame(cutoff = exp(seq(log(8), log(500), length.out = 30)))
torun <- data.frame(cutoff = exp(seq(log(8), log(500), length.out = 10)))
nrow(torun)
torun$cutoff

plan(multisession, workers = 10)
out <- furrr::future_pmap(torun, run_cv)
plan(sequential)
saveRDS(out, file = "output/01_loocv_temp_inlabru_mgcv.rds")
out <- readRDS("output/01_loocv_temp_inlabru_mgcv.rds")

out2 <- out |> bind_rows()
head(out2)

N <- nrow(haul)

theme_set(ggsidekick::theme_sleek())

out3 <- tidyr::pivot_longer(select(out2, cutoff, n, ll_test:ll_train_mgcv), cols = ll_test:ll_train_mgcv)

x <- out3 |>
  mutate(test = grepl("test", name)) |>
  mutate(test_char = ifelse(test, "Test", "Train")) |>
  mutate(name = gsub("ll_test$", "ll_test_INLA", name)) |>
  mutate(name = gsub("ll_train$", "ll_train_INLA", name)) |>
  mutate(model = stringr::str_remove(name, "ll_(test|train)_")) |>
  mutate(model = gsub("_", " ", model)) |>
  filter(value != 0) |>
  filter(!(value > -300 & test)) |>
  filter(!(value > -2350 & !test)) |>
  mutate(value = ifelse(test, value / N, value / (10 * N))) |>
  filter(cutoff < 200)

# x |>
#   ggplot(aes(n, value, colour = model)) +
#   facet_grid(~test_char, scales = "free_y") +
#   geom_line() +
#   scale_colour_brewer(palette = "Set2") +
#   labs(colour = "Model", y = "Log density") +
#   xlab("Knots") +
#   coord_cartesian(ylim = c(-2, NA)) +
#   theme(legend.position = "top")

mods <- c("sdmTMB", "sdmTMB noprior", "INLA", "mgcv")
cols <- RColorBrewer::brewer.pal(4, "Set2")
names(cols) <- mods

g1 <- x |>
  filter(!grepl("mgcv", model)) |>
  ggplot(aes(n, value, colour = model)) +
  facet_wrap(~test_char, scales = "free_y") +
  geom_line() +
  scale_colour_manual(values = cols, drop = FALSE) +
  labs(colour = "Model", y = "Log density") +
  xlab("Mesh vertices") +
  theme(legend.position = "top")

g2 <- x |>
  filter(!grepl("mgcv", model)) |>
  ggplot(aes(cutoff, value, colour = model)) +
  facet_wrap(~test_char, scales = "free_y") +
  geom_line() +
  scale_colour_manual(values = cols, drop = FALSE) +
  labs(colour = "Model", y = "Log density") +
  xlab("Cutoff distance") +
  theme(legend.position = "top")

g3 <- x |>
  filter(grepl("mgcv", model)) |>
  ggplot(aes(n / 3, value, colour = model)) +
  facet_wrap(~test_char, scales = "free_y") +
  geom_line() +
  scale_colour_manual(values = cols, drop = FALSE) +
  scale_x_continuous(breaks = seq(0, 300, 50)) +
  labs(colour = "Model", y = "Log density") +
  xlab("Smoother basis dimension (k)") +
  theme(legend.position = "top")

g2 / g1 / g3 +
  plot_layout(axes = "collect", guides = "collect") &
  theme(legend.position = "right")

ggsave("figures/temperature-inla-sdmTMB-mgcv.pdf", width = 7.5, height = 7)

# RMSE??
g00 <- x |>
  filter(!grepl("mgcv", model)) |>
  filter(test) |>
  mutate(value = value * N) |>
  mutate(value = value - max(value)) |>
  ggplot(aes(n, value, colour = model)) +
  facet_wrap(~test_char, scales = "free_y") +
  geom_line() +
  scale_colour_manual(values = cols, drop = FALSE) +
  labs(colour = "Model", y = "Relative log density") +
  xlab("Mesh vertices") +
  theme(legend.position = "top")

g0 <- tidyr::pivot_longer(out2, cols = ss_inla:ss_sdmTMB) |>
  mutate(rmse = sqrt(value / N)) |>
  ggplot(aes(n, rmse, colour = name)) +
  # scale_colour_manual(values = cols, drop = FALSE) +
  xlab("Mesh vertices") +
  labs(colour = "Model", y = "RMSE") +
  geom_line()

g4 <- out2 |>
  ggplot(aes(n, sigma_sdmTMB_noprior)) +
  xlab("Mesh vertices") +
  ylab("Observation error SD (sdmTMB)") +
  geom_line()

g00 / g0 / g4 + plot_layout(axes = "collect") # , guides = "collect")

out4 <- tidyr::pivot_longer(select(out2, cutoff, n, mgcv_edf:sdmTMB_cAIC_noprior), cols = mgcv_edf:sdmTMB_cAIC_noprior)

out4 |>
  filter(!grepl("noprior", name)) |>
  mutate(type = ifelse(grepl("edf", name), "EDF", "cAIC")) |>
  ggplot(aes(n, value, colour = name)) +
  geom_line() +
  facet_wrap(~type, scales = "free_y")
ggsave("figures/temperature-sdmTMB-mgcv-edf.pdf", width = 7.5, height = 3)

# -------------------------------------------------
# across a grid of PC priors

torun <- expand.grid(
  cutoff = exp(seq(log(8), log(150), length.out = 10)),
  range_gt = seq(20, 1200, length.out = 8),
  sigma_lt = seq(2, 12, length.out = 3),
  run_inla = FALSE,
  run_mgcv = FALSE,
  run_noprior = FALSE
)
nrow(torun)

plan(multisession, workers = 10)
out_pc <- furrr::future_pmap(torun, run_cv)
plan(sequential)
saveRDS(out_pc, file = "output/01_loocv_temp_pc_priors.rds")

torun2 <- expand.grid(
  cutoff = exp(seq(log(8), log(150), length.out = 10)),
  range_gt = NA,
  sigma_lt = NA,
  run_inla = FALSE,
  run_mgcv = FALSE,
  run_noprior = TRUE
)
nrow(torun2)
plan(multisession, workers = 10)
out_pc2 <- furrr::future_pmap(torun2, run_cv)
plan(sequential)
saveRDS(out_pc2, file = "output/01_loocv_temp_pc_priors2.rds")

out_pc <- readRDS("output/01_loocv_temp_pc_priors.rds")
out_pc2 <- readRDS("output/01_loocv_temp_pc_priors2.rds")

out_pc <- out_pc |>
  bind_rows() |>
  mutate(type = "priors")

out_pc2 <- out_pc2 |>
  bind_rows() |>
  mutate(type = "no priors") |>
  select(-ll_train_sdmTMB, -ll_test_sdmTMB) |>
  rename(ll_train_sdmTMB = ll_train_sdmTMB_noprior, ll_test_sdmTMB = ll_test_sdmTMB_noprior)

out_pc <- bind_rows(out_pc, out_pc2)

out_pc_long <- tidyr::pivot_longer(select(out_pc, type, cutoff, n, range_gt, sigma_lt, ll_test_sdmTMB:ll_train_sdmTMB), cols = ll_test_sdmTMB:ll_train_sdmTMB)

x <- out_pc_long |>
  mutate(test = grepl("test", name)) |>
  mutate(test_char = ifelse(test, "Test", "Train")) |>
  mutate(name = gsub("ll_test$", "ll_test_INLA", name)) |>
  mutate(name = gsub("ll_train$", "ll_train_INLA", name)) |>
  mutate(model = stringr::str_remove(name, "ll_(test|train)_")) |>
  mutate(model = gsub("_", " ", model)) |>
  filter(value != 0) |>
  filter(!(value > -300 & test)) |>
  filter(!(value > -2350 & !test)) |>
  mutate(value = ifelse(test, value / N, value / (10 * N))) |>
  mutate(prior = paste0("Range gt ", range_gt, ";", "Sigma lt ", sigma_lt))

x_noprior <- filter(x, type == "no priors")
x_prior <- filter(x, type == "priors")

g1 <- x_prior |>
  ggplot(aes(n, value, group = prior, colour = sigma_lt)) +
  facet_wrap(~test_char, scales = "free_y") +
  geom_line() +
  scale_colour_viridis_c(option = "G") +
  labs(colour = "Matérn GMRF PC prior:\nPr(Sigma < x) = 0.95", y = "Log density") +
  xlab("Knots") +
  theme(legend.position = "top") +
  geom_line(data = x_noprior, colour = "red", lwd = 1, lty = 2)

g2 <- x_prior |>
  ggplot(aes(n, value, group = prior, colour = range_gt)) +
  facet_wrap(~test_char, scales = "free_y") +
  geom_line() +
  scale_colour_viridis_c(option = "D") +
  labs(colour = "Matérn GMRF PC prior:\nPr(Range > x) = 0.95", y = "Log density") +
  xlab("Knots") +
  theme(legend.position = "top") +
  geom_line(data = x_noprior, colour = "red", lwd = 1, lty = 2)

g2 / g1

ggsave("figures/temperature-sdmTMB-pc-priors.pdf", width = 7.5, height = 7)
