library(sf)
library(sdmTMB)
library(lubridate)
library(dplyr)
library(future)
library(viridis)
library(ggplot2)
haul <- readRDS("data/haul_cleaned.rds")
catch <- readRDS("data/catch_cleaned.rds")
haul$trawl_id <- as.numeric(haul$trawl_id)

run_cv <- function(cutoff, bin_width, species, seed = NULL, do_full_fit = FALSE, parallel = TRUE) {
  catch_sub <- dplyr::filter(catch, common_name == species)

  print(species)
  print(cutoff)

  # Join catch and haul data
  joined_dat <- haul |>
    dplyr::filter(trawl_id %in% catch$trawl_id) |>
    left_join(catch_sub[, c("date", "cpue_kg_km2", "depth_m", "trawl_id")],
      by = "trawl_id"
    ) |>
    dplyr::filter(!is.na(cpue_kg_km2)) |>
    dplyr::mutate(present = ifelse(cpue_kg_km2 > 0, 1, 0))
  joined_dat$date <- lubridate::ymd(as.POSIXct.Date(joined_dat$date))
  joined_dat$year <- year(joined_dat$date)
  joined_dat$fyear <- as.factor(joined_dat$year)
  joined_dat$yday <- yday(joined_dat$date)

  if (!is.na(bin_width)) {
    breaks <- seq(min(joined_dat$Y), max(joined_dat$Y), by = bin_width)
    bins <- cut(joined_dat$Y, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    bins[which(is.na(bins))] <- max(bins, na.rm = TRUE) + 1
    bins_to_folds <- data.frame(bin = 1:max(bins), fold = rep(1:10, 1000)[1:max(bins)])
    joined_dat$bin <- bins
    joined_dat <- dplyr::left_join(joined_dat, bins_to_folds)
  } else {
    if (!is.null(seed)) set.seed(seed)
    folds <- sample(seq_len(10L), size = nrow(joined_dat), replace = TRUE)
    joined_dat$fold <- folds
    joined_dat$bin <- NA
  }

  sp::coordinates(joined_dat) <- c("X", "Y")

  # create boundary, same for all meshes
  boundary <- fmesher::fm_nonconvex_hull_inla(sp::coordinates(joined_dat),
    convex = -0.05 # Negative values are interpreted as fractions of the approximate initial set
  )

  inla_mesh <- fmesher::fm_mesh_2d_inla(
    loc = sp::coordinates(joined_dat),
    boundary = boundary,
    offset = c(-0.05, -0.05),
    max.n = 5000,
    max.n.strict = 5000,
    cutoff = cutoff,
    max.edge = c(100, 500)
  )
  joined_dat <- as.data.frame(joined_dat)
  mesh <- make_mesh(joined_dat, c("X", "Y"), mesh = inla_mesh)

  if (!do_full_fit) {
    cat("Cross validation model fit\n")
    fit_cv <- try(sdmTMB_cv(
      cpue_kg_km2 ~ -1 + as.factor(year) + poly(log_depth_scaled, 2),
      data = joined_dat,
      mesh = mesh,
      time = "year",
      spatial = "on",
      spatiotemporal = "off",
      anisotropy = TRUE,
      share_range = TRUE,
      family = tweedie(),
      parallel = parallel,
      k_folds = length(unique(joined_dat$fold)),
      fold_ids = joined_dat$fold
    ), silent = TRUE)
  }

  if (do_full_fit) {
    cat("Full model fit\n")
    fit_full <- try(sdmTMB(
      cpue_kg_km2 ~ -1 + as.factor(year) + poly(log_depth_scaled, 2),
      data = joined_dat,
      mesh = mesh,
      time = "year",
      spatial = "on",
      spatiotemporal = "off",
      anisotropy = TRUE,
      share_range = TRUE,
      family = tweedie()
    ), silent = TRUE)
    if (class(fit_full) != "try-error") {
      tidy_pars <- tidy(fit_full, "ran_pars")
      tidy_pars <- filter(tidy_pars, term != "range")
      aniso <- sdmTMB:::print_anisotropy(fit_full, return_dat = TRUE)
      tidy_pars <- bind_rows(
        tibble(term = "range_a", estimate = as.numeric(aniso$sp$a)),
        tidy_pars,
        tibble(term = "range_b", estimate = as.numeric(aniso$sp$b)),
        tibble(term = "angle", estimate = as.numeric(aniso$sp$degree))
      )
      s <- sanity(fit_full, silent = TRUE, gradient_thresh = 0.01)
      tidy_pars$sanity <- s$hessian_ok && s$eigen_values_ok && s$nlminb_ok && s$gradients_ok
      tidy_pars$cutoff <- cutoff
      tidy_pars$bin_width <- bin_width
      tidy_pars$species <- species
      tidy_pars$n <- mesh$mesh$n
      tidy_pars$loglike <- as.numeric(logLik(fit_full))
      tidy_pars$cAIC <- tryCatch(as.numeric(cAIC(fit_full, what = "cAIC")), error = function(e) NA)
      edf <- tryCatch(cAIC(fit_full, what = "EDF"), error = function(e) NA)
      if (all(!is.na(edf))) {
        # tidy_pars$edf_epsilon <- edf[["epsilon_st"]]
        tidy_pars$edf_omega <- edf[["omega_s"]]
      }
      return(tidy_pars)
    } else {
      return(
        data.frame(cutoff = cutoff, bin_width = bin_width, species = species, converged = FALSE, n = mesh$mesh$n)
      )
    }
  }

  if (class(fit_cv) != "try-error") {
    ll <- data.frame(cutoff = cutoff, bin_width = bin_width, species = species)
    ll$n <- mesh$mesh$n
    ll$converged <- fit_cv$converged
    ll$test_dens_ll_mean <- mean(fit_cv$data$cv_loglik)
    ll$test_dens_ll_sum <- sum(fit_cv$data$cv_loglik)
    ll$converged <- fit_cv$converged
    ll$seed <- seed
    ll_train <- list()
    for (ii in seq_len(max(joined_dat$fold))) {
      # predict to whole dataset, test and train
      pred <- predict(fit_cv$models[[ii]]) |> dplyr::filter(fold != ii)
      joined_dat_train <- dplyr::filter(joined_dat, fold != ii)
      p <- plogis(fit_cv$models[[ii]]$model$par[["thetaf"]]) + 1
      phi <- exp(fit_cv$models[[ii]]$model$par[["ln_phi"]])
      ll_train <- fishMod::dTweedie(y = joined_dat_train$cpue_kg_km2, p = p, mu = exp(pred$est), phi = phi)
    }
    ll$train_dens_ll_n_mean <- mean(unlist(ll_train))
    ll$train_dens_ll_sum <- sum(unlist(ll_train))
    return(ll)
  } else {
    return(
      data.frame(cutoff = cutoff, bin_width = bin_width, species = species, converged = FALSE, n = mesh$mesh$n)
    )
  }
}

# if (FALSE) {
#   df <- expand.grid(
#     cutoff = 40,
#     bin_width = 40,
#     species = "sablefish"
#   )
#   tictoc::tic()
#   plan(multisession)
#   out <- purrr::pmap(df[1, ], run_cv, parallel = TRUE)
#   tictoc::toc()
#   out2 <- purrr::pmap(df[1, ], run_cv, do_full_fit = TRUE, parallel = FALSE)
#
#   df <- expand.grid(
#     cutoff = 40,
#     bin_width = NA,
#     seed = 1,
#     species = "sablefish"
#   )
#   tictoc::tic()
#   plan(multisession)
#   out <- purrr::pmap(df[1, ], run_cv, parallel = TRUE)
#   tictoc::toc()
#   out2 <- purrr::pmap(df[1, ], run_cv, do_full_fit = TRUE, parallel = FALSE)
# }

# random k-fold: -----------------------------------------------------------------
df <- expand.grid(
  cutoff = round(exp(seq(log(10), log(175), length.out = 25))),
  bin_width = NA,
  seed = c(281, 92110, 27293, 8282, 812938),
  species = c("sablefish", "arrowtooth flounder", "petrale sole")
)
nrow(df)
f <- "output/gf-cv-random-out.rds"
if (!file.exists(f)) {
  plan(multicore, workers = 80L)
  out <- furrr::future_pmap(df, run_cv, parallel = FALSE)
  plan(sequential)
  saveRDS(out, file = f)
} else {
  out <- readRDS(f)
}

f <- "output/gf-cv-full-fit.rds"
if (!file.exists(f)) {
  plan(multicore, workers = 80L)
  df <- expand.grid(
    cutoff = round(exp(seq(log(10), log(175), length.out = 25))),
    bin_width = 10,
    seed = 123,
    species = c("sablefish", "arrowtooth flounder", "petrale sole")
  )
  out2 <- furrr::future_pmap(df, run_cv, do_full_fit = TRUE, parallel = FALSE)
  saveRDS(out2, file = f)
  plan(sequential)
} else {
  out2 <- readRDS(f)
}

if (FALSE) {
  # strip width blocked CV: ---------------------------------------------------------
  df <- expand.grid(
    cutoff = round(exp(seq(log(10), log(175), length.out = 25))),
    bin_width = seq(10, 130, by = 40),
    seed = 123,
    species = c("sablefish", "arrowtooth flounder", "petrale sole")
  )
  nrow(df)
  plan(multicore, workers = 80L)
  out <- furrr::future_pmap(df, run_cv, parallel = FALSE)
  saveRDS(out, file = "output/gf-cv-blocked-out.rds")
  plan(sequential)

  d <- readRDS("output/gf-cv-blocked-out.rds")
  d <- bind_rows(d)
  d |>
    ggplot(aes(n, test_dens_ll_sum, colour = factor(bin_width))) +
    geom_line() +
    facet_grid(species ~ bin_width, scales = "free_y") +
    geom_smooth(se = FALSE)
}

out_df <- bind_rows(out) |> filter(species != "yelloweye rockfish")

make_panel <- function(dat) {
  if (dat$species[[1]] == "lingcod") {
    dat <- filter(dat, test_dens_ll_sum > -100000)
  }

  dat |>
    select(n, converged, cutoff, species, seed, test_dens_ll_sum, train_dens_ll_sum) |>
    tidyr::pivot_longer(cols = c(test_dens_ll_sum, train_dens_ll_sum), names_to = "ll_type") |>
    filter(n > 80, converged) |>
    mutate(ll_type = ifelse(grepl("test", ll_type), "Out of sample", "In sample")) |>
    group_by(ll_type, species, seed) |>
    # mutate(
    #   value = value - value[n()]
    # ) |>
    group_by(n, cutoff, ll_type, species) |>
    summarize(est = mean(value), lwr = min(value), upr = max(value)) |>
    group_by(ll_type, species) |>
    # mutate(
    #   lwr = lwr - max(est),
    #   upr = upr - max(est),
    #   est = est - max(est)
    # ) |>
    ggplot(aes(n, est, colour = cutoff)) +
    scale_colour_viridis(direction = -1, end = 0.9, limits = c(min(out_df$cutoff), max(out_df$cutoff))) +
    geom_line() +
    # facet_wrap(species~ll_type, scales = "free_y", ncol = 2) +
    facet_wrap(~ll_type, scales = "free_y", ncol = 2) +
    # geom_point(aes(size = cutoff), pch = 21) +
    geom_point(pch = 21) +
    geom_linerange(aes(ymin = lwr, ymax = upr)) +
    geom_smooth(se = FALSE, method = "gam", colour = "grey50") +
    ylab("Log predictive density") +
    xlab("Mesh vertices") +
    ggsidekick::theme_sleek() +
    labs(colour = "Cutoff") +
    ggtitle(stringr::str_to_title(dat$species[1]))
  # geom_hline(yintercept = 0, lty = 2, col = "grey50")
}

d2 <- filter(out_df, species != "lingcod")
d2$species <- as.character(d2$species)
g <- lapply(split(d2, d2$species), make_panel)
library(patchwork)

g[[1]] <- g[[1]] + tagger::tag_facets(tag_prefix = "(") + theme(tagger.panel.tag.text = element_text(colour = "grey50"))
g[[2]] <- g[[2]] + tagger::tag_facets(tag_prefix = "(", tag_pool = letters[3:4]) + theme(tagger.panel.tag.text = element_text(colour = "grey50"))
g[[3]] <- g[[3]] + tagger::tag_facets(tag_prefix = "(", tag_pool = letters[5:6]) + theme(tagger.panel.tag.text = element_text(colour = "grey50"))

patchwork::wrap_plots(g, ncol = 1, axes = "collect", guides = "collect")
ggsave("figures/groundfish-cv-density.pdf", width = 7, height = 7)
ggsave("figures/groundfish-cv-density.png", width = 7, height = 7)

theta <- readRDS("output/gf-cv-full-fit.rds")
theta <- bind_rows(theta) |>
  filter(species != "yelloweye rockfish") |>
  filter(!isFALSE(converged))

theta |>
  tidyr::pivot_longer(cols = estimate) |>
  filter(value < 1000) |>
  ggplot(aes(n, value)) +
  geom_line() +
  facet_grid(term ~ species, scales = "free")

est <- theta |> tidyr::pivot_longer(cols = estimate, values_to = "est")
lwr <- theta |> tidyr::pivot_longer(cols = conf.low, values_to = "lwr")
upr <- theta |> tidyr::pivot_longer(cols = conf.high, values_to = "upr")
est <- dplyr::bind_cols(est, select(lwr, lwr))
est <- dplyr::bind_cols(est, select(upr, upr))

est |>
  filter(sanity) |>
  # filter(upr < 1000) |>
  filter(cutoff < 80) |>
  filter(term %in% c("phi", "sigma_O", "range_a", "range_b")) |>
  # filter(term %in% c("range_a")) |>
  ggplot(aes(n, est)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey90") +
  geom_line() +
  facet_grid(term ~ stringr::str_to_title(species), scales = "free") +
  ggsidekick::theme_sleek() +
  ylab("Estimate") +
  xlab("Mesh vertices")
ggsave("figures/groundfish-cv-parameters.pdf", width = 8, height = 6)
ggsave("figures/groundfish-cv-parameters.png", width = 8, height = 6)
