# remotes::install_github("inlabru-org/fmesher", ref = "stable")
library(sf)
# library(sp)
library(sdmTMB)
library(lubridate)
library(dplyr)
library(future)
library(viridis)
# plan(multisession)
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
      offset = "log_area_swept_ha",
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
  species = c("sablefish", "arrowtooth flounder", "petrale sole", "yelloweye rockfish")
)
nrow(df)
plan(multicore, workers = 80L)
out <- furrr::future_pmap(df, run_cv, parallel = FALSE)
plan(sequential)
saveRDS(out, file = "output/gf-cv-random-out.rds")

# strip width blocked CV: ---------------------------------------------------------
df <- expand.grid(
  cutoff = round(exp(seq(log(10), log(175), length.out = 25))),
  bin_width = seq(10, 130, by = 40),
  seed = 123,
  species = c("sablefish", "arrowtooth flounder", "petrale sole", "yelloweye rockfish")
)
nrow(df)
plan(multicore, workers = 80L)
out <- furrr::future_pmap(df, run_cv, parallel = FALSE)
saveRDS(out, file = "output/gf-cv-blocked-out.rds")
out2 <- furrr::future_pmap(df[1, ], run_cv, do_full_fit = TRUE, parallel = FALSE)
saveRDS(out2, file = "output/gf-cv-full-fit.rds")
plan(sequential)


if (FALSE) {
  out <- readRDS("output/gf-cv-out.rds")
  out_df <- bind_rows(out)

  out2 <- readRDS("output/gf-cv-out-arrowtooth.rds")
  out_df <- bind_rows(out_df, bind_rows(out2))

  out2 <- readRDS("output/gf-cv-out-petrale.rds")
  out_df <- bind_rows(out_df, bind_rows(out2))

  make_panel <- function(dat) {
    if (dat$species[[1]] == "lingcod") {
      dat <- filter(dat, test_dens_ll_sum > -100000)
    }

    dat |>
      select(n, converged, cutoff, species, test_dens_ll_sum, train_dens_ll_sum) |>
      tidyr::pivot_longer(cols = c(test_dens_ll_sum, train_dens_ll_sum), names_to = "ll_type") |>
      filter(n > 80, converged) |>
      mutate(ll_type = ifelse(grepl("test", ll_type), "Out of sample", "In sample")) |>
      group_by(n, cutoff, ll_type, species) |>
      # summarize(est = mean(exp(test_dens_ll_mean)), lwr = min(mean(exp(test_dens_ll_mean))), upr = max(exp(test_dens_ll_mean))) |>
      summarize(est = mean(value), lwr = min(value), upr = max(value)) |>
      group_by(ll_type, species) |>
      mutate(
        lwr = lwr - max(est),
        upr = upr - max(est),
        est = est - max(est)
      ) |>
      ggplot(aes(n, est, colour = cutoff)) +
      scale_colour_viridis(direction = -1, end = 0.9, limits = c(min(out_df$cutoff), max(out_df$cutoff))) +
      geom_line() +
      facet_wrap(~ll_type, scales = "free_y", ncol = 2) +
      # geom_point(aes(size = cutoff), pch = 21) +
      geom_point(pch = 21) +
      geom_linerange(aes(ymin = lwr, ymax = upr)) +
      geom_smooth(se = FALSE, method = "gam", colour = "grey50") +
      ylab("Relative log likelihood") +
      xlab("Mesh vertices") +
      ggsidekick::theme_sleek() +
      labs(colour = "Cutoff") +
      ggtitle(stringr::str_to_title(dat$species[1]))
  }

  d2 <- filter(out_df, species != "lingcod")
  d2$species <- as.character(d2$species)
  g <- lapply(split(d2, d2$species), make_panel)
  patchwork::wrap_plots(g, ncol = 1, axes = "collect", guides = "collect")

  1
}

##############

#
#
#
#
#
# df <- expand.grid(
#   # cutoff = 40,
#   # bin_width = 10,
#   cutoff = round(exp(seq(log(12), log(175), length.out = 20))),
#   bin_width = seq(10, 130, by = 40),
#   species = c("sablefish", "arrowtooth flounder", "lingcod", "petrale sole")
#   # create mesh
# )
# # out <- purrr::pmap(df, run_cv, parallel = T)
# nrow(df)
# cores <- 40
# nrow(df) * 560 / 60 / 60 / cores
# df
#
# unique(catch$common_name)
# set.seed(1234)
# plan(multicore, workers = 50L)
# tictoc::tic()
# # out <- purrr::pmap(df[1:3, ], run_cv, parallel = T)
# out <- furrr::future_pmap(df[1:3, ], run_cv, parallel = FALSE)
# tictoc::toc()
# saveRDS(out, "output/02-binomial-blockCV-2025-03-31.rds")
#
# tictoc::tic()
# # out2 <- purrr::pmap(df[1:2, ], run_cv, do_full_fit = TRUE, parallel = T)
# out2 <- filter(df, bin_width == 10) |>
#   furrr::future_pmap(run_cv, do_full_fit = TRUE, parallel = FALSE)
# tictoc::toc()
# plan(sequential)
# saveRDS(out2, "output/02-binomial-blockCV-pars-2025-03-28.rds")
#
# out <- readRDS("output/02-binomial-blockCV-2025-03-31.rds")
# out2 <- readRDS("output/02-binomial-blockCV-pars-2025-03-28.rds")
# out2 <- lapply(out2, function(x) {
#  if (!"loglik" %in% names(x)) {
#     x$loglik <- NA_real_
#   }
#   x
# })
# out2 <- lapply(out2, function(x) {
#   x$loglik <- as.numeric(x$loglik)
#   x
# })
# out <- dplyr::bind_rows(out)
# out2 <- dplyr::bind_rows(out2)
#
# head(out)
# head(out)
# head(out2)
# table(out$species)
# names(out)
#
# library(ggplot2)
# out |>
#   filter(present_dens_ll > -1000) |>
#   filter(bin_width < 100) |>
#   ggplot(aes(n, present_dens_ll)) +
#   geom_line() +
#   facet_grid(species~bin_width, scales = "free_y") +
#   geom_smooth(se = FALSE)
#
# unique(out2$species)
# tidyr::pivot_longer(out2, cols = estimate) |>
#   filter(value < 1000) |>
#   ggplot(aes(n, value)) + geom_line() +
#   facet_grid(term~species, scales = "free")
#
# unique(out2$species)
# tidyr::pivot_longer(out2, cols = loglik) |>
#   # filter(value < 1000) |>
#   ggplot(aes(n, value)) + geom_line() +
#   facet_wrap(~species, scales = "free")
#
# unique(out2$species)
# tidyr::pivot_longer(out2, cols = cAIC) |>
#   # filter(value < 1000) |>
#   ggplot(aes(n, value)) + geom_line() +
#   facet_wrap(~species, scales = "free")
#
# tidyr::pivot_longer(out2, cols = edf_omega) |>
#   # filter(value < 1000) |>
#   ggplot(aes(n, value)) + geom_line() +
#   facet_wrap(~species, scales = "free") +
#   geom_abline(intercept = 0, slope = 1, lty = 2)
#
# 1
#
# ## library(ggplot2)
# ## d <- readRDS("output/02_binomial_blockCV.rds")
# ##
# ## # filter out the number of species that don't converge enough
# ## d <- dplyr::filter(d, present_converged == TRUE) |>
# ##   dplyr::group_by(species) |>
# ##   dplyr::mutate(nobs = n()) |>
# ##   dplyr::filter(nobs >= 10) |>
# ##   dplyr::select(-nobs)
# ## d$species <- as.factor(as.character(d$species))
# ##
# ## # Bring in the random
# ## d_random <- readRDS("output/02_binomial_randomCV.rds")
# ## d_random <- dplyr::filter(d_random, present_converged == TRUE) |>
# ##   dplyr::filter(species %in% d$species) |>
# ##   dplyr::mutate(bin_width = NA)
# ##
# ## # Larger bins result in more widely spaced test regions, beyond the estimated
# ## # range. These regions are no longer correlated with the training data and predictions
# ## # become more uncertain / not as good
# ## dplyr::filter(d) |>
# ##   ggplot(aes(n, present_dens_ll, group = bin_width, col = bin_width)) +
# ##   geom_line() +
# ##   facet_wrap(~species, scale = "free_y")
# ##
# ##
# ## d |>
# ##   ggplot(aes(n, range, group = bin_width, col = bin_width)) +
# ##   geom_line() +
# ##   facet_wrap(~species, scale = "free_y") +
# ##   xlab("Mesh vertices") +
# ##   ylab("Estimated spatial range (km)") +
# ##   scale_color_viridis(option = "magma", begin = 0.2, end = 0.8, name = "Strip width (km)") +
# ##   theme_bw() +
# ##   theme(
# ##     strip.background = element_rect(fill = "white"),
# ##     strip.text = element_text(size = 5),
# ##     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
# ##   ) +
# ##   geom_point(data = d_random, aes(n, range), col = "black", alpha = 0.5)
# ## ggsave("figures/groundfish_range_v_n.png", height = 5, width = 7)
# ##
# ##
# ## library(scales)
# ## d |> # divide by 12670 to get average ll per obs
# ##   ggplot(aes(n, present_dens_ll / 12670, group = bin_width, col = bin_width)) +
# ##   geom_line() +
# ##   facet_wrap(~species, scale = "free_y") +
# ##   xlab("Mesh vertices") +
# ##   ylab("Predicted log likelihood") +
# ##   scale_color_viridis(option = "magma", begin = 0.2, end = 0.8, name = "Strip width (km)") +
# ##   theme_bw() +
# ##   theme(
# ##     strip.background = element_rect(fill = "white"),
# ##     strip.text = element_text(size = 5),
# ##     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
# ##     axis.text.y = element_text(size = 8)
# ##   ) +
# ##   geom_point(data = d_random, aes(n, present_dens_ll / 12670), col = "black", alpha = 0.5) # +
# ## # scale_y_continuous(labels = function(x) format(x, scientific=TRUE))
# ## ggsave("figures/groundfish_loglik_v_n.png", height = 5, width = 7)
