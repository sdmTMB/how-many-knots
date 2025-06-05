# remotes::install_github("inlabru-org/fmesher", ref = "stable")
library(sf)
library(sp)
library(inlabru)
library(INLA)
library(dplyr)
library(fmesher)
# https://www.maths.ed.ac.uk/~flindgre/2018/07/22/spatially-varying-mesh-quality/
library(blockCV)
library(future)
library(ggplot2)

theme_set(ggsidekick::theme_sleek() +
            theme(legend.position = "bottom"))

ggsave2 <- function(filename, ...) {
  ggsave(paste0(filename, ".png"), ...)
  ggsave(paste0(filename, ".pdf"), ...)
}
plan(multisession)

set.seed(2021)

haul <- readRDS("data/haul_cleaned.rds")
haul$date_formatted <- as.Date(as.numeric(haul$date_formatted), origin = "1970-01-01")
haul$year <- lubridate::year(haul$date_formatted)
haul$month <- lubridate::month(haul$date_formatted)
haul$yday <- lubridate::yday(haul$date_formatted)
haul$zday <- scale(haul$yday)
haul_new <- dplyr::filter(haul, year==2018,
                          !is.na(temperature_at_gear_c_der))
haul_new$log_depth_scaled <- as.numeric(haul_new$log_depth_scaled)
haul_new$log_depth_scaled2 <- as.numeric(haul_new$log_depth_scaled2)
haul_new$fold_id <- rep(1:10, length.out = nrow(haul_new))
haul_new$pred_train <- NA
haul <- haul_new

ns <- dplyr::group_by(haul, fold_id) |>
  dplyr::summarise(n = n(), n_other = nrow(haul)-n) |>
  dplyr::ungroup() |>
  dplyr::summarise(mean_n = mean(n), mean_other = mean(n_other))

#haul$fold <- sample(c(1,2), size=nrow(haul), replace=T, prob = c(0.1,0.9))
# initial loop over the cutoff values
df <- expand.grid(
  "cutoff" = exp(seq(log(8), log(500), length.out=200)),
  "n" = NA,
  "ll_train" = 0,
  "ll_test" = 0,
  "rmse_test" = NA,
  "rmse_train" = NA,
  "converged" = NA,
  "EDF" = NA
)

coordinates(haul) <- c("X", "Y")

# create boundary, same for all meshes
boundary <- inla.nonconvex.hull(coordinates(haul),
  convex = -0.05
)

for (i in nrow(df):1) {
  print(i)
  
  # Create mesh
  mesh <- inla.mesh.2d(
    loc = coordinates(haul),
    boundary = boundary,
    offset = c(-0.05, -0.05),
    max.n = 5000,
    max.n.strict = 5000,
    cutoff = df$cutoff[i],
    max.edge = c(100, 500)
  )
  df$n[i] <- mesh$n
  
  # Create SPDE model with PC priors
  matern <- inla.spde2.pcmatern(
    mesh,
    prior.sigma = c(5, 0.05),
    prior.range = c(20, 0.05)
  )
  
  components <- temperature_at_gear_c_der ~ log_depth_scaled + I(log_depth_scaled^2) +
    field(main = coordinates, model = matern)
  
  # Set log-likelihood to 0, and fold failure flag to FALSE
  df$ll_train[i] <- 0
  df$ll_test[i] <- 0
  fold_failed <- FALSE
  
  for (ii in 1:10) {
    message("  Fold ", ii)
    
    fit_train <- tryCatch({
      bru(
        components,
        data = haul[haul$fold_id != ii, , drop = FALSE],
        family = "gaussian"
      )
    }, error = function(e) {
      message("    -> bru() failed: ", e$message)
      return(NULL)
    })
    
    if (is.null(fit_train)) {
      fold_failed <- TRUE
      break
    }
    
    pred_train <- tryCatch({
      predict(fit_train,
              newdata = haul[haul$fold_id != ii, , drop = FALSE],
              ~ zday + I(zday^2) + field
      )
    }, error = function(e) {
      message("    -> pred_train failed: ", e$message)
      return(NULL)
    })
    
    pred_test <- tryCatch({
      predict(fit_train,
              newdata = haul[haul$fold_id == ii, , drop = FALSE],
              ~ zday + I(zday^2) + field
      )
    }, error = function(e) {
      message("    -> pred_test failed: ", e$message)
      return(NULL)
    })
    
    tau <- tryCatch({
      fit_train$summary.hyperpar$mean[1]
    }, error = function(e) {
      message("    -> failed to extract tau: ", e$message)
      return(NA)
    })
    
    # If any part failed, mark the fold as failed
    if (is.null(pred_train) || is.null(pred_test) || is.na(tau)) {
      fold_failed <- TRUE
      break
    }
    
    # Otherwise accumulate log-likelihoods
    df$ll_train[i] <- df$ll_train[i] +
      sum(dnorm(
        haul$temperature_at_gear_c_der[haul$fold_id != ii],
        mean = pred_train$mean,
        sd = sqrt(1 / tau),
        log = TRUE
      ))
    
    df$ll_test[i] <- df$ll_test[i] +
      sum(dnorm(
        haul$temperature_at_gear_c_der[haul$fold_id == ii],
        mean = pred_test$mean,
        sd = sqrt(1 / tau),
        log = TRUE
      ))
  }
  
  # If any fold failed, overwrite log-likelihoods as NA
  if (fold_failed) {
    message("  -> One or more folds failed; setting ll_test and ll_train to NA.")
    df$ll_train[i] <- NA
    df$ll_test[i] <- NA
  }
  
  # Save progress
  saveRDS(df, "output/01_loocv_temp.rds")
}

df <- readRDS("output/01_loocv_temp.rds") |>
  dplyr::filter(ll_test != 0, ll_train != 0) |>
  dplyr::filter(n != 341, n < 500) # problematic convergence
  
df_train <- dplyr::select(df, n, ll_train) |>
  dplyr::rename(ll = ll_train) |>
  dplyr::mutate(type = "Train")
# each dataset is used 9 times for training
df_train$ll <- df_train$ll / 9 

df_test <- dplyr::select(df, n, ll_test) |>
  dplyr::rename(ll = ll_test) |>
  dplyr::mutate(type = "Test")

# divide n by 10, because it's accumulated for each fold above and we want average across folds
ggplot(rbind(df_train, df_test), aes(n, ll, col = type)) + 
  geom_point() + 
  xlab("Mesh vertices") + 
  ylab("Log-likelihood")
ggsave2("figures/inlabru_temperature", height = 3, width = 6)
