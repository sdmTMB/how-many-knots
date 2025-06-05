# remotes::install_github("inlabru-org/fmesher", ref = "stable")
library(sf)
library(sp)
library(mgcv)
library(dplyr)
library(ggplot2)
library(dplyr)
library(viridis)
library(dplyr)
library(tidyr)
library(purrr)

theme_set(ggsidekick::theme_sleek() +
            theme(legend.position = "bottom"))

ggsave2 <- function(filename, ...) {
  ggsave(paste0(filename, ".png"), ...)
  ggsave(paste0(filename, ".pdf"), ...)
}
set.seed(2021)

haul <- readRDS("data/haul_cleaned.rds")
haul$date_formatted <- as.Date(as.numeric(haul$date_formatted), origin = "1970-01-01")
haul$year <- lubridate::year(haul$date_formatted)
haul$month <- lubridate::month(haul$date_formatted)
haul$yday <- lubridate::yday(haul$date_formatted)
haul$zday <- scale(haul$yday)
haul_new <- dplyr::filter(haul, year==2018,
                          !is.na(temperature_at_gear_c_der))
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
  "df_ratio" = seq(0.5, 1.8, by = 0.1),
  "n" = 0,
  "ll_train" = 0,
  "ll_test" = 0,
  "rmse_test" = NA,
  "rmse_train" = NA,
  "converged" = NA,
  "EDF" = NA
)


for (i in 1:nrow(df)) {
  print(i)
  
  # do cross validation here
  for(ii in 1:10) { # tensor, because variability may not be the same in X/Y directions
    fit_train <- gam(temperature_at_gear_c_der ~ zday + I(zday^2) + 
                       s(X, k = df$df_ratio[i]*10) + 
                       s(Y, k = df$df_ratio[i]*5) + 
                       ti(X, Y, k = df$df_ratio[i]*15),
                     data = haul[-which(haul$fold_id==ii), , drop = FALSE])
    pred_test <- predict(
      fit_train, haul[which(haul$fold_id==ii), , drop = FALSE],
    )
    pred_train <- predict(
      fit_train, haul[which(haul$fold_id!=ii), , drop = FALSE],
    )
    
    # df
    df$n[i] <- df$n[i] + sum(summary(fit_train)$s.table[, "edf"])
    
    # calculate total log density
    df$ll_train[i] <- df$ll_train[i] + 
      sum(dnorm(haul$temperature_at_gear_c_der[which(haul$fold_id!=ii)],
        mean = pred_train,
        sd = sqrt(fit_train$sig2),
        log = TRUE
      ))
    
    df$ll_test[i] <- df$ll_test[i] + 
      sum(dnorm(haul$temperature_at_gear_c_der[which(haul$fold_id==ii)],
                mean = pred_test,
                sd = sqrt(fit_train$sig2),
                log = TRUE
      ))
    
  }
  saveRDS(df, "output/01_loocv_temp_mgcv.rds")
}


df_train <- dplyr::select(df, n, ll_train) |>
  dplyr::rename(ll = ll_train) |>
  dplyr::mutate(type = "Train")
# each dataset is used 9 times for training
df_train$ll <- df_train$ll / 9 

df_test <- dplyr::select(df, n, ll_test) |>
  dplyr::rename(ll = ll_test) |>
  dplyr::mutate(type = "Test")

# divide n by 10, because it's accumulated for each fold above and we want average across folds
ggplot(rbind(df_train, df_test), aes(n/10, ll, col = type)) + 
  geom_line() + 
  xlab("Average EDF") + 
  ylab("Log-likelihood")
ggsave2("figures/mgcv_temperature", height = 3, width = 6)