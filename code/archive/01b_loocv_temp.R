# remotes::install_github("inlabru-org/fmesher", ref = "stable")
library(sf)
library(inlabru)
library(INLA)
library(dplyr)
# library(fmesher)
# https://www.maths.ed.ac.uk/~flindgre/2018/07/22/spatially-varying-mesh-quality/
library(blockCV)
library(future)
plan(multisession)

haul <- readRDS("data/june_bottom_temp.rds")
#haul[, c("X", "Y")] <- haul[, c("X", "Y")] / 1000
set.seed(2021)

coordinates(haul) <- c("X", "Y")

# initial loop over the cutoff values
df <- expand.grid(
  "cutoff" = seq(20, 160, by = 2.5),
  "holdout" = seq(1, nrow(haul)),
  "n" = NA,
  "dens_ll" = NA
)
# create boundary, same for all meshes
boundary <- inla.nonconvex.hull(coordinates(haul),
  convex = -0.05
)

for (i in 1:nrow(df)) {
  print(i)
  # create mesh
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

  # use PC prior for matern model
  matern <-
    inla.spde2.pcmatern(mesh,
      prior.sigma = c(5, 0.05),
      prior.range = c(20, 0.05)
    )
  components <- temperature_at_gear_c_der ~ Intercept(1) + log_depth_scaled + log_depth_scaled2 +
    field(
      main = coordinates,
      model = matern
    )

  # do cross validation here
  fit_train <- try(bru(components,
    haul[-df$holdout[i], , drop = FALSE],
    family = "gaussian"
  ), silent = TRUE)
  pred_test <- predict(
    fit_train, haul[df$holdout[i], , drop = FALSE],
    ~ Intercept + log_depth_scaled + log_depth_scaled2 + field
  )

  # calculate total log density
  tau <- fit_train$summary.hyperpar$mean[1]

  df$dens_ll[i] <-
    dnorm(haul$temperature_at_gear_c_der[df$holdout[i]],
      mean = pred_test$mean,
      sd = sqrt(1 / tau),
      log = TRUE
    )

  saveRDS(df, "output/01_loocv_temp.rds")
}
