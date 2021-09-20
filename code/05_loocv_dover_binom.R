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

dover <- readRDS("data/doversole_cleaned.rds")
# convert coordinates to km
dover$X <- dover$X / 1000
dover$Y <- dover$Y / 1000
# create occurrence field
dover$present <- ifelse(dover$cpue_kg_km2 > 0, 1, 0)

set.seed(2021)

coordinates(dover) <- c("X", "Y")

# initial loop over the cutoff values
df <- expand.grid(
  "cutoff" = c(10, 15, 20, 25, 30, 50, 75),
  "holdout" = seq(1, nrow(dover)),
  "n" = NA,
  "dens_ll" = NA
)
# create boundary, same for all meshes
boundary <- inla.nonconvex.hull(coordinates(dover),
  convex = -0.05
)

for (i in 1:nrow(df)) {
  print(i)
  # create mesh
  mesh <- inla.mesh.2d(
    loc = coordinates(dover),
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
  # components is equivalent to formula
  components <- present ~ Intercept + field(
    main = coordinates,
    model = matern
  )

  # do cross validation here
  fit_train <- try(bru(components,
    dover[-df$holdout[i], , drop = FALSE],
    family = "binomial"
  ), silent = TRUE)
  pred_test <- predict(fit_train, dover[df$holdout[i], , drop = FALSE], ~ Intercept + field)

  # calculate total log density
  df$dens_ll[i] <-
    dbinom(
      dover$present[df$holdout[i]],
      size = 1,
      prob = plogis(pred_test$mean),
      log = TRUE
    )

  saveRDS(df, "output/05_binomial_loocv.rds")
}
