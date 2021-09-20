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

catch <- readRDS("data/catch_cleaned_top50.rds")
catch <- dplyr::filter(catch, cpue_kg_km2 > 0)
haul <- readRDS("data/haul_cleaned.rds")

set.seed(2021)

# initial loop over the cutoff values
df <- expand.grid(
  "cutoff" = seq(5, 50, by = 2.5),
  "species" = unique(catch$common_name),
  "range" = 75,
  "folds" = 10,
  "n" = NA,
  "dens_ll" = NA
)


for (i in 1:nrow(df)) {
  catch_sub <- dplyr::filter(catch, common_name == df$species[i])

  # Join catch and haul data
  haul_new <- haul %>%
    left_join(catch_sub, by = "trawl_id") %>%
    select(
      trawl_id,
      X,
      Y,
      latitude = latitude_dd.x,
      longitude = longitude_dd.x,
      year = year,
      log_depth_scaled,
      log_depth_scaled2,
      cpue_kg_km2
    )
  # Set NA CPUEs to 0
  haul_new$cpue_kg_km2[which(is.na(haul_new$cpue_kg_km2))] <- 0
  haul_new <- dplyr::filter(haul_new, cpue_kg_km2 > 0)

  # convert coordinates to km
  haul_new$X <- haul_new$X / 1000
  haul_new$Y <- haul_new$Y / 1000
  # create occurrence field
  # haul_new$present <- ifelse(haul_new$cpue_kg_km2 > 0, 1, 0)
  coordinates(haul_new) <- c("X", "Y")

  # create boundary, same for all meshes
  boundary <- inla.nonconvex.hull(coordinates(haul_new),
    convex = -0.05
  )


  pa_data <-
    sf::st_as_sf(haul_new, coords = c("longitude", "latitude"))
  sb <- try(spatialBlock(
    speciesData = pa_data,
    theRange = df$range[i],
    k = df$folds[i],
    selection = "systematic",
    showBlocks = FALSE
  ),
  silent = TRUE
  )
  if (class(sb) != "try-error") {
    haul_new$fold <- sb$foldID
    df$blocks[i] <- nrow(as.data.frame(sb$blocks))

    # create mesh
    mesh <- inla.mesh.2d(
      loc = coordinates(haul_new),
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
    components <- cpue_kg_km2 ~ Intercept + field(
      main = coordinates,
      model = matern
    )

    # do cross validation here
    haul_new$pred <- NA
    fold_ll <- 0
    for (k in 1:max(haul_new$fold)) {
      fit_train <- try(bru(components,
        haul_new[which(haul_new$fold != k), ],
        family = "gamma"
      ), silent = TRUE)
      test_indx <- which(haul_new$fold == k)
      # if model didn't have problems
      if (class(fit_train)[1] == "bru") {
        if (fit_train$ok) {
          pred_test <- predict(fit_train,
            data = haul_new[test_indx, , drop = FALSE],
            formula = ~ Intercept + field
          )
          haul_new$pred[test_indx] <- pred_test$mean
        } else {
          # model had issues
          haul_new$pred[test_indx] <- NA
        }
      }
      # get precision parameter for the Gamma observations
      gamma_prec <- fit_train$summary.hyperpar$mean[1]
      # evaluate gamma likelihood for this fold
      # gamma_shape = gamma_prec
      # gamma_scale = mu / gamma_prec
      fold_ll[k] <- sum(
        dgamma(
          haul_new$cpue_kg_km2[test_indx],
          shape = gamma_prec,
          scale = exp(haul_new$pred[test_indx]) / gamma_prec,
          log = TRUE
        )
      )
    }
    # calculate total log density
    df$dens_ll[i] <- sum(fold_ll)
    saveRDS(df, "output/50_gamma_cv_df_top50.rds")
  }
}
