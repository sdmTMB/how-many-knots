# remotes::install_github("inlabru-org/fmesher", ref = "stable")
library(sf)
library(inlabru)
library(INLA)
library(dplyr)
# library(fmesher)
# https://www.maths.ed.ac.uk/~flindgre/2018/07/22/spatially-varying-mesh-quality/
library(future)
library(geoR)
plan(multisession)

# initial loop over the cutoff values
df <- expand.grid(
  "cutoff" = c(10,20,30,100,300,500),
  "folds" = 10,
  "range" = seq(25, 175, by = 50),
  "dens_ll" = NA
)

set.seed(2021)

for (i in 20:nrow(df)) {

  data(parana)
  # from inla book, https://www.paulamoraga.com/book-geospatial/sec-geostatisticaldatatheory.html#spatial-modeling-of-rainfall-in-paran%C3%A1-brazil
  parana = as.data.frame(parana)
  coo <- as.matrix(parana[,c("east","north")])
  bnd <- inla.nonconvex.hull(coo)
  # create custom mesh, given the cutoff sequence above
  meshb <- inla.mesh.2d(
    boundary = bnd, offset = c(50, 100),
    cutoff = df$cutoff[i], max.edge = c(30, 60)
  )
  
  matern <- inla.spde2.pcmatern(meshb, prior.sigma=c(5,0.05),
                                prior.range=c(10,0.05))
  components <- data ~ field(main=coordinates,model=matern)
  
  # add coordinates for dataframe
  parana <- as.data.frame(parana)
  
  # add folds
  pa_data <- sf::st_as_sf(parana, coords = c("east","north"))
  sb <- spatialBlock(
    speciesData = pa_data,
    species = "data",
    theRange = df$range[i],
    k = 10,
    selection = "random",
    showBlocks = FALSE
  )
  parana$fold <- sb$foldID
  
  coordinates(parana) = c("east","north")

  # do cross validation here
  parana$pred <- NA
  fold_ll <- 0
  for (k in 1:max(parana$fold)) {
    fit_train <- try(bru(components,
                         parana[which(parana$fold != k), ],
      family = "gaussian"
    ), silent = TRUE)
    test_indx <- which(parana$fold == k)
    # if model didn't have problems
    if (class(fit_train)[1] == "bru") {
      if (fit_train$ok) {
        pred_test <- predict(fit_train, 
          data = parana[test_indx, , drop = FALSE], 
          formula = ~ Intercept + field
        )
        parana$pred[test_indx] <- pred_test$mean
      } else {
        # model had issues
        parana$pred[test_indx] <- NA
      }
    }
    # get precision parameter for the Gamma observations
    tau <- fit_train$summary.hyperpar$mean[1]
    
    fold_ll[k] <- sum(dnorm(parana$data[test_indx],
      mean = parana$pred[test_indx],
      sd = sqrt(1/tau),
      log = TRUE
    ))
  }
  # calculate total log density
  df$dens_ll[i] <- sum(fold_ll)
  saveRDS(df, file = "output/14_parana.rds")
}

pdf("output/archive/parana.pdf")
ggplot(df, aes(cutoff, dens_ll)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~ range) + 
  xlab("Cutoff") + 
  ylab("10-fold predictive log density") + 
  theme_bw()
dev.off()

