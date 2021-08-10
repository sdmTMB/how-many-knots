#remotes::install_github("inlabru-org/fmesher", ref = "stable")
library(sf)
library(inlabru)
library(INLA)
library(dplyr)
#library(fmesher)
#https://www.maths.ed.ac.uk/~flindgre/2018/07/22/spatially-varying-mesh-quality/
  
dover <- readRDS("data/doversole_cleaned.rds")
# convert coordinates to km
dover$X <- dover$X / 1000
dover$Y <- dover$Y / 1000
# create occurrence field
dover$present <- ifelse(dover$cpue_kg_km2 > 0, 1, 0)

set.seed(2021)
n_folds <- 10
n_blocks <- 10
# first assign blocks
dover$block <- 1
for (jj in 1:n_blocks) {
  dover$block[which(dover$latitude < quantile(dover$latitude, 1 - jj * (1 / n_blocks)))] = jj + 1
}
# now assign folds
block_fold <- data.frame("block" = 1:n_blocks,
                        "fold" = rep(1:n_folds, n_blocks / n_folds))
dover <- dplyr::left_join(dover, block_fold)


coordinates(dover) <- c("X", "Y")

# initial loop over the cutoff values
df <- data.frame("cutoff" = seq(3, 120, by=6),
                 "n" = NA,
                "dens_ll" = NA)
# create boundary, same for all meshes
boundary <- inla.nonconvex.hull(coordinates(dover), 
                                convex = -0.05)

for (i in 1:nrow(df)) {
  # create mesh
  mesh <- inla.mesh.2d(loc = coordinates(dover),
                       boundary = boundary,
                       offset = c(-0.05,-0.05),
                       max.n = 5000,
                       max.n.strict = 5000,
                       cutoff = df$cutoff[i],
                       max.edge = c(100,500))
  df$n[i] = mesh$n

  # use PC prior for matern model
  matern <-
    inla.spde2.pcmatern(mesh,
                        prior.sigma = c(10, 0.01),
                        prior.range = c(1, 0.01))
  # components is equivalent to formula
  components <- present ~ Intercept + field(map = coordinates,
                                            model = matern)
  
  # do cross validation here
  dover$pred <- NA
  for (k in 1:max(dover$fold)) {
    fit_train <- bru(components,
                     dover[which(dover$fold != k),],
                     family = "binomial")
    test_indx <- which(dover$fold == k)
    pred_test <- predict(fit_train, dover[test_indx,])
    dover$pred[test_indx] <- pred_test$Intercept$mean + pred_test$field$mean
  }
  # calculate total log density
  df$dens_ll[i] <-
    sum(dbinom(
      dover$present,
      size = 1,
      prob = plogis(dover$pred),
      log = TRUE
    ))
  saveRDS(df, "output/02_binomial_df.rds")
}