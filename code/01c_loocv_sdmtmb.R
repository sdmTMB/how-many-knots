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
library(sdmTMB)
set.seed(1)
N = 1000
predictor_dat <- data.frame(X = runif(N), Y = runif(N), a1 = rnorm(N))
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.05)

dat <- sdmTMB_simulate(
  formula = ~ 1 + a1,
  data = predictor_dat,
  mesh = mesh,
  family = gaussian(),
  range = 0.1,
  sigma_O = 1,
  seed = 42,
  phi = 0.01,
  B = c(3, -0.1) # B0 = intercept, B1 = a1 slope
)

# initial loop over the cutoff values

cutoffs = c(0.02,0.05,0.1,0.2)
meshes = list()
for(i in 1:length(cutoffs)) {
  meshes[[i]] <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = cutoffs[i])
}

df <- expand.grid(
  "cutoff" = cutoffs,
  "holdout" = seq(1, nrow(dat)),
  "n" = NA,
  "dens_ll" = NA,
  "include_cov" = c(TRUE,FALSE)
)

for (i in 1:nrow(df)) {
  print(i)
  save.image("test.Rdata")
  mesh = meshes[[which(cutoffs==df$cutoff[i])]]
  weights = rep(1, nrow(dat))
  weights[df$holdout[i]] = 0
  
  if(df$include_cov[i]==TRUE) {
  fit_train <- try(sdmTMB(observed ~ 1 + a1, data = dat, 
                family = gaussian(),
                weights = weights,
                mesh = mesh),silent=TRUE)
  } else {
    fit_train <- try(sdmTMB(observed ~ 1, data = dat, 
                            family = gaussian(),
                            weights = weights,
                            mesh = mesh),silent=TRUE) 
  }
  
  if(class(fit_train) != "try-error") {
    # calculate total log density
    pred <- predict(fit_train, newdata = dat[df$holdout[c(i,i)],])
    
    df$dens_ll[i] <-
      dnorm(dat$observed[df$holdout[i]],
        mean = pred$est[1],
        sd = exp(fit_train$sd_report$par.fixed[["ln_phi"]]),
        log = TRUE
      )
  }
  saveRDS(df, "output/01_loocv_sim.rds")
}

dplyr::group_by(df, include_cov, cutoff) %>% 
  dplyr::summarize(s = sum(dens_ll))