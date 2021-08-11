# remotes::install_github("inlabru-org/fmesher", ref = "stable")
# library(sf)
library(inlabru)
library(INLA)
library(dplyr)
# library(fmesher)
# https://www.maths.ed.ac.uk/~flindgre/2018/07/22/spatially-varying-mesh-quality/

haul <- readRDS("data/haul_cleaned.rds")
catch <- readRDS("data/catch_cleaned.rds")

# species are dover sole, sablefish, petrale sole, lingcod, darkblotched rockfish
# g = group_by(catch, common_name) %>% dplyr::summarise(p = sum(present)) %>% dplyr::arrange(-p)

# initial loop over the cutoff values
df <- expand.grid(
  "cutoff" = seq(3, 120, by = 6),
  "species" = c(
    "Dover sole", "petrale sole", "darkblotched rockfish",
    "lingcod"),
  "n" = NA,
  "b"=NA,
  "b_sd"=NA,
  "b_lo"=NA,
  "b_hi"=NA,
  "range" = NA,
  "range_sd"=NA,
  "range_lo"=NA,
  "range_hi"=NA,
  "sd" = NA,
  "sd_sd"=NA,
  "sd_lo"=NA,
  "sd_hi"=NA)

set.seed(2021)

for (i in 1:nrow(df)) {
  catch_sub <- dplyr::filter(catch, common_name == df$species[i])

  # Join catch and haul data
  haul_new <- haul %>%
    left_join(catch_sub, by = "trawl_id") %>%
    select(trawl_id, X, Y,
      latitude = latitude_dd.x,
      longitude = longitude_dd.x,
      year = year,
      log_depth_scaled,
      log_depth_scaled2,
      cpue_kg_km2
    )
  # Set NA CPUEs to 0
  haul_new$cpue_kg_km2[which(is.na(haul_new$cpue_kg_km2))] <- 0

  # convert coordinates to km
  haul_new$X <- haul_new$X / 1000
  haul_new$Y <- haul_new$Y / 1000
  # create occurrence field
  haul_new$present <- ifelse(haul_new$cpue_kg_km2 > 0, 1, 0)
  coordinates(haul_new) <- c("X", "Y")

  # create boundary, same for all meshes
  boundary <- inla.nonconvex.hull(coordinates(haul_new),
    convex = -0.05
  )

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
  components <- present ~ Intercept + field(
    main = coordinates,
    model = matern
  )

  fit <- try(bru(components,
                       haul_new,
                       family = "binomial"
  ), silent = TRUE)
  if(fit$ok==TRUE) {
    df[i,c("range","range_sd","range_lo","range_hi")] <- 
      as.numeric(fit$summary.hyperpar[1,c("mean","sd","0.025quant","0.975quant")])
    df[i,c("sd","sd_sd","sd_lo","sd_hi")] <- 
      as.numeric(fit$summary.hyperpar[2,c("mean","sd","0.025quant","0.975quant")])
    df[i,c("b","b_sd","b_lo","b_hi")] <- 
      as.numeric(fit$summary.fixed[1,c("mean","sd","0.025quant","0.975quant")])
  }
  
  saveRDS(df, "output/09_binom_params_4species.rds")
}
