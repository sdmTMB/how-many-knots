# remotes::install_github("inlabru-org/fmesher", ref = "stable")
library(sf)
library(inlabru)
library(INLA)
library(dplyr)
# library(fmesher)
# https://www.maths.ed.ac.uk/~flindgre/2018/07/22/spatially-varying-mesh-quality/

set.seed(2021)
haul <- readRDS("data/june_bottom_temp.rds")

# initial loop over the cutoff values
df <- expand.grid(
  "cutoff" = seq(30, 160, by = 2),
  "blocks" = 10,
  "n" = NA,
  "log_cpo" = NA,
  "dic" = NA
)

for (i in 1:nrow(df)) {

  # Set NA CPUEs to 0
  # haul_new$cpue_kg_km2[which(is.na(haul_new$cpue_kg_km2))] <- 0
  haul_new <- dplyr::filter(
    haul,
    !is.na(temperature_at_gear_c_der),
  )
  # create occurrence field
  # haul_new$present <- ifelse(haul_new$cpue_kg_km2 > 0, 1, 0)
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
  components <- temperature_at_gear_c_der ~ Intercept(1) + log_depth_scaled + log_depth_scaled2 +
    field(
      main = coordinates,
      model = matern
    )

  # calculate cpo
  fit <- bru(components,
    data = haul_new,
    family = "gaussian",
    options = list(
      verbose = TRUE,
      control.compute = list(config = TRUE, cpo = TRUE, dic = TRUE, openmp.strategy = "huge")
    )
  )
  df$log_cpo[i] <- sum(log(fit$cpo$cpo))
  df$dic[i] <- fit$dic$dic

  saveRDS(df, "output/01_temp_cpo.rds")
}
