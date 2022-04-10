# remotes::install_github("inlabru-org/fmesher", ref = "stable")
d
# species are dover sole, sablefish, petrale sole, lingcod, darkblotched rockfish
# g = group_by(catch, common_name) %>% dplyr::summarise(p = sum(present)) %>% dplyr::arrange(-p)

# initial loop over the cutoff values
df <- expand.grid(
  "cutoff" = seq(10, 25),
  "species" = c(
    "Dover sole", "petrale sole", "darkblotched rockfish",
    "lingcod"
  ),
  "n" = NA,
  "mean" = NA, "sd" = NA, "lo" = NA, "hi" = NA
)

set.seed(2021)

for (i in 1:nrow(df)) {
  print(i)
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
  haul_new <- dplyr::filter(haul_new, cpue_kg_km2 > 0)
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
  components <- cpue_kg_km2 ~ Intercept + field(
    main = coordinates,
    model = matern
  )

  fit_train <- try(bru(components,
    haul_new,
    family = "Gamma"
  ), silent = TRUE)

  # calculate total log density
  df$mean[i] <- fit_train$summary.fixed$mean
  df$sd[i] <- fit_train$summary.fixed$sd
  df$lo[i] <- fit_train$summary.fixed$`0.025quant`
  df$hi[i] <- fit_train$summary.fixed$`0.975quant`
}

saveRDS(df, "output/df_gamma_catchrates_4sp.rds")
