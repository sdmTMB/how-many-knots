library(sf)
library(sdmTMB)
library(lubridate)
library(dplyr)
library(future)
library(viridis)
library(ggplot2)
haul <- readRDS("data/haul_cleaned.rds")
catch <- readRDS("data/catch_cleaned.rds")
haul$trawl_id <- as.numeric(haul$trawl_id)

catch_sub <- dplyr::filter(catch, common_name == "arrowtooth flounder")

# Join catch and haul data
joined_dat <- haul |>
  dplyr::filter(trawl_id %in% catch$trawl_id) |>
  left_join(catch_sub[, c("date", "cpue_kg_km2", "depth_m", "trawl_id")],
            by = "trawl_id"
  ) |>
  dplyr::filter(!is.na(cpue_kg_km2)) |>
  dplyr::mutate(present = ifelse(cpue_kg_km2 > 0, 1, 0))
joined_dat$date <- lubridate::ymd(as.POSIXct.Date(joined_dat$date))
joined_dat$year <- year(joined_dat$date)
joined_dat$fyear <- as.factor(joined_dat$year)
joined_dat$yday <- yday(joined_dat$date)


sp::coordinates(joined_dat) <- c("X", "Y")

joined_dat <- as.data.frame(joined_dat)          # instead of sp::coordinates<-
xy <- as.matrix(joined_dat[, c("X", "Y")])

boundary <- fmesher::fm_nonconvex_hull_inla(xy, convex = -0.05)

cutoffs <- c(25, 50, 150)
fits <- list()
pred <- list()
pred_knots <- list()
for (i in seq_along(cutoffs)) {
  inla_mesh <- fmesher::fm_mesh_2d_inla(
    loc = xy,
    boundary = boundary,
    offset = c(-0.05, -0.05),
    max.n = 5000,
    max.n.strict = 5000,
    cutoff = cutoffs[i],
    max.edge = c(100, 500)
  )
  mesh <- make_mesh(joined_dat, c("X", "Y"), mesh = inla_mesh)
  
  fits[[i]] <- sdmTMB(
    cpue_kg_km2 ~ -1 + as.factor(year) + poly(log_depth_scaled, 2),
    data = joined_dat,
    mesh = mesh,
    time = "year",
    spatial = "on",
    spatiotemporal = "off",
    anisotropy = TRUE,
    share_range = TRUE,
    family = tweedie()
  )
  
  pred[[i]] <- predict(fits[[i]]) |> dplyr::filter(year == 2003)
  knot_df <- as.data.frame(mesh$mesh$loc)
  names(knot_df) <- c("X","Y","Z")
  knot_df$year <- 2003
  knot_df$log_depth_scaled <- 0
  # could also use knot_df$omega_s <- fits[[i]]$tmb_obj$env$parList()$omega_s[, 1]
  pred_knots[[i]] <- predict(fits[[i]], newdata = knot_df)
}

# ---- Coastline (from temperature script) ----
map_data <- rnaturalearth::ne_countries(
  scale = "large", returnclass = "sf",
  country = c("canada", "united states of america", "mexico"))
coast <- suppressWarnings(suppressMessages(
  sf::st_crop(map_data, c(xmin = -130, ymin = 30, xmax = -107, ymax = 55))))
coast_proj <- sf::st_transform(coast, crs = 3157)  # UTM zone 10

.xlim <- c(235000, 1015000)
.ylim <- c(3586000, 5560000)

# ---- Plot ----
pred_all <- dplyr::bind_rows(pred, .id = "mesh") |>
  dplyr::mutate(cutoff = factor(paste0("Cutoff = ", cutoffs[as.integer(mesh)], " km"),
                                levels = paste0("Cutoff = ", cutoffs, " km")))

first_panel <- levels(pred_all$cutoff)[1]
lab_df <- data.frame(
  cutoff = factor(first_panel, levels = levels(pred_all$cutoff)),
  x = max(.xlim) - 310000,
  y = mean(.ylim) + c(100000, 1000000),
  label = c("United States", "Canada")
)
arrow_df <- data.frame(cutoff = factor(first_panel, levels = levels(pred_all$cutoff)))


p_omega <- ggplot(coast_proj) +
  geom_sf(fill = "grey80") +
  geom_point(data = pred_all,
             aes(x = X * 1000, y = Y * 1000, col = omega_s), size = 0.3) +
  scale_color_gradient2(name = expression(omega[s]),
                        low = "#2166AC", mid = "grey95", high = "#B2182B",
                        midpoint = 0) +
  geom_text(data = lab_df, aes(x = x, y = y, label = label),
            size = 3, colour = "grey10", inherit.aes = FALSE) +
  ggspatial::annotation_north_arrow(
    data = arrow_df,
    location = "bl", which_north = "true",
    pad_x = unit(0, "in"), pad_y = unit(0.05, "in"),
    height = unit(1.2, "cm"), width = unit(1.2, "cm"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"), line_col = "grey20")
  ) +
  facet_wrap(~cutoff, nrow = 1) +
  coord_sf(xlim = .xlim, ylim = .ylim) +
  scale_x_continuous(breaks = c(-126, -122, -118)) +
  labs(x = "", y = "")

p_omega

ggsave("figures/arrowtooth_omega_s.png", p_omega, width = 8, height = 5.5, bg = "white")

# ---- Plot ----
pred_all <- dplyr::bind_rows(pred_knots, .id = "mesh") |>
  dplyr::mutate(cutoff = factor(paste0("Cutoff = ", cutoffs[as.integer(mesh)], " km"),
                                levels = paste0("Cutoff = ", cutoffs, " km")))

land <- sf::st_union(coast_proj)

knots_sf <- pred_all |>
  dplyr::mutate(X_m = X * 1000, Y_m = Y * 1000) |>
  sf::st_as_sf(coords = c("X_m", "Y_m"), crs = sf::st_crs(coast_proj))

on_land <- lengths(sf::st_intersects(knots_sf, land)) > 0
pred_all <- pred_all[!on_land, ]

# ---- Labels / arrow on first panel ----
first_panel <- levels(pred_all$cutoff)[1]
lab_df <- data.frame(
  cutoff = factor(first_panel, levels = levels(pred_all$cutoff)),
  x = max(.xlim) - 310000,
  y = mean(.ylim) + c(100000, 1000000),
  label = c("United States", "Canada")
)
arrow_df <- data.frame(cutoff = factor(first_panel, levels = levels(pred_all$cutoff)))

p_omega <- ggplot(coast_proj) +
  geom_sf(fill = "grey80") +
  geom_point(data = pred_all,
             aes(x = X * 1000, y = Y * 1000, col = omega_s), size = 0.3) +
  scale_color_gradient2(name = expression(omega[s]),
                        low = "#2166AC", mid = "grey95", high = "#B2182B",
                        midpoint = 0) +
  geom_text(data = lab_df, aes(x = x, y = y, label = label),
            size = 3, colour = "grey10", inherit.aes = FALSE) +
  ggspatial::annotation_north_arrow(
    data = arrow_df,
    location = "bl", which_north = "true",
    pad_x = unit(0, "in"), pad_y = unit(0.05, "in"),
    height = unit(1.2, "cm"), width = unit(1.2, "cm"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"), line_col = "grey20")
  ) +
  facet_wrap(~cutoff, nrow = 1) +
  coord_sf(xlim = .xlim, ylim = .ylim) +
  scale_x_continuous(breaks = c(-124, -120)) +
  labs(x = "", y = "") +
  theme(panel.spacing.x = unit(1, "lines"))

p_omega

ggsave("figures/arrowtooth_omega_s_knots.png", p_omega, width = 8, height = 5.5, bg = "white")



