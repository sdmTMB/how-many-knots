# devtools::install_github("seananderson/ggsidekick")
library(ggsidekick)
library(sf)
library(dplyr)
library(sp)
library(inlabru)
library(fmesher)
library(lubridate)
library(patchwork)
library(ggplot2)
library(viridis)
library(tidyr)
library(raster)
library(scales)
library(cowplot)
# https://www.maths.ed.ac.uk/~flindgre/2018/07/22/spatially-varying-mesh-quality/

haul <- readRDS("data/haul_cleaned.rds")
haul$date_formatted <- as.Date(as.numeric(haul$date_formatted), origin = "1970-01-01")

haul$year <- lubridate::year(haul$date_formatted)
haul <- dplyr::filter(haul, year == 2018)
 
coordinates(haul) <- c("X", "Y")

# create boundary, same for all meshes
boundary <- fmesher::fm_nonconvex_hull_inla(sp::coordinates(haul),
  convex = -0.05
)

mesh1 <- fmesher::fm_mesh_2d_inla(
  loc = coordinates(haul),
  boundary = boundary,
  offset = c(-0.05, -0.05),
  max.n = 5000,
  max.n.strict = 5000,
  cutoff = 20,
  max.edge = c(100, 500)
)

mesh2 <- fmesher::fm_mesh_2d_inla(
  loc = coordinates(haul),
  boundary = boundary,
  offset = c(-0.05, -0.05),
  max.n = 5000,
  max.n.strict = 5000,
  cutoff = 80,
  max.edge = c(100, 500)
)

plot_df <- data.frame(coordinates(haul))

p3 <- ggplot() +
  gg(mesh2) +
  ggsidekick::theme_sleek() +
  xlab("") +
  ylab("") +
  geom_point(data = plot_df, aes(X, Y), alpha = 0.3, col = viridis(1), size = 0.4) + coord_fixed()

p2 <- ggplot() +
  gg(mesh1) +
  ggsidekick::theme_sleek() +
  xlab("") +
  ylab("") +
  geom_point(data = plot_df, aes(X, Y), alpha = 0.3, col = viridis(1), size = 0.4) + coord_fixed()


# Make the map
haul <- readRDS("data/haul_cleaned.rds")
haul$date_formatted <- as.Date(as.numeric(haul$date_formatted), origin = "1970-01-01")
haul$year <- lubridate::year(haul$date_formatted)
haul <- dplyr::filter(
  haul, year == 2018,
  !is.na(temperature_at_gear_c_der)
)
map_data <- rnaturalearth::ne_countries(
  scale = "large",
  returnclass = "sf", country = c("canada", "united states of america", "mexico"))
coast <- suppressWarnings(suppressMessages(
  sf::st_crop(map_data,
              c(xmin = -130, ymin = 30, xmax = -107, ymax = 55))))
coast_proj <- sf::st_transform(coast, crs = 3157) # zone 10

sf::st_bbox(coast_proj)

.xlim <- c(235000, 1015000)
.ylim <- c(3586000, 5560000)
# Map plot with explicit white background
p1 <- ggplot(coast_proj) +
  geom_sf(data = coast_proj, fill = "grey80") +
  geom_point(data = haul, aes(x = X*1000, y = Y*1000, col = temperature_at_gear_c_der), size = 0.3) +
  # scale_color_gradient2(name = "\u00B0C", midpoint = mean(haul$temperature_at_gear_c_der)) +
  scale_color_viridis_c(name = "\u00B0C", option = "G") +
  labs(x = "", y = "") +
  ggsidekick::theme_sleek() +
  coord_sf(xlim = .xlim, ylim = .ylim) +
  annotate(geom = "text", x = max(.xlim) - 310000, y = mean(.ylim) + 100000, label = "United States", size = 3, colour = "grey10", hjust = 0.5, vjust = 0.5) +
  annotate(geom = "text", x = max(.xlim) - 310000, y = mean(.ylim) + 1000000, label = "Canada", size = 3, colour = "grey10", hjust = 0.5, vjust = 0.5) +
  ggspatial::annotation_north_arrow(
    location = "bl", which_north = "true",
    pad_x = unit(0, "in"), pad_y = unit(0.05, "in"),
    height = unit(1.2, "cm"),
    width = unit(1.2, "cm"),
    style = ggspatial::north_arrow_nautical(
      fill = c("grey40", "white"),
      line_col = "grey20"
    )
  ) + 
  scale_x_continuous(breaks = c(-130, -125, -120)) +
  scale_y_continuous(breaks = c(30, 35, 40, 45, 50))
p1



# Convert 
coords_utm <- as.data.frame(mesh1$loc * 1000)
colnames(coords_utm) <- c("X", "Y")
coords_sf <- st_as_sf(coords_utm, coords = c("X", "Y"), crs = 3157)
coords_ll <- st_transform(coords_sf, crs = 4326)
coords_ll_mat <- st_coordinates(coords_ll)
mesh1_ll <- mesh1
mesh1_ll$loc <- coords_ll_mat

coords_utm2 <- as.data.frame(mesh2$loc * 1000)
colnames(coords_utm2) <- c("X", "Y")  # Required for st_as_sf()
coords_sf2 <- sf::st_as_sf(coords_utm2, coords = c("X", "Y"), crs = 3157)
coords_ll2 <- sf::st_transform(coords_sf2, crs = 4326)
coords_ll_mat2 <- sf::st_coordinates(coords_ll2)
mesh2_ll <- mesh2
mesh2_ll$loc <- coords_ll_mat2

# Multiply haul plot_df by 1000 to get meters
plot_df_m <- plot_df * 1000
haul_sf <- st_as_sf(plot_df_m, coords = c("X", "Y"), crs = 3157)
haul_ll <- st_transform(haul_sf, crs = 4326)
plot_df_ll <- as.data.frame(st_coordinates(haul_ll))

p3 <- ggplot() +
  inlabru::gg(mesh2_ll) +
  geom_point(data = plot_df_ll, aes(x = X, y = Y),
             alpha = 0.3, col = viridis::viridis(1), size = 0.4) +
  ggsidekick::theme_sleek() +
  coord_sf(crs = 4326) +
  xlab("") +
  ylab("") + 
  scale_x_continuous(breaks = c(-130, -125, -120)) +
  scale_y_continuous(breaks = c(30, 35, 40, 45, 50))

p2 <- ggplot() +
  inlabru::gg(mesh1_ll) +
  geom_point(data = plot_df_ll, aes(x = X, y = Y),
             alpha = 0.3, col = viridis::viridis(1), size = 0.4) +
  ggsidekick::theme_sleek() +
  coord_sf(crs = 4326) +
  xlab("") +
  ylab("") + 
  scale_x_continuous(breaks = c(-130, -125, -120)) +
  scale_y_continuous(breaks = c(30, 35, 40, 45, 50))

p1 + p2 + p3 + plot_layout(ncol = 3, axes = "collect", axis_titles = "collect") +
 plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
 theme(plot.tag = element_text(size = 10))

ggsave("figures/Figure_01.png", width = 8, height = 6)
ggsave("figures/Figure_01.pdf", width = 8, height = 6)
