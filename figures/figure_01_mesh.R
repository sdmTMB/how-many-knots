# devtools::install_github("seananderson/ggsidekick")
library(ggsidekick)
library(sf)
library(inlabru)
library(INLA)
library(dplyr)
library(sp)
# library(fmesher)
library(lubridate)
library(patchwork)
library(ggplot2)
library(viridis)
# https://www.maths.ed.ac.uk/~flindgre/2018/07/22/spatially-varying-mesh-quality/

haul <- readRDS("data/haul_cleaned.rds")
haul$date_formatted <- as.Date(as.numeric(haul$date_formatted), origin = "1970-01-01")

haul$year <- lubridate::year(haul$date_formatted)
haul <- dplyr::filter(haul, year == 2018)
 
coordinates(haul) <- c("X", "Y")

# create boundary, same for all meshes
boundary <- inla.nonconvex.hull(sp::coordinates(haul),
  convex = -0.05
)

mesh1 <- inla.mesh.2d(
  loc = coordinates(haul),
  boundary = boundary,
  offset = c(-0.05, -0.05),
  max.n = 5000,
  max.n.strict = 5000,
  cutoff = 20,
  max.edge = c(100, 500)
)

mesh2 <- inla.mesh.2d(
  loc = coordinates(haul),
  boundary = boundary,
  offset = c(-0.05, -0.05),
  max.n = 5000,
  max.n.strict = 5000,
  cutoff = 80,
  max.edge = c(100, 500)
)

plot_df <- data.frame(coordinates(haul))

p2 <- ggplot() +
  gg(mesh2) +
  ggsidekick::theme_sleek() +
  xlab("Eastings") +
  ylab("Northings") +
  geom_point(data = plot_df, aes(X, Y), alpha = 0.3, col = viridis(1), size = 0.4) + coord_fixed()

p1 <- ggplot() +
  gg(mesh1) +
  ggsidekick::theme_sleek() +
  xlab("Eastings") +
  ylab("Northings") +
  geom_point(data = plot_df, aes(X, Y), alpha = 0.3, col = viridis(1), size = 0.4) + coord_fixed()

p1 + p2 + plot_layout(ncol = 2, axes = "collect", axis_titles = "collect") +
 plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
 theme(plot.tag = element_text(size = 10))

ggsave("figures/Figure_01.png", width = 7, height = 6)
ggsave("figures/Figure_01.pdf", width = 7, height = 6)
