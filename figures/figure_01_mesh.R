library(sf)
#library(inlabru)
#library(INLA)
library(dplyr)
#library(fmesher)
library(lubridate)
library(cowplot)
# https://www.maths.ed.ac.uk/~flindgre/2018/07/22/spatially-varying-mesh-quality/

haul <- readRDS("data/haul_cleaned.rds")
haul$date_formatted <- as.Date(as.numeric(haul$date_formatted), origin = "1970-01-01")

haul$year <- lubridate::year(haul$date_formatted)
haul <- dplyr::filter(haul, year == 2018)

coordinates(haul) <- c("X", "Y")

# create boundary, same for all meshes
boundary <- inla.nonconvex.hull(coordinates(haul),
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
  theme_bw() +
  xlab("Eastings") +
  ylab("Northings") +
  geom_point(data = plot_df, aes(X, Y), alpha = 0.3, col = viridis(1), size = 0.4) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t = 5, r = 6, b = 20, l = 20))

p1 <- ggplot() +
  gg(mesh1) +
  theme_bw() +
  xlab("Eastings") +
  ylab("Northings") +
  geom_point(data = plot_df, aes(X, Y), alpha = 0.3, col = viridis(1), size = 0.4) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t = 5, r = 6, b = 20, l = 20))

combined_plot <- cowplot::plot_grid(
  p2, p1, 
  align = "hv",
  axis = "tblr",
  nrow=1
)

final_plot <- ggdraw() + 
  draw_plot(combined_plot, 0, 0, 1, 1) + 
  draw_label("Eastings", x = 0.5, y = -0.05, vjust = -3, size=12) + 
  draw_label("Northings", x = 0, y = 0, angle = 90, vjust = 1.5, hjust=-3.7, size=12)
  
ggsave(final_plot, filename = "figures/Figure_01.png", width=7, height = 6)

