library(sf)
library(inlabru)
library(INLA)
library(dplyr)
library(fmesher)
# https://www.maths.ed.ac.uk/~flindgre/2018/07/22/spatially-varying-mesh-quality/

dover <- readRDS("data/doversole_cleaned.rds")
# convert coordinates to km
dover$X <- dover$X / 1000
dover$Y <- dover$Y / 1000
# create occurrence field
dover$present <- ifelse(dover$cpue_kg_km2 > 0, 1, 0)
# filter out 2018
dover <- dover[which(dover$year == 2018), ]

set.seed(2021)
n_folds <- 10
n_blocks <- 10
# first assign blocks
dover$block <- 1
for (jj in 1:n_blocks) {
  dover$block[which(dover$latitude < quantile(dover$latitude, 1 - jj * (1 / n_blocks)))] <- jj + 1
}
# now assign folds
block_fold <- data.frame(
  "block" = 1:n_blocks,
  "fold" = rep(1:n_folds, n_blocks / n_folds)
)
dover <- dplyr::left_join(dover, block_fold)


coordinates(dover) <- c("X", "Y")

# initial loop over the cutoff values
df <- data.frame(
  "cutoff" = seq(3, 120, by = 6),
  "n" = NA,
  "dens_ll" = NA
)
# create boundary, same for all meshes
boundary <- inla.nonconvex.hull(coordinates(dover),
  convex = -0.05
)

mesh1 <- inla.mesh.2d(
  loc = coordinates(dover),
  boundary = boundary,
  offset = c(-0.05, -0.05),
  max.n = 5000,
  max.n.strict = 5000,
  cutoff = 20,
  max.edge = c(100, 500)
)

mesh2 <- inla.mesh.2d(
  loc = coordinates(dover),
  boundary = boundary,
  offset = c(-0.05, -0.05),
  max.n = 5000,
  max.n.strict = 5000,
  cutoff = 80,
  max.edge = c(100, 500)
)

plot_df <- data.frame(coordinates(dover))

p2 <- ggplot() +
  gg(mesh2) +
  theme_bw() +
  xlab("Eastings") +
  ylab("Northings") +
  geom_point(data = plot_df, aes(X, Y), alpha = 0.3, col = "purple", size = 0.4)

p1 <- ggplot() +
  gg(mesh1) +
  theme_bw() +
  xlab("Eastings") +
  ylab("Northings") +
  geom_point(data = plot_df, aes(X, Y), alpha = 0.3, col = "purple", size = 0.4)

pdf("figures/Figure_01.pdf")
gridExtra::grid.arrange(p2, p1, nrow = 1)
dev.off()

jpeg("figures/Figure_01.jpeg")
gridExtra::grid.arrange(p2, p1, nrow = 1)
dev.off()
