# devtools::install_github("seananderson/ggsidekick")
library(ggsidekick)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(raster)
library(scales)
library(patchwork)
library(cowplot)

df <- readRDS(file = "output/temp_model_df.rds")

haul <- readRDS("data/haul_cleaned.rds")
haul$date_formatted <- as.Date(as.numeric(haul$date_formatted), origin = "1970-01-01")
haul$year <- lubridate::year(haul$date_formatted)
haul <- dplyr::filter(
  haul, year == 2018,
  !is.na(temperature_at_gear_c_der)
)

# plot predictive in-sample density
df_train <- dplyr::select(df, -dens_ll_test) |>
  dplyr::rename(dens_ll = dens_ll_train) |>
  dplyr::mutate(name = "Train")
df_test <- dplyr::select(df, -dens_ll_train) |>
  dplyr::rename(dens_ll = dens_ll_test) |>
  dplyr::mutate(name = "Test")

# normalize values based on sample sizes -- in 01_bottom_temp script
df_train$dens_ll <- df_train$dens_ll / 590
df_test$dens_ll <- df_test$dens_ll / 65.5

df_plot <- rbind(df_train, df_test) |>
  dplyr::filter(converged == TRUE, n < 655, cutoff <= 300)

# Individual plots with explicit white backgrounds
p2a <- rbind(df_train, df_test) |>
  dplyr::filter(converged == TRUE, n < 655, cutoff <= 300) |>
  ggplot(aes(cutoff, dens_ll, group = name, col = name)) +
  geom_line() +
  xlab("Cutoff distance (km)") +
  ylab("Log density") +
  scale_color_viridis_d(option = "magma", begin = 0.2, end = 0.8, name = "Data") +
  ggsidekick::theme_sleek()

p2b <- rbind(df_train, df_test) |>
  dplyr::filter(converged == TRUE, n < 655, cutoff <= 300) |>
  ggplot(aes(n, dens_ll, group = name, col = name)) +
  geom_line() +
  xlab("Mesh vertices (n)") +
  ylab("Log density") +
  scale_color_viridis_d(option = "magma", begin = 0.2, end = 0.8, name = "Data") +
  ggsidekick::theme_sleek()

# get shoreline data
map_data <- rnaturalearth::ne_countries(
  scale = "large",
  returnclass = "sf", country = c("canada", "united states of america", "mexico"))
# Crop the polygon for plotting and efficiency:
# st_bbox(map_data) # find the rough coordinates
coast <- suppressWarnings(suppressMessages(
  sf::st_crop(map_data,
              c(xmin = -130, ymin = 30, xmax = -107, ymax = 55))))
coast_proj <- sf::st_transform(coast, crs = 3157) # zone 10

sf::st_bbox(coast_proj)
# Map plot with explicit white background
p1 <- ggplot(coast_proj) +
  geom_sf(data = coast_proj, fill = "grey80") +
  geom_point(data = haul, aes(x = X*1000, y = Y*1000, col = temperature_at_gear_c_der), size = 0.3) +
  # scale_color_gradient2(name = "\u00B0C", midpoint = mean(haul$temperature_at_gear_c_der)) +
  scale_color_viridis_c(name = "\u00B0C", option = "G") +
  labs(x = "Eastings", y = "Northings") +
  ggsidekick::theme_sleek() +
  coord_sf(xlim = c(235000, 1015000), ylim = c(3586000, 5560000))

p1

combined_plot <- (p2a / p2b) + plot_layout(guides = "collect")

design <- "
  12
"

p1 + combined_plot + 
  plot_layout(design = design) + 
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
 theme(plot.tag = element_text(size = 10))

ggsave(filename = "figures/Figure_02.png", width = 7, height = 6, bg = "white")
ggsave(filename = "figures/Figure_02.pdf", width = 7, height = 6)

# Make 2nd plot of parameters
df <- readRDS(file = "output/temp_model_df.rds")
lu <- dplyr::select(df, cutoff, n) |> distinct() # to join on n (vertices)

all_df <- readRDS(file = "output/temp_model_all_est.rds")
all_df <- dplyr::left_join(all_df, df[,c("cutoff","n")]) # join in mesh n
all_df$term[which(all_df$term == "(Intercept)")] <- "Intercept"
all_df$term[which(all_df$term == "I(zday^2)")] <- "day2"
all_df$term[which(all_df$term == "phi")] <- "Obs SD"
all_df$term[which(all_df$term == "range")] <- "Spatial range"
all_df$term[which(all_df$term == "sigma_O")] <- "Spatial SD"
all_df$term[which(all_df$term == "zday")] <- "day"
all_df$term <- factor(all_df$term, levels = c("Intercept", "day", "day2", "Spatial range", "Spatial SD", "Obs SD"))

# dplyr::filter(all_df, term == "Intercept") |> as.data.frame()

parsed_labels <- c(
  "Intercept" = "Intercept",
  "day" = "day",
  "day2" = "day^2",
  "Spatial range" = "Spatial~range~kappa",  
  "Spatial SD" = "Spatial~sigma",
  "Obs SD" = "Observation~sigma"
)

ggplot(all_df, aes(n, mean_estimate)) +
  geom_ribbon(aes(
    ymin = mean_estimate - 2 * mean_se,
    ymax = mean_estimate + 2 * mean_se
  ), alpha = 0.3) +
  geom_line() +
  facet_wrap(
    ~term,
    scales = "free_y",
    labeller = labeller(term = as_labeller(parsed_labels, label_parsed))
  ) +
  xlab("Mesh vertices") +
  ylab("Estimate") +
  scale_x_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  ggsidekick::theme_sleek()

ggsave(filename = "figures/Figure_02b.png", width = 7, height = 5)
ggsave(filename = "figures/Figure_02b.pdf", width = 7, height = 5)
