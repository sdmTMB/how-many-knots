library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(raster)
library(scales)
library(patchwork)
library(cowplot)

df = readRDS(file = "output/temp_model_df.rds")

haul <- readRDS("data/haul_cleaned.rds")
haul$date_formatted <- as.Date(as.numeric(haul$date_formatted), origin = "1970-01-01")
haul$year <- lubridate::year(haul$date_formatted)
haul <- dplyr::filter(haul, year==2018,
                          !is.na(temperature_at_gear_c_der))

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
  geom_line(linewidth = 1.5, alpha = 0.8) +
  xlab("Cutoff distance (km)") + 
  ylab("Log density") + 
  scale_color_viridis_d(option = "magma", begin = 0.2, end = 0.8, name = "Data") +
  theme_bw() + 
  theme(
    strip.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title.y = element_blank()
  )

p2b <- rbind(df_train, df_test) |> 
  dplyr::filter(converged == TRUE, n < 655, cutoff <= 300) |> 
  ggplot(aes(n, dens_ll, group = name, col = name)) +
  geom_line(linewidth = 1.5, alpha = 0.8) +
  xlab("Mesh vertices (n)") + 
  ylab("Log density") + 
  scale_color_viridis_d(option = "magma", begin = 0.2, end = 0.8, name = "Data") +
  theme_bw() + 
  theme(
    strip.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title.y = element_blank()
  )

# Combined plot with patchwork
combined_plot <- (p2a / p2b) +
  plot_layout(guides = "collect") &
  theme(
    plot.margin = margin(3, 3, 3, 5),
    plot.background = element_rect(fill = "white", color = NA) # Ensure white background
  )

final_plot <- ggdraw() +
  draw_label("Average log-likelihood", x = 0.03, y = 0.5, angle = 90, vjust = 0.3, size = 12) +
  draw_plot(combined_plot, x = 0.05, y = 0, width = 0.9, height = 1) +
  theme(plot.background = element_rect(fill = "white", color = NA)) # Ensure white background

# Map plot with explicit white background
p1 <- ggplot(coast_proj) +
  geom_sf(fill = "beige") +
  geom_point(data = haul, aes(x = X * 1000, y = Y * 1000, col = temperature_at_gear_c_der), size = 0.2) +
  scale_color_gradient2(name = "\u00B0C", midpoint = mean(haul$temperature_at_gear_c_der)) +
  xlim(295039.2, 1005243.2) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "grey30", color = NA), # Ensure white panel
    panel.grid = element_blank()
  )

# Combine map and final plot with cowplot
g <- cowplot::plot_grid(
  p1, final_plot,
  nrow = 1
)

# Save plot with white background
ggsave(plot = g, filename = "figures/Figure_02.png", width = 7, height = 6, bg = "white")


# Make 2nd plot of parameters 
all_df <- readRDS(file = "output/temp_model_all_est.rds")
all_df$term[which(all_df$term=="(Intercept)")] <- "Intercept"
all_df$term[which(all_df$term=="I(zday^2)")] <- "day2"
all_df$term[which(all_df$term=="phi")] <- "Obs SD"
all_df$term[which(all_df$term=="range")] <- "Spatial range"
all_df$term[which(all_df$term=="sigma_O")] <- "Spatial SD"
all_df$term[which(all_df$term=="zday")] <- "day"
all_df$term <- factor(all_df$term, levels = c("Intercept","day","day2","Spatial range","Spatial SD","Obs SD"))

ggplot(all_df, aes(cutoff, mean_estimate)) + 
  geom_ribbon(aes(ymin = mean_estimate - 2*mean_se, ymax = mean_estimate + 2*mean_se), fill=viridis(1),alpha=0.3) + 
  geom_line(col=viridis(1)) + 
  facet_wrap(~term, scale="free_y") +
  xlab("Cutoff distance (km)") + 
  ylab("Estimate") + 
  theme_bw() + 
  theme(strip.background = element_rect(fill="white"))
ggsave(filename = "figures/Figure_02b.png", width = 7, height = 5)




