library(ggplot2)
library(dplyr)
theme_set(ggsidekick::theme_sleek() +
    theme(legend.position = "bottom"))

ggsave2 <- function(filename, ...) {
  ggsave(paste0(filename, ".png"), ...)
  ggsave(paste0(filename, ".pdf"), ...)
}

index <- readRDS("output/04_index_estimates.rds")

index$cutoff <- as.factor(index$cutoff)
# index <- dplyr::filter(index, species == "arrowtooth flounder")
capitalize_first <- function(x) {
  paste0(toupper(substring(x, 1, 1)), substring(x, 2))
}
index$species <- capitalize_first(index$species)

ggplot(index, aes(year, log_est, group = cutoff, col = cutoff)) +
  geom_pointrange(aes(ymin = log_est - 1.96 * se, ymax = log_est + 1.96 * se), position = position_dodge(0.5), alpha = 0.7, fatten = 0.1) +
  xlab("Year") +
  ylab("Ln estimate") +
  facet_wrap(~species, scale = "free_y", ncol = 5) +
  scale_color_viridis_d(option = "magma", begin = 0.2, end = 0.8, name = "Cutoff (km)")
ggsave2("figures/SI_figure_indices", height = 7, width = 9)

# do caterpillar plot of average error by species
indx_se <- dplyr::group_by(index, species, cutoff) |>
  dplyr::summarize(mean_se = mean(se))

indx_se$species <- factor(indx_se$species, levels = rev(unique(indx_se$species)))

ggplot(indx_se, aes(species, mean_se, col = cutoff)) +
  xlab("") +
  ylab("Mean standard error of indices") +
  geom_point(alpha = 0.7, position = position_dodge(0.4)) +
  coord_flip() +
  theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.2)) +
  scale_color_viridis_d(option = "magma", begin = 0.2, end = 0.8, name = "Cutoff (km)") + theme(legend.position = "right") +
  scale_y_continuous(lim = c(0, NA), expand = expansion(mult = c(0, 0.05)))
ggsave2("figures/mean_index_se", height = 5, width = 5.25)
