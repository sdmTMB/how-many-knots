library(ggplot2)
library(dplyr)
library(viridis)
library(dplyr)
library(tidyr)
library(purrr)

theme_set(ggsidekick::theme_sleek() +
    theme(legend.position = "bottom"))

ggsave2 <- function(filename, ...) {
  ggsave(paste0(filename, ".png"), ...)
  ggsave(paste0(filename, ".pdf"), ...)
}

index <- readRDS("output/04_index_estimates.rds")
# filter out species - cutoff combinations that don't pass sanity()
df <- readRDS("output_index_df.rds")
df <- dplyr::group_by(df, species) |>
  dplyr::mutate(n_conv = length(which(converged==TRUE))) 

index <- left_join(index, df[,c("cutoff","species","converged", "n_conv")])
index <- dplyr::filter(index, n_conv >= 3, converged == TRUE,
                       species != "stripetail rockfish") # filter out species that have 0 - 1 models

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

p1 <- ggplot(indx_se, aes(species, mean_se, col = cutoff)) +
  xlab("") +
  ylab("Mean standard error of indices") +
  geom_point(alpha = 0.7, position = position_dodge(0.4)) +
  coord_flip() +
  theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.2)) +
  scale_color_viridis_d(option = "magma", begin = 0.2, end = 0.8, name = "Cutoff (km)") + theme(legend.position = "right") +
  scale_y_continuous(lim = c(0, NA), expand = expansion(mult = c(0, 0.05)))
p1
ggsave2("figures/mean_index_se", height = 5, width = 5.25)


index_wide <- index |>
  select(species, year, cutoff, est) |>
  pivot_wider(names_from = cutoff, values_from = est, names_prefix = "cutoff_")

cor_results <- index_wide |>
  group_by(species) |>
  summarise(
    cor50 = cor(cutoff_15, cutoff_50, use = "complete.obs"),
    cor100 = cor(cutoff_15, cutoff_100, use = "complete.obs"),
    ratio50 = mean(cutoff_15 / cutoff_50),
    ratio100 = mean(cutoff_15 / cutoff_100),
    .groups = "drop"
  )

cor_results_long <- cor_results |>
  pivot_longer(
    cols = starts_with("cor"),
    names_to = "comparison",
    values_to = "correlation"
  ) |>
  dplyr::rename(cutoff = comparison)
cor_results_long$cutoff[which(cor_results_long$cutoff == "cor50")] <- 50
cor_results_long$cutoff[which(cor_results_long$cutoff == "cor100")] <- 100
cor_results_long$cutoff <- factor(cor_results_long$cutoff, levels = c(15, 50, 100))


my_cols <- viridis(n = 3, option = "magma", begin = 0.2, end = 0.8)

cor_results_long$species <- factor(cor_results_long$species, levels = rev(unique(cor_results_long$species)))
p2 <- ggplot(cor_results_long, aes(species, correlation, col = cutoff)) +
  xlab("") +
  ylab("Correlation with 15km cutoff") +
  geom_point(alpha = 0.7, position = position_dodge(0.4)) +
  coord_flip() +
  theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.2)) +
  scale_color_manual(
    values = my_cols[2:3],
    name = "Cutoff",
    labels = c("50", "100")
  ) + 
  theme(legend.position = "right") +
  scale_y_continuous(lim = c(0.5, NA), expand = expansion(mult = c(0, 0.05)))

ratio_results_long <- cor_results |>
  pivot_longer(
    cols = starts_with("ratio"),
    names_to = "comparison",
    values_to = "ratio"
  ) |>
  dplyr::rename(cutoff = comparison)
ratio_results_long$cutoff[which(ratio_results_long$cutoff == "ratio50")] <- 50
ratio_results_long$cutoff[which(ratio_results_long$cutoff == "ratio100")] <- 100
ratio_results_long$cutoff <- factor(ratio_results_long$cutoff, levels = c(15, 50, 100))
ratio_results_long$species <- factor(ratio_results_long$species, levels = rev(unique(ratio_results_long$species)))

p3 <- ggplot(ratio_results_long, aes(species, ratio, col = cutoff)) +
  xlab("") +
  ylab("Ratio with 15km cutoff") +
  geom_point(alpha = 0.7, position = position_dodge(0.4)) +
  coord_flip() +
  theme(panel.grid.major = element_line(colour = "grey90", linewidth = 0.2)) +
  scale_color_manual(
    values = my_cols[2:3],
    name = "Cutoff",
    labels = c("50", "100")
  ) + 
  theme(legend.position = "right") +
  scale_y_continuous(lim = c(0, NA), expand = expansion(mult = c(0, 0.05)))


gridExtra::grid.arrange(p1, p2, p3, nrow=1)
ggsave2("figures/combo_gfish", height = 7, width = 7)