library(ggplot2)
library(viridis)

d <- readRDS("output/02_binomial_blockCV.rds")

capitalize_first <- function(x) {
  paste0(toupper(substring(x, 1, 1)), substring(x, 2))
}
d$species <- capitalize_first(d$species)
# filter out the number of species that don't converge enough
d <- dplyr::filter(d, present_converged == TRUE) |>
  dplyr::group_by(species) |>
  dplyr::mutate(nobs = n()) |>
  dplyr::filter(nobs >= 22) |> # this results in top 20 
  dplyr::select(-nobs)
d$species <- as.factor(as.character(d$species))

# Bring in the random 
d_random <- readRDS("output/02_binomial_randomCV.rds")
d_random$species <- capitalize_first(d_random$species)
d_random <- dplyr::filter(d_random, present_converged == TRUE) |>
  dplyr::filter(species %in% d$species) |>
  dplyr::mutate(bin_width = NA)

# Larger bins result in more widely spaced test regions, beyond the estimated
# range. These regions are no longer correlated with the training data and predictions
# become more uncertain / not as good
dplyr::filter(d) |>
  ggplot(aes(n, present_dens_ll, group = bin_width, col = bin_width)) + 
  geom_line() + 
  facet_wrap(~species, scale = "free_y")



d |>
  ggplot(aes(cutoff, range, group = bin_width, col = bin_width)) + 
  geom_line() + 
  facet_wrap(~species, scale = "free_y") + 
  xlab("Cutoff distance (km)") + ylab("Estimated spatial range (km)") + 
  scale_color_viridis(option="magma", begin = 0.2, end = 0.8, name = "Strip width (km)") + 
  theme_bw() + 
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size=4.5),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) + 
  geom_point(data = d_random, aes(cutoff, range), col="black", alpha=0.5)
ggsave("figures/groundfish_range_v_cutoff.png", height = 5, width = 7)


library(scales)
d |> # divide by 12670 to get average ll per obs
  ggplot(aes(cutoff, present_dens_ll / 12670, group = bin_width, col = bin_width)) + 
  geom_line() + 
  facet_wrap(~species, scale = "free_y") + 
  xlab("Cutoff distance (km)") + ylab("Predicted log likelihood") + 
  scale_color_viridis(option="magma", begin = 0.2, end = 0.8, name = "Strip width (km)") + 
  theme_bw() + 
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size=4.5),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=8)) + 
  geom_point(data = d_random, aes(cutoff, present_dens_ll / 12670), col="black", alpha=0.5)# + 
#scale_y_continuous(labels = function(x) format(x, scientific=TRUE))
ggsave("figures/groundfish_loglik_v_cutoff.png", height = 5, width = 7)


# load in the parameters estimates
ggplot(d, aes(cutoff, range, group=bin_width, col = bin_width)) + 
  geom_point() + 
  facet_wrap(~ species, scale="free_y") + 
  xlab("Cutoff distance (km)") +
  ylab("Estimated spatial range (km)") + 
  scale_color_viridis(option="magma", begin = 0.2, end = 0.8, name = "Strip width (km)") +
  theme_bw() + 
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size=4.5),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=8)) + 
  geom_point(data = d_random, aes(cutoff, range), col="black", alpha=0.5)
ggsave("figures/groundfish_range_v_cutoff.png", height = 5, width = 7)

ggplot(d, aes(cutoff, sigma_O, group=bin_width, col = bin_width)) + 
  geom_point() + 
  facet_wrap(~ species, scale="free_y") + 
  xlab("Cutoff distance (km)") +
  ylab(expression("Estimated spatial " * sigma)) + 
  scale_color_viridis(option="magma", begin = 0.2, end = 0.8, name = "Strip width (km)") +
  theme_bw() + 
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size=4.5),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=8)) + 
  geom_point(data = d_random, aes(cutoff, sigma_O), col="black", alpha=0.5)
ggsave("figures/groundfish_sigmaO_v_cutoff.png", height = 5, width = 7)

ggplot(d, aes(cutoff, sigma_E, group=bin_width, col = bin_width)) + 
  geom_point() + 
  facet_wrap(~ species, scale="free_y") + 
  xlab("Cutoff distance (km)") +
  ylab(expression("Estimated spatiotemporal " * sigma)) + 
  scale_color_viridis(option="magma", begin = 0.2, end = 0.8, name = "Strip width (km)") +
  theme_bw() + 
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(size=4.5),
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=8)) + 
  geom_point(data = d_random, aes(cutoff, sigma_E), col="black", alpha=0.5)
ggsave("figures/groundfish_sigmaE_v_cutoff.png", height = 5, width = 7)
