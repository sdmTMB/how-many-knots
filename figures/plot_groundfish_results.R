# devtools::install_github("seananderson/ggsidekick")
library(ggsidekick)
library(ggplot2)
library(viridis)
library(dplyr)
library(scales)
library(patchwork)
theme_set(ggsidekick::theme_sleek() +
    theme(legend.position = "bottom"))

d <- readRDS("output/02_binomial_blockCV.rds")

ggsave2 <- function(filename, ...) {
  ggsave(paste0(filename, ".png"), ...)
  ggsave(paste0(filename, ".pdf"), ...)
}

stripwidth_scale <- scale_color_viridis(option = "magma", begin = 0.2, end = 0.8, name = "Cross validation strip width (km)")

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

d_random <- d_random |>
  mutate(species_wrapped = gsub(" ", "\\\n", species))

d <- d |>
  mutate(species_wrapped = gsub(" ", "\\\n", species))

###### Spatial range plots
# d |>
#   ggplot(aes(cutoff, range, group = bin_width, col = bin_width)) +
#   geom_line() +
#   facet_wrap(~species, scale = "free_y") +
#   xlab("Cutoff distance (km)") +
#   ylab("Estimated spatial range (km)") +
#   stripwidth_scale +
#   geom_point(data = d_random, aes(cutoff, range), col = "black", alpha = 0.5)
# ggsave2("figures/groundfish_range_v_cutoff", height = 6, width = 8)

sub <- dplyr::filter(d, species %in% c("Sablefish", "Widow rockfish", "Lingcod", "Arrowtooth flounder"),
                     bin_width == 50)
sub$species <- as.factor(as.character(sub$species))

p1 <- sub |>
  ggplot(aes(n, range, group = bin_width, col = bin_width)) +
  geom_line() +
  facet_wrap(~species, scale = "free_y", ncol = 1) +
  xlab("Mesh vertices (n)") +
  ylab("Estimated spatial range (km)") +
  stripwidth_scale +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) + 
  theme(legend.position = "none")

p2 <- sub |>
  ggplot(aes(n, present_dens_ll, group = bin_width, col = bin_width)) +
  geom_line() +
  facet_wrap(~species, scale = "free_y", ncol = 1) +
  xlab("Mesh vertices (n)") +
  ylab("Predicted log likelihood") +
  stripwidth_scale +
  theme(legend.position = "none")
# Remove xlab from p1?
p1_clean <- p1 + xlab(NULL)
# Bind the range and LL figs in columns
(p1_clean | p2_clean) +
  plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") &
  xlab("Mesh vertices (n)")
ggsave2("figures/groundfish_range_and_ll_v_n", height = 6, width = 8)

d |>
  ggplot(aes(n, range, group = bin_width, col = bin_width)) +
  geom_line() +
  facet_wrap(~species, scale = "free_y") +
  xlab("Mesh vertices (n)") +
  ylab("Estimated spatial range (km)") +
  stripwidth_scale +
  geom_point(data = d_random, aes(n, range), col = "black", alpha = 0.5) +
	scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  geom_line(data = d_random, aes(n, cutoff), col = "black")
ggsave2("figures/groundfish_range_v_n", height = 6, width = 8)

d |> # divide by 12670 to get average ll per obs
  ggplot(aes(n, present_dens_ll / 12670, group = bin_width, col = bin_width)) +
  geom_line() +
  facet_wrap(~species, scale = "free_y") +
  xlab("Mesh vertices (n)") +
  ylab("Predicted log likelihood") +
  stripwidth_scale +
  geom_point(data = d_random, aes(n, present_dens_ll / 12670), col = "black", alpha = 0.5) # +
# scale_y_continuous(labels = function(x) format(x, scientific=TRUE))
ggsave2("figures/groundfish_loglik_v_n", height = 6, width = 9)

###### Log likelihood plots
d |> # divide by 12670 to get average ll per obs
  ggplot(aes(cutoff, present_dens_ll / 12670, group = bin_width, col = bin_width)) +
  geom_line() +
  facet_wrap(~species, scale = "free_y") +
  xlab("Cutoff distance (km)") +
  ylab("Predicted log likelihood") +
  stripwidth_scale +
  geom_point(data = d_random, aes(cutoff, present_dens_ll / 12670), col = "black", alpha = 0.5) # +
# scale_y_continuous(labels = function(x) format(x, scientific=TRUE))
ggsave2("figures/groundfish_loglik_v_cutoff", height = 6, width = 9)

d |> # divide by 12670 to get average ll per obs
  ggplot(aes(n, present_dens_ll / 12670, group = bin_width, col = bin_width)) +
  geom_line() +
  facet_wrap(~species, scale = "free_y") +
  xlab("Mesh vertices (n)") +
  ylab("Predicted log likelihood") +
  stripwidth_scale +
  geom_point(data = d_random, aes(n, present_dens_ll / 12670), col = "black", alpha = 0.5) # +
# scale_y_continuous(labels = function(x) format(x, scientific=TRUE))
ggsave2("figures/groundfish_loglik_v_vertices", height = 6, width = 9)


ggplot(d, aes(n, sigma_O, group = bin_width, col = bin_width)) +
  geom_point() +
  geom_line() +
  facet_wrap(~species, scale = "free_y") +
  xlab("Mesh vertices (n)") +
  ylab(expression("Estimated spatial " * sigma)) +
  stripwidth_scale +
  geom_point(data = d_random, aes(n, sigma_O), col = "black", alpha = 0.5)
ggsave2("figures/groundfish_sigmaO_v_vertices", height = 6, width = 8.5)

ggplot(d, aes(n, sigma_E, group = bin_width, col = bin_width)) +
  geom_point() +
  geom_line() +
  facet_wrap(~species, scale = "free_y") +
  xlab("Mesh vertices (n)") +
  ylab(expression("Estimated spatiotemporal " * sigma)) +
  stripwidth_scale +
  geom_point(data = d_random, aes(n, sigma_E), col = "black", alpha = 0.5)
ggsave2("figures/groundfish_sigmaE_v_vertices", height = 6, width = 8.75)

ggplot(d, aes(n, sigma_O / sigma_E, group = bin_width, col = bin_width)) +
  geom_point() +
  geom_line() +
  facet_wrap(~species, scale = "free_y") +
  xlab("Mesh vertices (n)") +
  ylab(expression("Ratio of spatial to spatiotemporal " * sigma)) +
  stripwidth_scale +
  geom_point(data = d_random, aes(n, sigma_O / sigma_E), col = "black", alpha = 0.5)
ggsave2("figures/groundfish_sigmaRatio_v_vertices", height = 6, width = 8.7)
