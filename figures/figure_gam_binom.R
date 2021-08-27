library(ggplot2)
df_bin <- readRDS("output/10_binom_dens_4species.rds")

df_bin$species <- as.character(df_bin$species)
df_bin$species[which(df_bin$species == "dover sole")] <- "Dover sole"
df_bin$species[which(df_bin$species == "petrale sole")] <- "Petrale sole"
df_bin$species[which(df_bin$species == "darkblotched rockfish")] <- "Darkblotched rockfish"
df_bin$species[which(df_bin$species == "lingcod")] <- "Lingcod"

g1 <- ggplot(dplyr::filter(df_bin, n <= 500), aes(n, dens_ll)) +
  geom_point(size = 2, col = "black", alpha = 0.8) +
  geom_smooth(se = FALSE, colour = "black") +
  theme_bw() +
  facet_wrap(~species, scale = "free_y") +
  xlab("Knots") +
  ylab("Binomial predictive density") +
  theme(strip.background = element_rect(fill = "white"))
g1

ggsave("figures/Figure_binom_TMB.pdf", width = 8, height = 6)
ggsave("figures/Figure_binom_TMB.jpeg", width = 8, height = 6)
