library(ggplot2)
df_gamma <- readRDS("output/03_gamma_dens_4species.rds")
df_gamma$species <- as.character(df_gamma$species)
df_gamma$species[which(df_gamma$species == "dover sole")] <- "Dover Sole"
df_gamma$species[which(df_gamma$species == "petrale sole")] <- "Petrale Sole"
df_gamma$species[which(df_gamma$species == "darkblotched rockfish")] <- "Darkblotched Rockfish"
df_gamma$species[which(df_gamma$species == "lingcod")] <- "Lingcod"

g1 <- ggplot(df_gamma, aes(cutoff, dens_ll_train)) +
  geom_point(size = 3, col = "darkblue", alpha = 0.5) +
  # geom_smooth() +
  theme_bw() +
  facet_wrap(~species, scale = "free_y", nrow = 1) +
  xlab("") +
  ylab("In sample log density") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text.x = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

g2 <- ggplot(df_gamma, aes(cutoff, dens_ll)) +
  geom_point(size = 3, col = "darkblue", alpha = 0.5) +
  # geom_smooth() +
  theme_bw() +
  facet_wrap(~species, scale = "free_y", nrow = 1) +
  xlab("Knots") +
  ylab("Out of sample log density") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text.x = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
# 
# pdf("figures/Figure_S1.pdf")
# gridExtra::grid.arrange(g1, g2)
# dev.off()
# 
# jpeg("figures/Figure_S1.jpeg")
# gridExtra::grid.arrange(g1, g2)
# dev.off()
