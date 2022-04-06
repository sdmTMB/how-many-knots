library(ggplot2)
df_bin <- readRDS("output/02_binom_dens_4species.rds")

df_bin$species <- as.character(df_bin$species)
df_bin$species[which(df_bin$species == "dover sole")] <- "Dover Sole"
df_bin$species[which(df_bin$species == "petrale sole")] <- "Petrale Sole"
df_bin$species[which(df_bin$species == "darkblotched rockfish")] <- "Darkblotched Rockfish"
df_bin$species[which(df_bin$species == "lingcod")] <- "Lingcod"

g1 <- ggplot(df_bin, aes(cutoff, dens_ll_train)) +
  geom_point(size = 3, col = "darkblue", alpha = 0.5) +
  # geom_smooth() +
  theme_bw() +
  facet_wrap(~species, scale = "free_y", nrow = 1) +
  xlab("") +
  ylab("In sample log density") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text.x = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

g2 <- ggplot(df_bin, aes(cutoff, dens_ll)) +
  geom_point(size = 3, col = "darkblue", alpha = 0.5) +
  # geom_smooth() +
  theme_bw() +
  facet_wrap(~species, scale = "free_y", nrow = 1) +
  xlab("Cutoff distance (km)") +
  ylab("Out of sample log density") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text.x = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#pdf("figures/Figure_02.pdf")
#g = gridExtra::grid.arrange(g1, g2)
#dev.off()
# 
# jpeg("figures/Figure_02.jpeg")
# gridExtra::grid.arrange(g1, g2)
# dev.off()
