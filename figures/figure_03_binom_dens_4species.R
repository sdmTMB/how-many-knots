library(ggplot2)
library(ggpubr)
df_bin <- readRDS("output/02_binom_dens_4species.rds")

df_bin$species <- as.character(df_bin$species)
df_bin$species[which(df_bin$species == "dover sole")] <- "Dover Sole"
df_bin$species[which(df_bin$species == "petrale sole")] <- "Petrale Sole"
df_bin$species[which(df_bin$species == "darkblotched rockfish")] <- "Darkblotched"
df_bin$species[which(df_bin$species == "lingcod")] <- "Lingcod"
df_bin$species[which(df_bin$species == "sablefish")] <- "Sablefish"
df_bin$species[which(df_bin$species == "Pacific ocean perch")] <- "POP"

df_bin$Blocks = as.factor(df_bin$blocks)

# make all values relative to finest mesh
df_bin = dplyr::group_by(df_bin, species, Blocks) %>%
  dplyr::mutate(dens_ll = dens_ll - dens_ll[which(cutoff==3)],
                dens_ll_train = dens_ll_train - dens_ll_train[which(cutoff==3)])

g1 <- ggplot(df_bin, aes(cutoff, dens_ll_train, group=Blocks, col=Blocks)) +
  geom_point(size = 2,alpha = 0.5) +
  # geom_smooth() +
  theme_bw() +
  facet_wrap(~species, nrow = 1) +
  xlab("") +
  ylab("In sample log density") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text.x = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  scale_color_viridis_d(end=0.7)

g2 <- ggplot(df_bin, aes(cutoff, dens_ll, group=Blocks, col=Blocks)) +
  geom_point(size = 2, alpha = 0.5) +
  # geom_smooth() +
  theme_bw() +
  facet_wrap(~species, nrow = 1) +
  xlab("Cutoff distance (km)") +
  ylab("Out of sample log density") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text.x = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  scale_color_viridis_d(end=0.7)

#pdf("test.pdf")
g = ggarrange(g1, g2, ncol=1,common.legend=TRUE)
#dev.off()
#pdf("figures/Figure_02.pdf")
#g = gridExtra::grid.arrange(g1, g2)
#dev.off()
# 
# jpeg("figures/Figure_02.jpeg")
# gridExtra::grid.arrange(g1, g2)
# dev.off()
