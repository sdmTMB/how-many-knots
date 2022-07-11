library(ggplot2)
library(ggpubr)

# random cross validation
df_gamma <- readRDS("output/03_gamma_dens_4species_random.rds")
df_gamma$blocks = "Random"
# block cross validation
df_gamma_block <- readRDS("output/03_gamma_dens_4species.rds")
df_gamma = rbind(df_gamma, df_gamma_block)

df_gamma = dplyr::filter(df_gamma, n <= 500, dens_ll > -10000)
df_gamma$species <- as.character(df_gamma$species)
df_gamma$species[which(df_gamma$species == "dover sole")] <- "Dover Sole"
df_gamma$species[which(df_gamma$species == "petrale sole")] <- "Petrale Sole"
df_gamma$species[which(df_gamma$species == "darkblotched rockfish")] <- "Darkblotched Rockfish"
df_gamma$species[which(df_gamma$species == "lingcod")] <- "Lingcod"
df_gamma$species[which(df_gamma$species == "sablefish")] <- "Sablefish"
df_gamma$species[which(df_gamma$species == "Pacific ocean perch")] <- "POP"

df_gamma$Blocks = as.factor(df_gamma$blocks)

g1 <- ggplot(df_gamma, aes(cutoff, dens_ll_train, group=Blocks, col=Blocks)) +
  geom_point(size = 2,alpha = 0.5) +
  # geom_smooth() +
  theme_bw() +
  facet_wrap(~species, nrow = 1,scale="free") +
  xlab("") +
  ylab("Log density (train)") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text.x = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  scale_color_viridis_d(end=0.7) + 
  theme(axis.text.y = element_text(angle = 90)) + 
  geom_smooth(se=FALSE,size=0.3)

g2 <- ggplot(df_gamma, aes(cutoff, dens_ll,group=Blocks, col=Blocks)) +
  geom_point(size = 2,alpha = 0.5) +
  # geom_smooth() +
  theme_bw() +
  facet_wrap(~species, nrow = 1,scale="free") +
  xlab("Cutoff distance (km)") +
  ylab("Log density (test)") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text.x = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  scale_color_viridis_d(end=0.7) + 
  theme(axis.text.y = element_text(angle = 90)) + 
  geom_smooth(se=FALSE,size=0.3)
g = ggarrange(g1, g2, ncol=1,common.legend=TRUE)

pdf("figures/Figure_04.pdf")
g
dev.off()

jpeg("figures/Figure_04.jpeg")
g
dev.off()
