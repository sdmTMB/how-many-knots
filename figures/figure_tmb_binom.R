library(ggplot2)
library(viridis)
library(RColorBrewer)

df_bin <- readRDS("output/08_binom_dens_4species_TMB.rds")

df_bin$species <- as.character(df_bin$species)
df_bin$species[which(df_bin$species == "dover sole")] <- "Dover sole"
df_bin$species[which(df_bin$species == "petrale sole")] <- "Petrale sole"
df_bin$species[which(df_bin$species == "darkblotched rockfish")] <- "Darkblotched rockfish"
df_bin$species[which(df_bin$species == "lingcod")] <- "Lingcod"
df_bin$species[which(df_bin$species == "Pacific ocean perch")] <- "POP"
df_bin$species[which(df_bin$species == "sablefish")] <- "Sablefish"

df_bin$Folds <- paste0(df_bin$blocks, " bands")

#rand_col <- grey(0.5)
band_col <- magma(2, begin=0.4, end=0.8)
#block_col <- viridis(length(unique(df_block$Folds)),begin=0.2,end=0.7)
#cols <- c(band_col)
names(band_col) = levels(df_bin$Folds)
custom_cols <- scale_color_manual(name = "levels", values=band_col)

df_demean = dplyr::group_by(df_bin, species, Folds) %>%
  dplyr::mutate(dens_ll = dens_ll - mean(dens_ll,na.rm=T))

g1 <- ggplot(dplyr::filter(df_demean, n <= 800), 
             aes(cutoff, dens_ll, group=Folds, col=Folds)) +
  #geom_point(size = 2, alpha = 0.8) +
  geom_smooth(se = FALSE, method="loess",span=0.5) +
  theme_bw() +
  facet_wrap(~species, scale = "free_y") +
  xlab("Cutoff distance") +
  ylab("Log density (test)") +
  theme(strip.background = element_rect(fill = "white")) + 
  custom_cols
g1

#ggsave("figures/Figure_gamma_TMB.pdf", width = 8, height = 6)
#ggsave("figures/Figure_gamma_TMB.jpeg", width = 8, height = 6)
