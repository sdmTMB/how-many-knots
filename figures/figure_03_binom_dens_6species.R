library(ggplot2)
library(ggpubr)

df_random <- readRDS("output/02_binom_dens_6species_random.rds")
df_random$Folds = "Random"
df_band <- readRDS("output/02_binom_dens_6species_block.rds")
df_band$Folds = paste(df_band$blocks, "bands")
df_block <- readRDS("output/02_binom_dens_6species_blockCV.rds")
df_block$Folds = paste("Range", df_block$range)

df <- rbind(df_band[,c("cutoff","species","n","dens_ll","Folds")],
            df_block[,c("cutoff","species","n","dens_ll","Folds")],
            df_random[,c("cutoff","species","n","dens_ll","Folds")])

df$species <- as.character(df$species)
df$species[which(df$species == "dover sole")] <- "Dover Sole"
df$species[which(df$species == "petrale sole")] <- "Petrale Sole"
df$species[which(df$species == "darkblotched rockfish")] <- "Darkblotched"
df$species[which(df$species == "lingcod")] <- "Lingcod"
df$species[which(df$species == "sablefish")] <- "Sablefish"
df$species[which(df$species == "Pacific ocean perch")] <- "POP"

df$Folds = factor(df$Folds, levels = c("Random", "10 bands", "20 bands",
                                       "Range 25", "Range 50", "Range 75",
                                       "Range 100", "Range 125", "Range 150"))

# make all values relative to finest mesh
df = dplyr::group_by(df, species, Folds) %>%
  dplyr::mutate(dens_ll = dens_ll - mean(dens_ll,na.rm=T))

# g1 <- ggplot(df_bin, aes(cutoff, dens_ll_train, group=Blocks, col=Blocks)) +
#   geom_point(size = 2,alpha = 0.5) +
#   # geom_smooth() +
#   theme_bw() +
#   facet_wrap(~species, nrow = 1,scale="free") +
#   xlab("") +
#   ylab("Log density (train)") +
#   theme(strip.background = element_rect(fill = "white")) +
#   theme(strip.text.x = element_text(size = 7)) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
#   scale_color_viridis_d(end=0.7) + 
#   theme(axis.text.y = element_text(angle = 90)) + 
#   geom_smooth(se=FALSE,size=0.3)

# custom colors
library(RColorBrewer)
rand_col <- grey(0.3)
band_col <- magma(2, begin=0.4, end=0.8)
block_col <- viridis(6,end=0.8)
cols <- c(rand_col, band_col, block_col)
names(cols) = levels(df$Folds)
custom_cols <- scale_color_manual(name = "Levels", values=cols)

g <- ggplot(dplyr::filter(df, Folds %in% c("Range 50","Range 100","Range 150")== FALSE), aes(cutoff, dens_ll, group=Folds, col=Folds)) +
  geom_point(size = 2, alpha = 0.5) +
  # geom_smooth() +
  theme_bw() +
  facet_wrap(~species, nrow = 2, scale="free") +
  xlab("Cutoff distance (km)") +
  ylab("Log density (test)") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text.x = element_text(size = 7)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  #scale_color_viridis_d(end=0.7) + 
  theme(axis.text.y = element_text(angle = 90)) + 
  geom_smooth(se=FALSE,size=0.3, method="loess",span=0.1) + 
  custom_cols
#pdf("test.pdf")
#g = ggarrange(g1, g2, ncol=1,common.legend=TRUE)
#dev.off()
# pdf("figures/Figure_03.pdf")
# g
# dev.off()
# 
# jpeg("figures/Figure_03.jpeg")
# g
# dev.off()
