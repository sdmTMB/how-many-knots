library(ggplot2)
library(dplyr)

df_bin <- readRDS("output/09_binom_params_4species.rds")

# join in predictive density
df_bin_dens <- readRDS("output/02_binom_dens_4species.rds")
df_bin <- left_join(df_bin, df_bin_dens[, c("cutoff", "n", "species", "dens_ll")])

df_bin$species <- as.character(df_bin$species)
df_bin$species[which(df_bin$species == "Dover sole")] <- "Dover Sole"
df_bin$species[which(df_bin$species == "petrale sole")] <- "Petrale Sole"
df_bin$species[which(df_bin$species == "darkblotched rockfish")] <- "Darkblotched Rockfish"
df_bin$species[which(df_bin$species == "lingcod")] <- "Lingcod"

g1 <- ggplot(df_bin, aes(n, range)) +
  # geom_pointrange(aes(ymin=range_lo,ymax=range_hi),col="darkblue",alpha=0.5) +
  geom_point() +
  facet_wrap(~species, scale = "free", nrow = 1) +
  theme_bw() +
  xlab("Knots") +
  ylab(expression(paste("Matérn range ", kappa, " (km)"))) +
  theme(strip.background = element_rect(fill = "white"))

g2 <- ggplot(df_bin, aes(n, sd)) +
  # geom_pointrange(aes(ymin=sd_lo,ymax=sd_hi),col="darkblue",alpha=0.5) +
  geom_point() +
  facet_wrap(~species, scale = "free", nrow = 1) +
  theme_bw() +
  xlab("Knots") +
  ylab(expression(paste("Matérn ", sigma[omega]))) +
  theme(strip.background = element_rect(fill = "white"))


g3 <- ggplot(df_bin, aes(n, range_sd)) +
  # geom_pointrange(aes(ymin=range_lo,ymax=range_hi),col="darkblue",alpha=0.5) +
  geom_point() +
  facet_wrap(~species, scale = "free", nrow = 1) +
  theme_bw() +
  xlab("Knots") +
  ylab(expression(paste("Standard deviation of Matérn range"))) +
  theme(strip.background = element_rect(fill = "white"))

g4 <- ggplot(df_bin, aes(n, sd_sd)) +
  # geom_pointrange(aes(ymin=sd_lo,ymax=sd_hi),col="darkblue",alpha=0.5) +
  geom_point() +
  facet_wrap(~species, scale = "free", nrow = 1) +
  theme_bw() +
  xlab("Knots") +
  ylab(expression(paste("Standard deviation of Matérn ", sigma[omega]))) +
  theme(strip.background = element_rect(fill = "white"))


pdf("figures/Figure_05.pdf")
gridExtra::grid.arrange(g1, g2, g3, g4, ncol = 1)
dev.off()

jpeg("figures/Figure_05.jpeg")
gridExtra::grid.arrange(g1, g2, g3, g4, ncol = 1)
dev.off()
#
# for(i in 1:length(unique(df_bin$species))) {
#   indx = which(df_bin$species==unique(df_bin$species)[i])
#   #print(cor(df_bin[indx,c("n", "b","sd","range","dens_ll")]))
# }
