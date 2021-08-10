library(ggplot2)
df_bin = readRDS("output/02_binom_dens_4species.rds")

df_bin$species = as.character(df_bin$species)
df_bin$species[which(df_bin$species=="dover sole")] = "Dover sole"
df_bin$species[which(df_bin$species=="petrale sole")] = "Petrale sole"
df_bin$species[which(df_bin$species=="darkblotched rockfish")] = "Darkblotched rockfish"
df_bin$species[which(df_bin$species=="lingcod")] = "Lingcod"

g1 = ggplot(dplyr::filter(df_bin,n<=500), aes(n, dens_ll)) + 
  geom_point(size = 3, col="darkblue",alpha=0.5) +
  #geom_smooth() + 
  theme_bw() + 
  facet_wrap(~species, scale="free_y") + 
  xlab("Knots")+
  ylab("Binomial predictive density") + 
  theme(strip.background =element_rect(fill="white"))

pdf("figures/Figure_02.pdf")
g1
dev.off()

jpeg("figures/Figure_02.jpeg")
g1
dev.off()



