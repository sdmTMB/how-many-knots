library(ggplot2)
df_gamma = readRDS("output/03_gamma_dens_4species.rds")
df_gamma$species = as.character(df_gamma$species)
df_gamma$species[which(df_gamma$species=="dover sole")] = "Dover Sole"
df_gamma$species[which(df_gamma$species=="petrale sole")] = "Petrale Sole"
df_gamma$species[which(df_gamma$species=="darkblotched rockfish")] = "Darkblotched Rockfish"
df_gamma$species[which(df_gamma$species=="lingcod")] = "Lingcod"

g1 = ggplot(df_gamma, aes(n, dens_ll)) + 
  geom_point(size = 3, col="darkblue",alpha=0.5) +
  #geom_smooth() + 
  theme_bw() + 
  facet_wrap(~species, scale="free_y") + 
  xlab("Knots")+
  ylab("Gamma predictive density") + 
  theme(strip.background =element_rect(fill="white"))

pdf("figures/Figure_S1.pdf")
g1
dev.off()

jpeg("figures/Figure_S1.jpeg")
g1
dev.off()



