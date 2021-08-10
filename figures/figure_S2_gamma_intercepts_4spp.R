library(ggplot2)
df = readRDS("output/df_gamma_catchrates_4sp.rds")

df$species = as.character(df$species)
df$species[which(df$species=="dover sole")] = "Dover sole"
df$species[which(df$species=="petrale sole")] = "Petrale sole"
df$species[which(df$species=="darkblotched rockfish")] = "Darkblotched rockfish"
df$species[which(df$species=="lingcod")] = "Lingcod"


g1 = ggplot(dplyr::filter(df,!is.na(sd),sd<0.6), aes(n, mean)) + 
  geom_pointrange(aes(ymin=lo,ymax=hi), col="darkblue",alpha=0.5) +
  #geom_smooth() + 
  theme_bw() + 
  facet_wrap(~species, scale="free") + 
  xlab("Knots")+
  ylab("Intercept") + 
  theme(strip.background =element_rect(fill="white"))

pdf("figures/Figure_S2.pdf")
g1
dev.off()

jpeg("figures/Figure_S2.jpeg")
g1
dev.off()



