library(ggplot2)
library(dplyr)

df_bin = readRDS("output/09_binom_params_4species.rds")

df_bin$species = as.character(df_bin$species)
df_bin$species[which(df_bin$species=="Dover sole")] = "Dover sole"
df_bin$species[which(df_bin$species=="petrale sole")] = "Petrale sole"
df_bin$species[which(df_bin$species=="darkblotched rockfish")] = "Darkblotched rockfish"
df_bin$species[which(df_bin$species=="lingcod")] = "Lingcod"

g1 = ggplot(df_bin, aes(n, range)) + 
  geom_pointrange(aes(ymin=range_lo,ymax=range_hi),col="darkblue",alpha=0.5) + 
  facet_wrap(~species,scale="free") + 
  theme_bw() + 
  xlab("Knots")+
  ylab("Estimated range (km)") + 
  theme(strip.background =element_rect(fill="white"))
  
pdf("figures/Figure_05.pdf")
g1
dev.off()

jpeg("figures/Figure_05.jpeg")
g1
dev.off()



