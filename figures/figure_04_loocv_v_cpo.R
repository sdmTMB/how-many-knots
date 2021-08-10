library(ggplot2)
library(dplyr)
cpo = readRDS("output/06_binomial_cpo_df.rds")
df_bin = readRDS("output/05_binomial_loocv.rds")

# aggregate loocv by cutoff
loocv = dplyr::group_by(df_bin, holdout) %>% 
  dplyr::mutate(tempsum = sum(dens_ll)) %>% 
  dplyr::filter(!is.na(tempsum)) %>% 
  dplyr::group_by(cutoff) %>%
  dplyr::summarize(dens = sum(dens_ll),n=n[1])

cpo = dplyr::filter(cpo, cutoff>=15, cutoff<=75) %>%
  dplyr::select(-dic) %>% 
  dplyr::rename(dens = log_cpo)

# bind together
cpo$Estimate = "CPO"
loocv$Estimate = "LOOCV"
d = rbind(cpo,loocv)

g1 = ggplot(d, aes(n, dens, color = Estimate, group=Estimate)) + 
  geom_point(size = 3,alpha=0.5) +
  #geom_smooth() + 
  theme_bw() + 
  xlab("Knots")+
  ylab("Binomial predictive density")
  
pdf("figures/Figure_04.pdf")
g1
dev.off()

jpeg("figures/Figure_04.jpeg")
g1
dev.off()



