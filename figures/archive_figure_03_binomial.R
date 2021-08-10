library(ggplot2)
df_bin = readRDS("output/04_binomial_cv_df.rds")
df_bin2 = readRDS("output/04_binomial_cv_df2.rds")

df_bin = rbind(df_bin,df_bin2)
df_bin = dplyr::arrange(df_bin, blocks,cutoff)

df_bin$blocks = as.factor(df_bin$blocks)
g1 = ggplot(df_bin, aes(n, dens_ll,group=blocks,col=blocks)) + 
  geom_point(size = 2,alpha=0.5) +
  #geom_smooth() + 
  theme_bw() + 
  xlab("Knots")+
  ylab("Binomial predictive density")

pdf("figures/Figure_03.pdf")
g1
dev.off()




