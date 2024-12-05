library(ggplot2)
df_bin <- readRDS("output/04_binomial_cv_df.rds")

# how do range and folds affect

# first use the 10-fold mdoels to look at block effect
df_bin$blocks <- as.factor(df_bin$blocks)
g1 <- ggplot(dplyr::filter(df_bin, folds == 10), aes(cutoff, dens_ll, group = blocks, col = blocks)) +
  geom_point(size = 2, alpha = 0.5) +
  # geom_line() +
  geom_smooth(se = FALSE, method = "loess", alpha = 0.5) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  xlab("Cutoff distance (km)") +
  labs(colour = "Blocks") +
  ylab("Log density (test)")

g2 <- ggplot(dplyr::filter(df_bin, folds == 40), aes(n, dens_ll, group = blocks, col = blocks)) +
  geom_point(size = 2, alpha = 0.5) +
  # geom_line() +
  geom_smooth(se = FALSE, method = "loess", alpha = 0.5) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  xlab("Cutoff distance (km)") +
  labs(colour = "Blocks") +
  ylab("Log density (test)")

# pdf("figures/Figure_03.pdf")
# g1
# dev.off()
# 
# jpeg("figures/Figure_03.jpeg")
# g1
# dev.off()
# 
# 
# pdf("figures/Figure_S3.pdf")
# g2
# dev.off()
# 
# jpeg("figures/Figure_S3.jpeg")
# g2
# dev.off()
