library(ggplot2)
df_bin <- readRDS("output/04_binomial_cv_df.rds")
folds_blocks <- readRDS("output/folds_blocks.rds")
df_bin <- dplyr::left_join(df_bin, folds_blocks) %>%
  dplyr::rename(blocks = b)

# how do range and folds affect

# first use the 10-fold mdoels to look at block effect
df_bin$blocks <- as.factor(df_bin$blocks)
g1 <- ggplot(df_bin, aes(n, dens_ll, group = blocks, col = blocks)) +
  geom_point(size = 2, alpha = 0.5) +
  # geom_line() +
  geom_smooth(se = FALSE, method = "loess") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  xlab("Knots") +
  labs(colour = "Blocks") +
  ylab("Binomial predictive density")
g1

pdf("figures/Figure_03.pdf")
g1
dev.off()

jpeg("figures/Figure_03.jpeg")
g1
dev.off()
