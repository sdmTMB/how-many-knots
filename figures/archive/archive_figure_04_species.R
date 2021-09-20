library(ggplot2)
df_bin <- readRDS("output/05_binomial_species_df.rds")
dover <- readRDS("output/02_binomial_df.rds")
dover$blocks <- 10
dover$species <- "dover sole"

df_bin <- rbind(dover, df_bin)

df_bin <- dplyr::filter(df_bin, species != "sablefish")
g1 <- ggplot(df_bin, aes(n, dens_ll)) +
  geom_point(size = 3, col = "darkblue", alpha = 0.5) +
  # geom_smooth() +
  theme_bw() +
  facet_wrap(~species, scale = "free_y") +
  xlab("Knots") +
  ylab("Binomial predictive density") +
  theme(strip.background = element_rect(fill = "white"))

pdf("figures/Figure_02.pdf")
g1
dev.off()

jpeg("figures/Figure_02.jpeg")
g1
dev.off()
