library(ggplot2)
library(dplyr)

df <- readRDS("output/15_binom_dens_4species.rds")

# first plot is for the tail probability of sigma
d <- dplyr::filter(df, sensitivity == "sigma_cutoff")

d$sensitivity <- as.factor(d$sensitivity)
levels(d$sensitivity) <- c("sigma[0]")
d$prior.sigma_cutoff <- as.factor(d$prior.sigma_cutoff)

p1 <- ggplot(
  d,
  aes(n, dens_ll, col = prior.sigma_cutoff)
) +
  geom_line() +
  #facet_wrap(~sensitivity, labeller = label_parsed) +
  xlab("Knots") +
  ylab("Predictive density") +
  labs(col = expression(paste(sigma[0]))) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) + 
  scale_color_viridis(discrete=TRUE,end=0.8)

# second plot is for the cutoff of sigma probabilites
d <- dplyr::filter(df, sensitivity == "sigma_p")

d$sensitivity <- as.factor(d$sensitivity)
levels(d$sensitivity) <- c("")#c("Pr(sigma)")
d$prior.sigma_thresh <- as.factor(d$prior.sigma_thresh)

p2 <- ggplot(
  d,
  aes(n, dens_ll, col = prior.sigma_thresh)
) +
  geom_line(alpha=0.7) +
  #facet_wrap(~sensitivity, labeller = label_parsed) +
  xlab("Knots") +
  ylab("Predictive density") +
  #labs(col = expression(paste("Pr(", sigma, ") < ", sigma[0]))) +
  labs(col = expression(theta[sigma])) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) + 
  scale_color_viridis(discrete=TRUE,end=0.8)


d <- dplyr::filter(df, sensitivity == "range_cutoff")

d$sensitivity <- as.factor(d$sensitivity)
levels(d$sensitivity) <- c("kappa[0]")
d$prior.range_cutoff <- as.factor(d$prior.range_cutoff)

p3 <- ggplot(
  d,
  aes(n, dens_ll, col = prior.range_cutoff)
) +
  geom_line() +
  #facet_wrap(~sensitivity, labeller = label_parsed) +
  xlab("Knots") +
  ylab("Predictive density") +
  labs(col = expression(paste(kappa[0]))) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) + 
  scale_color_viridis(discrete=TRUE,end=0.8)
