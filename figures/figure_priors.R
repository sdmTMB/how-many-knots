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
  facet_wrap(~sensitivity, labeller = label_parsed) +
  xlab("Knots") +
  ylab("Predictive density") +
  labs(col = expression(paste(sigma))) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"))


# second plot is for the cutoff of sigma probabilites
d <- dplyr::filter(df, sensitivity == "sigma_p")

d$sensitivity <- as.factor(d$sensitivity)
levels(d$sensitivity) <- c("Pr(sigma)")
d$prior.sigma_thresh <- as.factor(d$prior.sigma_thresh)

p2 <- ggplot(
  d,
  aes(n, dens_ll, col = prior.sigma_thresh)
) +
  geom_line() +
  facet_wrap(~sensitivity, labeller = label_parsed) +
  xlab("Knots") +
  ylab("Predictive density") +
  labs(col = expression(paste("Pr(", sigma, ") < ", sigma[0]))) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"))
