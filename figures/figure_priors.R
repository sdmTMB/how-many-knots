library(ggplot2)
library(dplyr)
library(viridis)
library(ggpubr)
df <- readRDS("output/15_binom_dens_4species.rds")

# Sean's function in sdmTMB -- 
inla_docs_log <- function(range = 1, sigma = 0.6, matern_range = 5,
                          range_prob = 0.05, matern_SD = 2, SD_prob = 0.05) {
  d <- 2
  dhalf <- d / 2
  lam1 <- -log(range_prob)*matern_range^dhalf
  lam2 <- -log(SD_prob)/matern_SD
  log(dhalf) + log(lam1) + log(range^(-1 - dhalf)) - lam1 * range^-dhalf + 
    log(lam2) - lam2 * sigma
}
ranges <- seq(0.001, 100, length.out = 1001) # adjust as needed
sigmas <- seq(0.001, 10, length.out = 1001) # adjust as needed
id <- data.frame(Var2 = seq_along(ranges), range = ranges)
id1 <- data.frame(Var1 = seq_along(sigmas), sigma = sigmas)

sample_matern <- function(matern_range = 1, matern_SD = 1.5,
                          range_prob = 0.05, SD_prob = 0.05, samples=10000) {
  
  out <- sapply(ranges, function(.range) {
    sapply(sigmas, function(.sigma) {
      inla_docs_log(.range, .sigma, matern_range = matern_range, matern_SD = matern_SD, # adjust as needed
                    range_prob = range_prob, SD_prob = SD_prob)  # adjust as needed
    })
  })
  priors = reshape2::melt(out) %>%
    dplyr::left_join(id, by = "Var2") %>%
    dplyr::left_join(id1, by = "Var1") 
  s = sample(seq(1,nrow(priors)), size=samples, replace=FALSE, prob = exp(priors$value))
  return(priors[s,])
}

#######################################################################
# first plot is for the tail probability of sigma
#######################################################################
d <- dplyr::filter(df, sensitivity == "sigma_cutoff")

d$sensitivity <- as.factor(d$sensitivity)
levels(d$sensitivity) <- c("sigma[0]")
d$prior.sigma_cutoff <- as.factor(d$prior.sigma_cutoff)

p1 <- ggplot(
  d,
  aes(n, dens_ll, col = prior.sigma_cutoff)
) +
  geom_line(alpha=0.4) +
  #facet_wrap(~sensitivity, labeller = label_parsed) +
  xlab("Knots") +
  ylab("") +
  labs(col = expression(paste(sigma[0]))) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) + 
  scale_color_viridis(discrete=TRUE,end=0.8)

# now show the 3 priors 
r1 = sample_matern(matern_range = 20, matern_SD = 1, range_prob = 0.05, SD_prob = 0.05)
r2 = sample_matern(matern_range = 20, matern_SD = 5, range_prob = 0.05, SD_prob = 0.05)
r3 = sample_matern(matern_range = 20, matern_SD = 10, range_prob = 0.05, SD_prob = 0.05)
r1$prior.sigma_cutoff = 1
r2$prior.sigma_cutoff = 5
r3$prior.sigma_cutoff = 10

g1 = ggplot(rbind(r1,r2,r3), aes(sigma, 
                            fill=as.factor(prior.sigma_cutoff),
                            col=as.factor(prior.sigma_cutoff))) + 
  geom_density(alpha=0.4,col=NA) + 
  theme_bw() + 
  theme(legend.position = "none") + 
  xlab(expression(sigma)) + 
  ylab("") + 
  scale_fill_viridis(discrete=TRUE,end=0.8)
# ap = cowplot::align_plots(g1, p1, align="h",axis="tb")
# row_1 = cowplot::plot_grid(ap[[1]], ap[[2]],rel_widths = c(1, 1.6))

###############################################################
# second plot is for the cutoff of sigma probabilites
###############################################################
d <- dplyr::filter(df, sensitivity == "sigma_p")

d$sensitivity <- as.factor(d$sensitivity)
levels(d$sensitivity) <- c("")#c("Pr(sigma)")
d$prior.sigma_thresh <- as.factor(d$prior.sigma_thresh)

p2 <- ggplot(
  d,
  aes(n, dens_ll, col = prior.sigma_thresh)
) +
  geom_line(alpha=0.4) +
  #facet_wrap(~sensitivity, labeller = label_parsed) +
  xlab("Knots") +
  ylab("") +
  #labs(col = expression(paste("Pr(", sigma, ") < ", sigma[0]))) +
  labs(col = expression(theta[sigma])) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) + 
  scale_color_viridis(discrete=TRUE,end=0.8)

# now show the 3 priors 
r1 = sample_matern(matern_range = 20, matern_SD = 5, range_prob = 0.05, SD_prob = 0.01)
r2 = sample_matern(matern_range = 20, matern_SD = 5, range_prob = 0.05, SD_prob = 0.10)
r3 = sample_matern(matern_range = 20, matern_SD = 5, range_prob = 0.05, SD_prob = 0.40)
r1$sigma_prior = 0.01
r2$sigma_prior = 0.1
r3$sigma_prior = 0.4

g2 = ggplot(rbind(r1,r2,r3), aes(sigma, 
                                 fill=as.factor(sigma_prior),
                                 col=as.factor(sigma_prior))) + 
  geom_density(alpha=0.4,col=NA) + 
  theme_bw() + 
  theme(legend.position = "none") + 
  xlab(expression(sigma)) + 
  ylab("") + 
  scale_fill_viridis(discrete=TRUE,end=0.8)
# ap = cowplot::align_plots(g2, p2, align="h",axis="tb")
# row_2 = cowplot::plot_grid(ap[[1]], ap[[2]],rel_widths = c(1, 1.6))

############################################################
# Third plot is threshold for range / kappa
############################################################
d <- dplyr::filter(df, sensitivity == "range_cutoff")

d$sensitivity <- as.factor(d$sensitivity)
levels(d$sensitivity) <- c("kappa[0]")
d$prior.range_cutoff <- as.factor(d$prior.range_cutoff)

p3 <- ggplot(
  d,
  aes(n, dens_ll, col = prior.range_cutoff)
) +
  geom_line(alpha=0.4) +
  #facet_wrap(~sensitivity, labeller = label_parsed) +
  xlab("Knots") +
  ylab("") +
  labs(col = expression(paste(kappa[0]))) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) + 
  scale_color_viridis(discrete=TRUE,end=0.8)

# now show the 3 priors 
r1 = sample_matern(matern_range = 5, matern_SD = 5, range_prob = 0.05, SD_prob = 0.05)
r2 = sample_matern(matern_range = 20, matern_SD = 5, range_prob = 0.05, SD_prob = 0.05)
r3 = sample_matern(matern_range = 40, matern_SD = 5, range_prob = 0.05, SD_prob = 0.05)
r1$kappa = 5
r2$kappa = 20
r3$kappa = 40

g3 = ggplot(rbind(r1,r2,r3), aes(range, 
                                 fill=as.factor(kappa),
                                 col=as.factor(kappa))) + 
  geom_density(alpha=0.4,col=NA) + 
  theme_bw() + 
  theme(legend.position = "none") + 
  xlab(expression(kappa)) + 
  ylab("") + 
  scale_fill_viridis(discrete=TRUE,end=0.8)

# ap = cowplot::align_plots(g3, p3, align="h",axis="tb")
# row_3 = cowplot::plot_grid(ap[[1]], ap[[2]],rel_widths = c(1, 1.6))

############################################################
# Fourth plot is threshold for range / kappa
############################################################
d <- dplyr::filter(df, sensitivity == "range_p")

d$sensitivity <- as.factor(d$sensitivity)
levels(d$sensitivity) <- c("kappa")
d$range_p <- as.factor(d$prior.range_thresh)

p4 <- ggplot(
  d,
  aes(n, dens_ll, col = range_p)
) +
  geom_line(alpha=0.4) +
  #facet_wrap(~sensitivity, labeller = label_parsed) +
  xlab("Knots") +
  ylab("") +
  labs(col = expression(paste(theta[kappa]))) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) + 
  scale_color_viridis(discrete=TRUE,end=0.8)

# now show the 3 priors 
r1 = sample_matern(matern_range = 20, matern_SD = 5, range_prob = 0.01, SD_prob = 0.05)
r2 = sample_matern(matern_range = 20, matern_SD = 5, range_prob = 0.1, SD_prob = 0.05)
r3 = sample_matern(matern_range = 20, matern_SD = 5, range_prob = 0.4, SD_prob = 0.05)
r1$kappa_prob = 0.01
r2$kappa_prob = 0.1
r3$kappa_prob = 0.4

g4 = ggplot(rbind(r1,r2,r3), aes(range, 
                                 fill=as.factor(kappa_prob),
                                 col=as.factor(kappa_prob))) + 
  geom_density(alpha=0.4,col=NA) + 
  theme_bw() + 
  theme(legend.position = "none") + 
  xlab(expression(kappa)) + 
  ylab("") + 
  scale_fill_viridis(discrete=TRUE,end=0.8)

# ap = cowplot::align_plots(g4, p4, align="h",axis="tb")
# row_4 = cowplot::plot_grid(ap[[1]], ap[[2]],rel_widths = c(1, 1.6))

# final_fig = gridExtra::grid.arrange(row_1,row_2,row_3,row_4, nrow=4)
# annotate_figure(final_fig,
#                 left = text_grob("Tooth length", color = "green", rot = 90)
# )
#https://github.com/kassambara/ggpubr/issues/78
col_1 <- ggpubr::annotate_figure(gridExtra::arrangeGrob(g1,g2,g3,g4,ncol=1),
                                 left = text_grob("Prior density", rot = 90,vjust = 2)
)
col_2 <- ggpubr::annotate_figure(gridExtra::arrangeGrob(p1,p2,p3,p4,ncol=1),
                                 left = text_grob("Log density (test)", rot = 90, vjust=2)
)

final_fig <- gridExtra::arrangeGrob(col_1,col_2,ncol=2)
# 
# pdf("figures/Figure_S5.pdf")
# final_fig
# dev.off()
# 
# jpeg("figures/Figure_S5.jpeg")
# final_fig
# dev.off()
