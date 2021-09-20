library(ggplot2)
library(sdmTMB)
library(dplyr)
library(future)
library(sf)
library(inlabru)
library(INLA)

x <- seq(-1, 1, length.out = 150)
y <- seq(-1, 1, length.out = 150)
loc <- expand.grid(x = x, y = y)
time_steps <- 1
N <- nrow(loc)
X <- model.matrix(~1, data.frame(x1 = rnorm(N * time_steps)))
mesh <- make_mesh(loc, xy_cols = c("x", "y"), cutoff = 0.1)
plot(mesh$mesh, asp = 1, main = "")
s <- sdmTMB_sim(
  x = loc$x, y = loc$y, mesh = mesh, X = X,
  betas = c(0.2, 0.3), time_steps = time_steps,
  phi = 0.2, thetaf = 1.3, range = 0.2, sigma_O = 0.4,
  seed = 1234, family = gaussian()
)

true <-
  ggplot(s, aes(x, y, fill = mu)) +
  geom_raster() +
  scale_fill_viridis_c()

# split

# Sample 500 points to observe:
set.seed(1029)
s$fold <- sample(1:2, size = nrow(s), replace = T, prob = c(0.9, 0.1))
sampled <- sample(seq_len(nrow(loc)), 500)

# true_b <- ggplot(s, aes(x, y, fill = mu))+
#   geom_raster() +
#   scale_fill_viridis_c() +
#   geom_point(data=d, aes(x,y), shape=1, alpha=0.5)

# Fit an cross validation across a sequence of INLA mesh
# using 90/10 split `cutoff` values i.e.,
# smaller cutoff = higher resolution mesh:

# assign folds randomly, using a 90/10 split for train/test
cutoffs <- c(0.35, 0.18, 0.1, 0.056, 0.015)
test_ll <- 0

for (i in 4:length(cutoffs)) {
  d <- s[sampled, ]
  # return to SpatialPointsDataFrame
  coordinates(d) <- c("x", "y")
  mesh <- inla.mesh.2d(
    loc = coordinates(d),
    cutoff = cutoffs[i], max.n = 1000
  )

  matern <-
    inla.spde2.pcmatern(mesh,
      prior.sigma = c(2, 0.05),
      prior.range = c(3, 0.05)
    )

  # components is equivalent to formula
  components <- observed ~ Intercept + field(
    main = coordinates,
    model = matern
  )

  # do cross validation here

  fit_train <- try(bru(components,
    d[which(d$fold == 1), ],
    family = "gaussian"
  ), silent = TRUE)

  # calculate LL for out of sample data
  pred_test <- predict(fit_train,
    data = d[which(d$fold == 2), ],
    formula = ~ Intercept + field
  )

  tau <- fit_train$summary.hyperpar$mean[1]
  test_ll[i] <- sum(dnorm(d$observed[which(d$fold == 2)],
    mean = pred_test$mean,
    sd = sqrt(1 / tau),
    log = TRUE
  ))

  # predict the whole df for pred-obs plots
  pred_test <- predict(fit_train,
    data = d,
    formula = ~ Intercept + field
  )
  d$pred <- pred_test$mean

  # also predict a field for mapping
  nd <- expand.grid(
    x = seq(-1, 1, length.out = 150),
    y = seq(-1, 1, length.out = 150)
  )
  pred_test <- predict(fit_train,
    data = nd,
    formula = ~ Intercept + field
  )
  nd$pred <- pred_test$mean

  nd$cutoff <- d$cutoff <- cutoffs[i]
  nd$knots <- d$knots <- mesh$n

  if (i == 1) {
    df_all <- d
    nd_all <- nd
  } else {
    df_all <- rbind(df_all, d)
    nd_all <- rbind(nd_all, nd)
  }
}



# Predictions from 1st fold:

ggplot(p, aes(x, y, fill = est)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~meta, nrow = 1) +
  coord_fixed(expand = FALSE) +
  ggtitle("Predicted value") +
  theme_bw() +
  xlab("Longitude") +
  ylab("Latitude")
