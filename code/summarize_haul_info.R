library(nwfscSurvey)
library(dplyr)
library(spatstat.geom)
library(sdmTMB)

d <- pull_haul()

d <- sdmTMB::add_utm_columns(d, ll_names = c("longitude_dd","latitude_dd"))

haul_ranges <- dplyr::group_by(d, year) |>
  dplyr::summarise(minX = min(X), maxX = max(X), minY = min(Y), maxY = max(Y)) |>
  as.data.frame()

nn_by_year <- d |>
  dplyr::group_by(year) |>
  dplyr::summarise(
    mean_nn_dist_km = if(n() > 1) mean(nndist(X,Y)) else NA_real_,
    n_points = n(),
    .groups = "drop"
  )

summaries <- dplyr::left_join(haul_ranges, nn_by_year)
