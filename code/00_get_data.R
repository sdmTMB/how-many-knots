remotes::install_github("pbs-assess/sdmTMB", "b146eb3")
remotes::install_github("pfmc-assessments/nwfscSurvey", "0386065")

library(sf)
library(inlabru)
library(INLA)
library(dplyr)
library(sp)

# Prepare data
# haul data includes environmental covariates with location information
haul <- nwfscSurvey::pull_haul(
  #YearRange = c(2018),
  survey = "NWFSC.Combo"
)

# project lat/lon to UTM, after removing missing values and unsatisfactory hauls
haul <- haul %>%
  filter(
    !is.na(longitude_dd), !is.na(latitude_dd),
    performance == "Satisfactory"
  ) %>%
  dplyr::select(trawl_id, latitude_dd, longitude_dd, depth_hi_prec_m, temperature_at_gear_c_der, date_formatted)

haul <- sdmTMB::add_utm_columns(haul, ll_names = c("longitude_dd","latitude_dd"))

# center and scale depth, removing NAs
haul <- dplyr::filter(haul, !is.na(depth_hi_prec_m))
haul$log_depth_scaled <- scale(log(haul$depth_hi_prec_m))
haul$log_depth_scaled2 <- haul$log_depth_scaled^2

# catch data includes catch, effort, etc. This takes a few minutes to grab all ~ 900 spp
catch <- nwfscSurvey::pull_catch(
  #YearRange = c(2018),
  survey = "NWFSC.Combo"
)
# format to later join catch and haul
names(catch) <- tolower(names(catch))
catch$trawl_id <- as.numeric(catch$trawl_id)
catch <- sdmTMB::add_utm_columns(catch, ll_names = c("longitude_dd","latitude_dd"))

saveRDS(haul, "data/haul_cleaned.rds")

dplyr::group_by(catch, common_name) %>% 
  dplyr::summarize(n = length(which(cpue_kg_km2>0))) %>%
  dplyr::arrange(-n)

# get url from 
url <- "https://raw.githubusercontent.com/pfmc-assessments/indexwc/refs/heads/main/data-raw/configuration.csv"
config <- read.csv(url)

catch$common_name <- tolower(catch$common_name)

# petrale sole", "darkblotched rockfish","lingcod", "sablefish
sub <- dplyr::filter(catch, common_name %in% config$species)
names(sub) <- tolower(names(sub))
saveRDS(sub, "data/catch_cleaned.rds")

#spp_names <- dplyr::group_by(sub, common_name) |>
#  dplyr::summarise(sci_name = scientific_name[1])
#write.csv(spp_names, "species_names.csv", row.names = FALSE)

# create temperature dataset
trawlid_date <- dplyr::filter(sub, year==2018) |>
  dplyr::group_by(trawl_id) %>%
  dplyr::summarize(date = date[1])
trawlid_date$month <- substr(trawlid_date$date, 6, 8)
trawlid_date <- dplyr::filter(trawlid_date, month == "Jun")
haul <- dplyr::filter(haul, trawl_id %in% trawlid_date$trawl_id)
saveRDS(haul, "data/june_bottom_temp.rds")

