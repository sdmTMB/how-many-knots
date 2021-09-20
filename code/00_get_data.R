library(sf)
library(inlabru)
library(INLA)
library(dplyr)
library(sp)

# Prepare data
# haul data includes environmental covariates with location information
haul <- nwfscSurvey::PullHaul.fn(
  YearRange = c(2018),
  SurveyName = "NWFSC.Combo"
)

# project lat/lon to UTM, after removing missing values and unsatisfactory hauls
haul <- haul %>%
  filter(
    !is.na(longitude_dd), !is.na(latitude_dd),
    performance == "Satisfactory"
  ) %>%
  dplyr::select(trawl_id, latitude_dd, longitude_dd, depth_hi_prec_m, temperature_at_gear_c_der)

haul_trans <- haul
coordinates(haul_trans) <- c("longitude_dd", "latitude_dd")
proj4string(haul_trans) <- CRS("+proj=longlat +datum=WGS84")
newproj <- paste("+proj=utm +zone=10 ellps=WGS84 +datum=WGS84")

haul_trans <- spTransform(haul_trans, CRS(newproj))
haul_trans <- as.data.frame(haul_trans)
haul$X <- haul_trans$longitude_dd
haul$Y <- haul_trans$latitude_dd

# center and scale depth, removing NAs
haul <- dplyr::filter(haul, !is.na(depth_hi_prec_m))
haul$log_depth_scaled <- scale(log(haul$depth_hi_prec_m))
haul$log_depth_scaled2 <- haul$log_depth_scaled^2

# catch data includes catch, effort, etc. This takes a few minutes to grab all ~ 900 spp
catch <- nwfscSurvey::PullCatch.fn(
  YearRange = c(2018),
  SurveyName = "NWFSC.Combo"
)
# format to later join catch and haul
names(catch) <- tolower(names(catch))
catch$trawl_id <- as.numeric(catch$trawl_id)

dover <- dplyr::filter(catch, common_name == "Dover sole")

# Join catch and haul data
haul_new <- haul %>%
  left_join(dover, by = "trawl_id") %>%
  select(trawl_id, X, Y,
    latitude = latitude_dd.x,
    longitude = longitude_dd.x,
    year = year,
    log_depth_scaled,
    log_depth_scaled2,
    cpue_kg_km2,
    temperature_at_gear_c_der
  )
# Set NA CPUEs to 0
haul_new$cpue_kg_km2[which(is.na(haul_new$cpue_kg_km2))] <- 0

saveRDS(haul_new, "data/doversole_cleaned.rds")

saveRDS(haul, "data/haul_cleaned.rds")

# petrale sole", "darkblotched rockfish","lingcod", "sablefish
sub <- dplyr::filter(catch, common_name %in% c(
  "Dover sole",
  "petrale sole",
  "darkblotched rockfish",
  "lingcod"
))
saveRDS(sub, "data/catch_cleaned.rds")



# top_spp <-
#   dplyr::group_by(catch, common_name) %>%
#   dplyr::summarise(m = sum(cpue_kg_km2)) %>%
#   dplyr::arrange(-m) %>%
#   dplyr::filter(!is.na(common_name))
#
# sub = dplyr::filter(catch, common_name %in% top_spp$common_name[1:30]) %>%
#   dplyr::filter(common_name %in% c("crushed urchin","Mud Urchin","Red Star")==FALSE)

# saveRDS(sub, "data/catch_cleaned_top50.rds")
