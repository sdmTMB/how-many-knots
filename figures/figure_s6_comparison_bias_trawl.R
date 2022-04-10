library(sf)
library(inlabru)
library(INLA)
library(dplyr)
library(fmesher)

df = readRDS(file = "output/temp_model_df.rds")
all_df = readRDS(file = "output/temp_model_all_est.rds")

all_df = as.data.frame(all_df)
all_df$Data = ifelse(all_df$fold==2,"Train","Test")

g <- dplyr::filter(as.data.frame(all_df), fold==1) %>%
  ggplot(aes(cutoff, abs(pred-temperature_at_gear_c_der), group=trawl_id)) + 
  geom_line() + 
  facet_wrap(~trawl_id) + 
  xlab("Cutoff distance (km)") +
  ylab("Mean absolute error") + 
  theme_bw()

