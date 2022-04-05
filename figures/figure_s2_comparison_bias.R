library(sf)
library(inlabru)
library(INLA)
library(dplyr)
library(fmesher)

df = readRDS(file = "output/temp_model_df.rds")
all_df = readRDS(file = "output/temp_model_all_est.rds")

all_df = as.data.frame(all_df)
all_df$Data = ifelse(all_df$fold==2,"Train","Test")


# make plot for cutoff distances of 40/80/120
g <- dplyr::group_by(all_df, cutoff, Data) %>% 
  dplyr::summarise(mean_pred = mean(abs(pred-temperature_at_gear_c_der))) %>%
  ggplot(aes(cutoff,mean_pred,col=Data,group=Data)) + 
  geom_point() + geom_line()+
  xlab("Cutoff distance (km)") + 
  ylab("Mean absolute error") + 
  #facet_wrap(~Data,scale="free")+
  theme_bw()
