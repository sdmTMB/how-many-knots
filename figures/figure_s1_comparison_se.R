library(sf)
library(inlabru)
library(INLA)
library(dplyr)
library(fmesher)

df = readRDS(file = "output/temp_model_df.rds")
all_df = readRDS(file = "output/temp_model_all_est.rds")

all_df = as.data.frame(all_df)
all_df$Data = ifelse(all_df$fold==2,"Train","Test")

g <- dplyr::group_by(all_df, cutoff, Data) %>% 
  dplyr::summarise(mean_se = mean(se)) %>%
  ggplot(aes(cutoff,mean_se,col=Data,group=Data)) + 
  geom_point() + geom_line()+
  theme_bw() + xlab("Cutoff distance (km)") + 
  ylab("Mean standard error")
