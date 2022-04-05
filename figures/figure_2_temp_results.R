library(sf)
library(inlabru)
library(INLA)
library(dplyr)
library(fmesher)
library(tidyr)
df = readRDS(file = "output/temp_model_df.rds")
all_df = readRDS(file = "output/temp_model_all_est.rds")

all_df = as.data.frame(all_df)
all_df$Data = ifelse(all_df$fold==2,"train","test")

haul <- readRDS("data/june_bottom_temp.rds")
haul[, c("X", "Y")] <- haul[, c("X", "Y")] / 1000

meshes <- readRDS(file = "output/temp_meshes.rds")

# plot predictive in-sample density
df_long <- pivot_longer(df, cols = 4:5)
df_long$Data = ifelse(df_long$name=="dens_ll_train","Train","Test")

haul$Temp = haul$temperature_at_gear_c_der - mean(haul$temperature_at_gear_c_der,na.rm=T)
p1 <- dplyr::filter(haul, !is.na(Temp)) %>%
ggplot(aes(X,Y,col=Temp)) + 
  geom_point(alpha=0.7) + 
  xlab("Eastings") + 
  ylab("Northings") + 
  theme_bw() + 
  scale_color_gradient2()

p2 <- ggplot(df_long, aes(cutoff, value)) +  
  geom_point(col="darkblue",alpha=0.7) + geom_line(col="darkblue")+
  xlab("Cutoff distance (km)") + 
  ylab("Log predictive density") + 
  facet_wrap(~Data,scale="free",ncol=1) + 
  theme_bw() + 
  theme(strip.background =element_rect(fill="white"))
  
g <- cowplot::plot_grid(p1,p2)

