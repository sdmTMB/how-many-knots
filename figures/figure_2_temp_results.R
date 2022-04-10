library(sf)
library(inlabru)
library(INLA)
library(dplyr)
library(fmesher)
library(tidyr)
library(ggplot2)
library(viridis)
library(raster)
options("rgdal_show_exportToProj4_warnings"="none")
library(rgdal)
library(scales)

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

# grab coastline
shore_gcs <- rnaturalearth::ne_countries(continent = "north america", scale = "medium", returnclass = "sp")
shore <- sp::spTransform(shore_gcs, CRS = CRS(SRS_string='EPSG:32610'))
shore <- fortify(shore)
unit_scale <- 1000 # to change units from m to km
shore$long <- shore$long/unit_scale
shore$lat <- shore$lat/unit_scale

haul$Temp = haul$temperature_at_gear_c_der# - mean(haul$temperature_at_gear_c_der,na.rm=T)
p1 <- dplyr::filter(haul, !is.na(Temp)) %>%
ggplot(aes(X,Y,col=Temp)) + 
  geom_point(alpha=0.7) + 
  xlab("Eastings") + 
  ylab("Northings") + 
  theme_bw() + 
  theme(panel.background=element_rect(fill = "grey20", colour = "grey20"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())+
  scale_color_gradient2(midpoint = mean(haul$temperature_at_gear_c_der,na.rm=T)) + 
  annotation_map(shore, color = "black", fill = "grey70",size=0.2)

p2 <- ggplot(df_long, aes(cutoff, value)) +  
  geom_point(col="darkblue",alpha=0.7) + geom_line(col="darkblue")+
  xlab("Cutoff distance (km)") + 
  ylab("Log density") + 
  facet_wrap(~Data,scale="free",ncol=1) + 
  theme_bw() + 
  theme(strip.background =element_rect(fill="white"))
  
g <- cowplot::plot_grid(p1,p2)

pdf("figures/Figure_02.pdf")
g
dev.off()

jpeg("figures/Figure_02.jpeg")
g
dev.off()