library(sf)
library(inlabru)
library(INLA)
library(dplyr)
library(fmesher)
library(tidyr)

d = readRDS("output/01_loocv_temp.rds")
d = dplyr::group_by(d,cutoff) %>%
  dplyr::summarise(tot_ll = sum(dens_ll,na.rm=T),
                   n = n[1])

# bring in cpo data
cpo = readRDS("output/01_temp_cpo.rds")
cpo = dplyr::filter(cpo,cutoff %in% d$cutoff)
d = dplyr::left_join(d, cpo[,c("cutoff","log_cpo")])

d = pivot_longer(d,cols=c(2,4))
d$metric = ifelse(d$name=="tot_ll","LOOCV predictive density","CPO")

g = ggplot(d, aes(cutoff, value)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~metric,scale="free",ncol=1) + 
  xlab("Cutoff distance (km)") + 
  ylab("Value") + 
  theme_bw() + 
  theme(strip.background =element_rect(fill="white"))



