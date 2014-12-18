#Download delle mappe utili per RgoogleMaps e i suoi grafici
library(RgoogleMaps)

#Secondo lo smoothing principale
#> range(xVeneto)
#[1] 10.63369 13.08902
#> range(yVeneto)
#[1] 44.80249 46.66845

RangeVenetoLon<-c(10.58, 13.13)
RangeVenetoLat<-c(44.75, 46.71)
MapVeneto<-GetMap.bbox(lonR=RangeVenetoLon, latR=RangeVenetoLat, size = c(640,640),MINIMUMSIZE=T)
PlotOnStaticMap(MapVeneto)

save(file="MapVeneto.RData",MapVeneto)
