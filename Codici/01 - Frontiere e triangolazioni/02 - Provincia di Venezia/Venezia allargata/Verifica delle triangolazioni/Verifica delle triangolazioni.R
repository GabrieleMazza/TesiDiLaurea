#DOWNLOAD DELLE MAPPE PER RGOOGLEMAPS
library(RgoogleMaps)
load("Prova.RData")
##### VENETO #####

x.r<-c(12.28,12.4)
y.r<-c(45.4,45.46)

Map<-GetMap.bbox(lonR=x.r,latR=y.r,size = c(640,640),MINIMUMSIZE=TRUE)
PlotOnStaticMap(Map,lon=xVeneto,lat=yVeneto,type='l')

x.r<-c(12.0,12.6)
y.r<-c(45.0,45.5)
Map<-GetMap.bbox(lonR=x.r,latR=y.r,size = c(640,640),MINIMUMSIZE=TRUE)
PlotOnStaticMap(Map,lon=xVeneto,lat=yVeneto,type='l')
