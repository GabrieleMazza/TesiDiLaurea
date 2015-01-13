#DOWNLOAD DELLE MAPPE PER RGOOGLEMAPS
library(RgoogleMaps)

##### VENETO #####

x.range<-c(10.63369,13.08902)
y.range<-c(44.80249,46.66845)
#Li allargo
x<-c(10.53,13.18)
y<-c(44.7,46.76)

MapVeneto<-GetMap.bbox(lonR=x,latR=y,size = c(640,640),MINIMUMSIZE=TRUE)
png(filename="Veneto.png")
PlotOnStaticMap(MapVeneto)
dev.off()

##### VENEZIA #####

x.range<-c(11.99096,13.08623)
y.range<-c(45.31084,45.84894)
#Li allargo
x<-c(11.89,13.18)
y<-c(45.01,46.14)

MapVenezia<-GetMap.bbox(lonR=x,latR=y,size = c(640,640),MINIMUMSIZE=TRUE)
png(filename="Venezia.png")
PlotOnStaticMap(MapVenezia)
dev.off()

save(file="Maps.RData",MapVeneto,MapVenezia)
