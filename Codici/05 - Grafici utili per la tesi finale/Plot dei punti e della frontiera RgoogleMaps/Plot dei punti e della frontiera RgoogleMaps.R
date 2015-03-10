load("Territorio.RData")
library(RgoogleMaps)

xcom=x[is.na(Codici)==FALSE]
ycom=y[is.na(Codici)==FALSE]
xbound
ybound

xlim=c(min(x,xbound)-0.1,max(x,xbound)+0.1)
ylim=c(min(y,ybound)-0.1,max(y,ybound)+0.1)
Map<-GetMap.bbox(xlim, ylim, size = c(640, 640),MINIMUMSIZE = TRUE)

png(filename="Comuni e frontiera.png")
PlotOnStaticMap(Map,lon=xcom,lat=ycom,col="red",pch=20)
PlotOnStaticMap(Map,lon=xbound,lat=ybound,col="red",type='l',lwd=1,add=TRUE)
dev.off()

#Zoom sul lido di venezia
xlim=c(12.2,12.5)
ylim=c(45.2,45.5)
Map<-GetMap.bbox(xlim, ylim, size = c(640, 640),MINIMUMSIZE = TRUE)

png(filename="Comuni e frontiera zoom.png")
PlotOnStaticMap(Map,lon=xcom,lat=ycom,col="red",pch=16)
PlotOnStaticMap(Map,lon=xbound,lat=ybound,col="red",type='l',add=TRUE)
dev.off()