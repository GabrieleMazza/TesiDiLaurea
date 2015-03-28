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
PlotOnStaticMap(Map,lon=xcom,lat=ycom,col="black",pch=16)
PlotOnStaticMap(Map,lon=xcom[Codici==414 | Codici==581],lat=ycom[Codici==414 | Codici==581],col="blue",pch=16,add=T)
PlotOnStaticMap(Map,lon=xbound,lat=ybound,col="red",type='l',lwd=1,add=TRUE)
dev.off()