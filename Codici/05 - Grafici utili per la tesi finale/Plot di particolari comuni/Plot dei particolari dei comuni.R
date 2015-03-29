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
PlotOnStaticMap(Map,lon=xcom[Codici==414 | Codici==581],lat=ycom[Codici==414 | Codici==581],col="blue",pch=16,cex=2,add=T)
PlotOnStaticMap(Map,lon=xbound,lat=ybound,col="red",type='l',lwd=1,add=TRUE)
dev.off()


png(filename="Comuni replicati2.png")
PlotOnStaticMap(Map,lon=xcom,lat=ycom,col="black",pch=16)
PlotOnStaticMap(Map,lon=xcom[Codici==581],lat=ycom[Codici==581],col="blue",pch=16,cex=2,add=T)
PlotOnStaticMap(Map,lon=xcom[Codici==586],lat=ycom[Codici==586],col="green",pch=16,cex=2,add=T)
PlotOnStaticMap(Map,lon=xbound,lat=ybound,col="black",type='l',lwd=1,add=TRUE)
dev.off()

png(filename="Selezione.png")
PlotOnStaticMap(Map,lon=xcom[Codici==581 | Codici==586 | Codici==389 | Codici==425 | Codici==388 | Codici==402],lat=ycom[Codici==581 | Codici==586 | Codici==389 | Codici==425 | Codici==388 | Codici==402],col="blue",pch=16,cex=2)
PlotOnStaticMap(Map,lon=xbound,lat=ybound,col="black",type='l',lwd=1,add=TRUE)
dev.off()
#Zoom sul lido di venezia
xlim=c(12.2,12.5)
ylim=c(45.2,45.5)
Map<-GetMap.bbox(xlim, ylim, size = c(640, 640),MINIMUMSIZE = TRUE)


png(filename="Comuni replicati.png")
PlotOnStaticMap(Map,lon=xcom,lat=ycom,col="black",pch=16)
PlotOnStaticMap(Map,lon=xcom[Codici==582 | Codici==583 | Codici==584 | Codici==585],lat=ycom[Codici==582 | Codici==583 | Codici==584 | Codici==585],col="red",pch=16,cex=2,add=T)
PlotOnStaticMap(Map,lon=xbound,lat=ybound,col="black",type='l',lwd=1,add=TRUE)
dev.off()
