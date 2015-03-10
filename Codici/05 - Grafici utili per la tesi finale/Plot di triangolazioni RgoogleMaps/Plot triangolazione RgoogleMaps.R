load("Territorio.RData")
library(RgoogleMaps)

xcom=x[is.na(Codici)==FALSE]
ycom=y[is.na(Codici)==FALSE]
xbound
ybound

xlim=c(min(x,xbound)-0.1,max(x,xbound)+0.1)
ylim=c(min(y,ybound)-0.1,max(y,ybound)+0.1)
Map<-GetMap.bbox(xlim, ylim, size = c(640, 640),MINIMUMSIZE = TRUE)

png(filename="Triangolazione.png")
PlotOnStaticMap(Map)
for (ne in 1:dim(Triang)[1])
{
    PlotOnStaticMap(Map,c(x[Triang[ne,1]],x[Triang[ne,2]]),lat=c(y[Triang[ne,1]],y[Triang[ne,2]]),col="red",type='l',lwd=1,add=TRUE)
    PlotOnStaticMap(Map,c(x[Triang[ne,2]],x[Triang[ne,3]]),lat=c(y[Triang[ne,2]],y[Triang[ne,3]]),col="red",type='l',lwd=1,add=TRUE)
    PlotOnStaticMap(Map,c(x[Triang[ne,3]],x[Triang[ne,1]]),lat=c(y[Triang[ne,3]],y[Triang[ne,1]]),col="red",type='l',lwd=1,add=TRUE)
}
dev.off()


#Zoom sul lido di venezia
xlim=c(12.2,12.5)
ylim=c(45.2,45.5)
Map<-GetMap.bbox(xlim, ylim, size = c(640, 640),MINIMUMSIZE = TRUE)

png(filename="Triangolazione zoom.png")
PlotOnStaticMap(Map)
for (ne in 1:dim(Triang)[1])
{
    PlotOnStaticMap(Map,c(x[Triang[ne,1]],x[Triang[ne,2]]),lat=c(y[Triang[ne,1]],y[Triang[ne,2]]),col="red",type='l',lwd=1,add=TRUE)
    PlotOnStaticMap(Map,c(x[Triang[ne,2]],x[Triang[ne,3]]),lat=c(y[Triang[ne,2]],y[Triang[ne,3]]),col="red",type='l',lwd=1,add=TRUE)
    PlotOnStaticMap(Map,c(x[Triang[ne,3]],x[Triang[ne,1]]),lat=c(y[Triang[ne,3]],y[Triang[ne,1]]),col="red",type='l',lwd=1,add=TRUE)
}
dev.off()