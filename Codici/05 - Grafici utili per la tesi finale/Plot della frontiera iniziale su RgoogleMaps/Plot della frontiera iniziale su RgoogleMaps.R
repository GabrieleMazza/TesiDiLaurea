load("Venezia.RData")
library(RgoogleMaps)

tmp<- slot(venezia, 'polygons')
sub.tmp <- slot(tmp[[1]],'Polygons') 
order<-slot(tmp[[1]],'plotOrder')

xbound<-NULL
ybound<-NULL
labels<-NULL
for(i in 1:101)
{
    xbound<-c(xbound,sub.tmp[[order[i]]]@coords[,1])
    ybound<-c(ybound,sub.tmp[[order[i]]]@coords[,2])
    labels<-c(labels,rep(i,length(sub.tmp[[order[i]]]@coords[,1])))
}

length(xbound)
length(ybound)
length(labels)

xlim=c(min(xbound)-0.1,max(xbound)+0.1)
ylim=c(min(ybound)-0.1,max(ybound)+0.1)
Map<-GetMap.bbox(xlim, ylim, size = c(640, 640),MINIMUMSIZE = TRUE)

png(filename="Frontiera iniziale.png")
PlotOnStaticMap(Map)
for(i in 1:101)
{
    x<-xbound[labels==i]
    y<-ybound[labels==i]
    for(j in 1:length(x))
    {
        PlotOnStaticMap(Map,lon=c(x[j],x[j+1]),lat=c(y[j],y[j+1]),col="red",type='l',lwd=1,add=TRUE)
    }
    PlotOnStaticMap(Map,lon=c(x[j],x[1]),lat=c(y[j],y[1]),col="red",type='l',lwd=1,add=TRUE)
}
dev.off()

png(filename="Frontiera iniziale ridotta.png")
PlotOnStaticMap(Map)
for(i in c(1,2,5,6,10,15,4,11))
{
    x<-xbound[labels==i]
    y<-ybound[labels==i]
    for(j in 1:length(x))
    {
        PlotOnStaticMap(Map,lon=c(x[j],x[j+1]),lat=c(y[j],y[j+1]),col="red",type='l',lwd=1,add=TRUE)
    }
    PlotOnStaticMap(Map,lon=c(x[j],x[1]),lat=c(y[j],y[1]),col="red",type='l',lwd=1,add=TRUE)
}
dev.off()


#Zoom sul lido di venezia
xlim=c(12.2,12.5)
ylim=c(45.2,45.5)
Map<-GetMap.bbox(xlim, ylim, size = c(640, 640),MINIMUMSIZE = TRUE)

png(filename="Frontiera iniziale zoom.png")
PlotOnStaticMap(Map)
for(i in 1:101)
{
    x<-xbound[labels==i]
    y<-ybound[labels==i]
    for(j in 1:length(x))
    {
        PlotOnStaticMap(Map,lon=c(x[j],x[j+1]),lat=c(y[j],y[j+1]),col="red",type='l',lwd=1,add=TRUE)
    }
    PlotOnStaticMap(Map,lon=c(x[j],x[1]),lat=c(y[j],y[1]),col="red",type='l',lwd=1,add=TRUE)
}
dev.off()

png(filename="Frontiera iniziale ridotta zoom.png")
PlotOnStaticMap(Map)
for(i in c(1,2,5,6,10,15,4,11))
{
    x<-xbound[labels==i]
    y<-ybound[labels==i]
    for(j in 1:length(x))
    {
        PlotOnStaticMap(Map,lon=c(x[j],x[j+1]),lat=c(y[j],y[j+1]),col="red",type='l',lwd=1,add=TRUE)
    }
    PlotOnStaticMap(Map,lon=c(x[j],x[1]),lat=c(y[j],y[1]),col="red",type='l',lwd=1,add=TRUE)
}
dev.off()