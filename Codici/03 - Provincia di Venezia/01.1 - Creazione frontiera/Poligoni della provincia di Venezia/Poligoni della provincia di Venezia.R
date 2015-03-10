# POLIGONI DEL VENETO

#Scarico la frontiera
#require(raster)
#veneto =  subset(getData('GADM', country='ITA', level=1), NAME_1=="Veneto")
#plot(veneto)
#save(file="Veneto.RData",veneto)
load("Venezia.RData")

#La frontiera è memorizzata come una unione di 130 poligoni (a causa delle isole della)
#laguna veneta)
#Prendo solo la regione, ed è formata da un poligono di 20000 vertici.. lo riduco in base
#alla tolleranza indicata
#pIù bassa è la tolleranza, più punti sono tenuti per il poligono di confine

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

for(i in 1:101)
{
    png(filename=paste(i,".png",sep=""))
    plot(xbound[labels==i],ybound[labels==i],type='l',lwd=2,main=paste("Poligono ",i,", ",length(ybound[labels==i])," vertici",sep=""),xlab="Longitudine",ylab="Latitudine")
    dev.off()
}
