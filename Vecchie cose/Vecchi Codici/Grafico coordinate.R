#Genero il grafico per la visualizzazione dei dati del bordo del veneto e dei comuni
#interni. Inoltre estraggo tutti i punti con cui è definito il poligono di bordo


##### ESTRAZIONE DEI VERTICI DEL CONFINE #####

#Da Internet:
require(raster)
veneto =  subset(getData('GADM', country='ITA', level=1), NAME_1=="Veneto")
plot(veneto)

#Da RData
load("Veneto.RData")

#E' internamente memorizzato come poligono. Devo estrarne i vertici
#Esploriamo l'oggetto...
#slotNames(veneto)
#str(veneto)

tmp<- slot(veneto, 'polygons')
sub.tmp <- slot(tmp[[1]],'Polygons') 
order<-slot(tmp[[1]],'plotOrder')

xbound<-NULL
ybound<-NULL
labels<-NULL
for(i in 1:131)
{
  xbound<-c(xbound,sub.tmp[[order[i]]]@coords[,1])
  ybound<-c(ybound,sub.tmp[[order[i]]]@coords[,2])
  labels<-c(labels,rep(i,length(sub.tmp[[order[i]]]@coords[,1])))
}

length(xbound)
length(ybound)
length(labels)

save(file="Boundaries.RData",xbound,ybound,labels)

##### ESTRAZIONE DEI PUNTI DEI COMUNI #####

Data<-read.table("coordinates.txt",header=T)
attach(Data)
save(file="Comuni.RData",xcom,ycom,IDcom)
detach(Data)

##### GRAFICI #####

#Punti di confine
load("Boundaries.RData")
#Punti interni
load("Comuni.RData")

toplotx<-NULL
toploty<-NULL
for(i in 1:length(xbound))
{
  if(labels[i]==1)
  {
    toplotx<-c(toplotx,xbound[i])
    toploty<-c(toploty,ybound[i])
  }  
}

#Per il grafico completo:
plot(toplotx,toploty,xlim=c(10.62467, 13.10042), ylim=c(44.79208, 46.68047),main="Veneto", xlab="Longitudine", ylab="Latitudine", type='l')

#Per il grafico ingrandito
#plot(toplotx,toploty,xlim=c(12.2, 12.7), ylim=c(44.9, 45.6),main="Veneto", xlab="Longitudine", ylab="Latitudine", type='l')

points(xcom,ycom,pch=20,col="red")

#Venezia è nel poligono 6
#chioggia è nel poligono 16

for (j in 2:131)
{
  toplotx<-NULL
  toploty<-NULL
  for(i in 1:length(xbound))
  {
    if(labels[i]==j)
    {
      toplotx<-c(toplotx,xbound[i])
      toploty<-c(toploty,ybound[i])
    }  
  }
  points(toplotx,toploty,type='l')
}
