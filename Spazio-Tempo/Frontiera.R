##### DOWNLOAD DELLA FRONTIERA #####

#Scarico la frontiera
require(raster)
veneto =  subset(getData('GADM', country='ITA', level=1), NAME_1=="Veneto")
plot(veneto)

#La frontiera è memorizzata come una unione di 130 poligoni (a causa delle isole della)
#laguna veneta)
#Prendo solo la regione, ed è formata da un poligono di 20000 vertici.. lo riduco in base
#alla tolleranza indicata
#pIù bassa è la tolleranza, più punti sono tenuti per il poligono di confine

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

##### RIDUZIONE DELLA FRONTIERA #####

#Funzioni di appoggio
source("Functions.R")
#Punti di confine
#Leggo i punti di comune
Coord<-read.table(file="Coordinate.txt",header=TRUE)
names(Coord)


## VENETO ##
#Ora devo studiare il caso del Veneto
#Il Veneto è nel poligono 1

dist=0
tol<-0.005

xVeneto_tot=xbound[labels==1]
yVeneto_tot=ybound[labels==1]

xVeneto=xVeneto_tot[1]
yVeneto=yVeneto_tot[1]

for(i in 2:length(xVeneto_tot))
{
    #Aggiorno la distanza
    dist<-dist+sqrt((xVeneto_tot[i]-xVeneto_tot[i-1])^2+(yVeneto_tot[i]-yVeneto_tot[i-1])^2)
    if(dist>tol)
    {
        xVeneto<-c(xVeneto,xbound[i])
        yVeneto<-c(yVeneto,ybound[i])
        dist<-0
    }  
}
#NUMERO DI PUNTI CON CUI DESCRIVO LA FRONTIERA DEL VENETO
length(xVeneto)

#Controllo se ci sono intersezioni
#Intersect<-Intersections(xVeneto,yVeneto)
#Intersect
#Se sono segnalate intersezioni, allora devo intervenire, togliendo alcuni punti

plot(xVeneto,yVeneto,type='l',xlim=c(12.1,12.4),ylim=c(45,45.5))
points(Coord$Longitudine,Coord$Latitudine,pch=16,col="red")
#identify(xVeneto,yVeneto,pos=T,plot=T)

#Casi analizzati
#TOLLERANZA 0.05
IDdelete<-c(182,185,186,192)
#TOLLERANZA 0.005
IDdelete<-c(1714,1715)

xVeneto=xVeneto[-IDdelete]
yVeneto=yVeneto[-IDdelete]

#Tutti i comuni per ora stanno all'interno della regione?
#Devo provare tutto tranne venezia e chioggia
IDdelete<-c(391,425)
ROWdelete<-NULL

for (i in 1:length(Coord$Codice))
{
    for (j in 1:length(IDdelete))
    {
        if (Coord$Codice[i]==IDdelete[j])
        {
            ROWdelete<-c(ROWdelete,i)
        }
    }
}

tmpLon<-Coord$Longitudine[-ROWdelete]
tmpLat<-Coord$Latitudine[-ROWdelete]

#LA FRONTIERA RIDOTTA DESCRIVE UN POLIGONO CHE RACCHIUDE TUTTI I COMUNI?
library(SDMTools)
PolyPoints<-cbind(xVeneto,yVeneto)
if(sum(pnt.in.poly(cbind(tmpLon,tmpLat),PolyPoints)$pip)==length(tmpLon))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}

#Se qualcuno non sta dentro, allora devo abbassare la tolleranza


## VENEZIA ##
#Ora devo studiare il caso di Venezia
#Venezia è nel poligono 6

dist=0
tol<-0.01

xVenezia_tot=xbound[labels==6]
yVenezia_tot=ybound[labels==6]

xVenezia=xVenezia_tot[1]
yVenezia=yVenezia_tot[1]

for(i in 2:length(xVenezia_tot))
{
    #Aggiorno la distanza
    dist<-dist+sqrt((xVenezia_tot[i]-xVenezia_tot[i-1])^2+(yVenezia_tot[i]-yVenezia_tot[i-1])^2)
    if(dist>tol)
    {
        xVenezia<-c(xVenezia,xVenezia_tot[i])
        yVenezia<-c(yVenezia,yVenezia_tot[i])
        dist<-0
    }
}
length(xVenezia)
points(xVenezia,yVenezia,type="l")

#Controllo se ci sono intersezioni
Intersect<-Intersections(xVenezia,yVenezia)
Intersect

#Controllo che Venezia stia dentro
library(SDMTools)
PolyPoints<-cbind(xVenezia,yVenezia)
if(pnt.in.poly(cbind(12.327500,45.438056),PolyPoints)$pip==1)
{
    print("Venezia è dentro")
} else
{
    print("Venezia è fuori")
}


## CHIOGGIA ##
#Ora devo studiare il caso di Chioggia
#Chioggia è nel poligono 16

dist=0
tol<-0.01

xChioggia_tot=xbound[labels==16]
yChioggia_tot=ybound[labels==16]

xChioggia=xChioggia_tot[1]
yChioggia=yChioggia_tot[1]

for(i in 2:length(xChioggia_tot))
{
    #Aggiorno la distanza
    dist<-dist+sqrt((xChioggia_tot[i]-xChioggia_tot[i-1])^2+(yChioggia_tot[i]-yChioggia_tot[i-1])^2)
    if(dist>tol)
    {
        xChioggia<-c(xChioggia,xChioggia_tot[i])
        yChioggia<-c(yChioggia,yChioggia_tot[i])
        dist<-0
    }
}
length(xChioggia)
points(xChioggia,yChioggia,type="l")

#Controllo se ci sono intersezioni
Intersect<-Intersections(xChioggia,yChioggia)
Intersect

#Controllo che Venezia stia dentro
library(SDMTools)
PolyPoints<-cbind(xChioggia,yChioggia)
if(pnt.in.poly(cbind(12.279444,45.220556),PolyPoints)$pip==1)
{
    print("Chioggia è dentro")
} else
{
    print("Chioggia è fuori")
}



##### SALVATAGGIO DEI RISULTATI #####
save(file="Frontiera.RData",xVeneto,yVeneto,xVenezia,yVenezia,xChioggia,yChioggia)
