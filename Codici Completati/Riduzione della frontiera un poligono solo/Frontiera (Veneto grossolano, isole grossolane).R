#Creazione della frontiera in poligono unico
#Tolleranza per il Veneto: 0.05
#Tolleranza per le isole: 0.01
#Totale punti: 407
#AVVERTENZA! SE SI USA QUESTA FRONTIERA, UN COMUNE RESTA FUORI!
#SI TRATTA DI ARIANO NEL POLESINE
#CODICE 533

#Funzioni di appoggio
source("Functions.R")

##### DOWNLOAD DELLA FRONTIERA #####

#Scarico la frontiera
require(raster)
veneto =  subset(getData('GADM', country='ITA', level=1), NAME_1=="Veneto")
#plot(veneto)

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

## VENETO ##
#Ora devo studiare il caso del Veneto
#Il Veneto è nel poligono 1

dist=0
tol<-0.05

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

#Controllo se ci sono intersezioni
#Intersect<-Intersections(xVeneto,yVeneto)
#Intersect

plot(xVeneto,yVeneto,type='l',xlim=c(12.1,12.4),ylim=c(45,45.5))
#identify(xVeneto,yVeneto,pos=T,plot=T)

xVeneto=xVeneto[-c(182,185,186,192)]
yVeneto=yVeneto[-c(182,185,186,192)]

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
points(xVenezia,yVenezia,type="l")

#Controllo se ci sono intersezioni
#Intersect<-Intersections(xVenezia,yVenezia)
#Intersect

## CHIOGGIA ##
#Ora devo studiare il caso di Chioggia
#Chioggia è nel poligono 16 che va unito al poligono 4

dist=0
tol<-0.01

#POLIGONO 16
xChioggia_tot=xbound[labels==16]
yChioggia_tot=ybound[labels==16]

xChioggia1=xChioggia_tot[1]
yChioggia1=yChioggia_tot[1]

for(i in 2:length(xChioggia_tot))
{
    #Aggiorno la distanza
    dist<-dist+sqrt((xChioggia_tot[i]-xChioggia_tot[i-1])^2+(yChioggia_tot[i]-yChioggia_tot[i-1])^2)
    if(dist>tol)
    {
        xChioggia1<-c(xChioggia1,xChioggia_tot[i])
        yChioggia1<-c(yChioggia1,yChioggia_tot[i])
        dist<-0
    }
}
points(xChioggia1,yChioggia1,type="l")

#POLIGONO 4
xChioggia_tot=xbound[labels==4]
yChioggia_tot=ybound[labels==4]

xChioggia2=xChioggia_tot[1]
yChioggia2=yChioggia_tot[1]

for(i in 2:length(xChioggia_tot))
{
    #Aggiorno la distanza
    dist<-dist+sqrt((xChioggia_tot[i]-xChioggia_tot[i-1])^2+(yChioggia_tot[i]-yChioggia_tot[i-1])^2)
    if(dist>tol)
    {
        xChioggia2<-c(xChioggia2,xChioggia_tot[i])
        yChioggia2<-c(yChioggia2,yChioggia_tot[i])
        dist<-0
    }
}
points(xChioggia2,yChioggia2,type="l")

#Iniziamo ad intersecare questi due
plot(xChioggia1,yChioggia1,type="l",xlim=c(12.28,12.29),ylim=c(45.215,45.230))
points(xChioggia2,yChioggia2,type="l")

#identify(xChioggia1,yChioggia1)
#identify(xChioggia2,yChioggia2)

xChioggia<-c(xChioggia2[1:31],12.2875,12.2865,12.286,12.2841,12.284,12.283,xChioggia1[4:length(xChioggia1)],xChioggia1[1:2],12.2832,12.2877,xChioggia2[32:length(xChioggia2)])
yChioggia<-c(yChioggia2[1:31],45.217,45.217,45.2160,45.215,45.217,45.217,yChioggia1[4:length(yChioggia1)],yChioggia1[1:2],45.2175,45.2175,yChioggia2[32:length(yChioggia2)])
plot(xChioggia,yChioggia,type='l',xlim=c(12.28,12.29),ylim=c(45.215,45.230))

#Controllo se ci sono intersezioni
#Intersect<-Intersections(xChioggia,yChioggia)
#Intersect


##### UNIONE DELLE REGIONI #####
#Cerco di unire Venezia e CHioggia e di fare un poligono unico.
#Inizio da Venezia

#Ingrandisco nella zona di Venezia
plot(xVeneto,yVeneto,type="l",xlim=c(12.27,12.31),ylim=c(45.44,45.48))
points(xVenezia,yVenezia,type='l')

#identify(xVeneto,yVeneto)
#identify(xVenezia,yVenezia)

xBoundaries<-c(xVeneto[1:121],12.280411,xVenezia,12.2795,xVeneto[122:length(xVeneto)])
yBoundaries<-c(yVeneto[1:121],45.465299,yVenezia,45.464772,yVeneto[122:length(yVeneto)])

plot(xBoundaries,yBoundaries,type="l",xlim=c(12.27,12.31),ylim=c(45.44,45.48))

#Ingrandisco nella zona di Chioggia
plot(xBoundaries,yBoundaries,type="l",xlim=c(12.25,12.27),ylim=c(45.18,45.19))
points(xChioggia,yChioggia,type='l')

#identify(xBoundaries,yBoundaries)
#identify(xChioggia,yChioggia)

# I due poligoni si toccano, devo smussare in due punti
#Smusso nel primo punto
xBoundaries2<-c(xBoundaries[1:187],12.250,xBoundaries[188:length(xBoundaries)])
yBoundaries2<-c(yBoundaries[1:187],45.1818,yBoundaries[188:length(yBoundaries)])
    
plot(xBoundaries2,yBoundaries2,type="l",xlim=c(12.25,12.27),ylim=c(45.18,45.19))
points(xChioggia,yChioggia,type='l')

#Smusso di nuovo nel secondo punto
plot(xBoundaries2,yBoundaries2,type="l",xlim=c(12.26,12.32),ylim=c(45.176,45.184))
points(xChioggia,yChioggia,type='l')

#identify(xBoundaries2,yBoundaries2)
#identify(xChioggia,yChioggia)

xBoundaries3<-c(xBoundaries2[1:189],12.280,xBoundaries2[190:length(xBoundaries2)])
yBoundaries3<-c(yBoundaries2[1:189],45.177,yBoundaries2[190:length(yBoundaries2)])

plot(xBoundaries3,yBoundaries3,type="l",xlim=c(12.26,12.32),ylim=c(45.176,45.184))
points(xChioggia,yChioggia,type='l')

#Ora che ho smussato, sostituisco
xBoundaries<-xBoundaries3
yBoundaries<-yBoundaries3

#Ora devo unire a Chioggia
#Ingrandisco nella zona di Chioggia
plot(xBoundaries,yBoundaries,type="l",xlim=c(12.25,12.27),ylim=c(45.18,45.19))
points(xChioggia,yChioggia,type='l')

#identify(xBoundaries,yBoundaries)
#identify(xChioggia,yChioggia)

xBoundaries<-c(xBoundaries[1:188],12.251,xChioggia[19:length(xChioggia)],xChioggia[1:18],12.25203,12.2515,xBoundaries[189:length(xBoundaries)])
yBoundaries<-c(yBoundaries[1:188],45.182,yChioggia[19:length(yChioggia)],yChioggia[1:18],45.18244,45.182,yBoundaries[189:length(yBoundaries)])

plot(xBoundaries,yBoundaries,type="l",xlim=c(12.25,12.27),ylim=c(45.18,45.19))


##### CONTROLLI INTERSEZIONI E COMUNI #####
#Ho generato intersezioni?
#Intersect<-Intersections(xBoundaries,yBoundaries)
#Intersect
#Si, le ho generate
#Devo toglierle
plot(xBoundaries,yBoundaries,type="l",xlim=c(12.41,12.43),ylim=c(45.01,45.025))
#identify(xBoundaries,yBoundaries)
#Si risolve se tolgo il 267
plot(xBoundaries,yBoundaries,type="l",xlim=c(12.46,12.5),ylim=c(44.9,45))
#identify(xBoundaries,yBoundaries)
#Si risolve se tolgo il 273
xBoundaries<-xBoundaries[-c(267,273)]
yBoundaries<-yBoundaries[-c(267,273)]
Intersect<-Intersections(xBoundaries,yBoundaries)
Intersect
#Nessuna intersezione
plot(xBoundaries,yBoundaries,type="l")


#Tutti i comuni ora stanno all'interno della regione?
#Leggo i punti di comune
Coord<-read.table(file="Coordinate.txt",header=TRUE)
names(Coord)

#LA FRONTIERA RIDOTTA DESCRIVE UN POLIGONO CHE RACCHIUDE TUTTI I COMUNI?
library(SDMTools)
PolyPoints<-cbind(xBoundaries,yBoundaries)
if(sum(pnt.in.poly(cbind(Coord$Longitudine,Coord$Latitudine),PolyPoints)$pip)==length(Coord$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}
#Ci sono comuni esterni. Quali?

row<-pnt.in.poly(cbind(Coord$Longitudine,Coord$Latitudine),PolyPoints)$pip==0
Coord$Comune[row]
Coord$Codice[row]

##### SALVATAGGIO DEI RISULTATI #####
save(file="Frontiera (Veneto grossolano, isole grossolane).RData",xBoundaries,yBoundaries)
