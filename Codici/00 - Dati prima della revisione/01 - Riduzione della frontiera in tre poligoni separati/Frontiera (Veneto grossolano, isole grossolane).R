#Creazione della frontiera in tre poligoni
#Tolleranza per il Veneto: 0.05
#Tolleranza per le isole: 0.01
#Totale punti Veneto: 330
#Totale punti Venezia: 27
#Totale punti Chioggia: 47
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


##### MODIFICHE DELLE REGIONI #####

#Ingrandisco nella zona di Chioggia
plot(xVeneto,yVeneto,type="l",xlim=c(12.25,12.27),ylim=c(45.18,45.19))
points(xChioggia,yChioggia,type='l')

#identify(xVeneto,yVeneto)
#identify(xChioggia,yChioggia)

# I due poligoni si toccano, devo smussare in due punti
#Smusso nel primo punto
xVeneto2<-c(xVeneto[1:158],12.250,xVeneto[159:length(xVeneto)])
yVeneto2<-c(yVeneto[1:158],45.1818,yVeneto[159:length(yVeneto)])
    
plot(xVeneto2,yVeneto2,type="l",xlim=c(12.25,12.27),ylim=c(45.18,45.19))
points(xChioggia,yChioggia,type='l')

#Smusso di nuovo nel secondo punto
plot(xVeneto2,yVeneto2,type="l",xlim=c(12.26,12.32),ylim=c(45.176,45.184))
points(xChioggia,yChioggia,type='l')

#identify(xVeneto2,yVeneto2)
#identify(xChioggia,yChioggia)

xVeneto3<-c(xVeneto2[1:160],12.280,xVeneto2[161:length(xVeneto2)])
yVeneto3<-c(yVeneto2[1:160],45.177,yVeneto2[161:length(yVeneto2)])

plot(xVeneto3,yVeneto3,type="l",xlim=c(12.26,12.32),ylim=c(45.176,45.184))
points(xChioggia,yChioggia,type='l')

#Ora che ho smussato, sostituisco
xVeneto<-xVeneto3
yVeneto<-yVeneto3

##### SALVATAGGIO DEI RISULTATI #####
save(file="Frontiera (Veneto grossolano, isole grossolane).RData",xVeneto,yVeneto,xVenezia,yVenezia,xChioggia,yChioggia)
