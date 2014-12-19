#Creazione della frontiera in poligono unico usamdo una tecnica di smoothing
#La provincia di 



#Funzioni di appoggio
source("Functions.R")
library(fda)
library(KernSmooth)
library(geosphere)
library(RTriangle)
library(SDMTools)




##### DOWNLOAD DELLA FRONTIERA #####

#Scarico la frontiera
#require(raster)
#venezia =  subset(getData('GADM', country='ITA', level=2), NAME_2=="Venezia")
#plot(venezia)
#save(file="Venezia.RData",venezia)
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
for(i in 1:131)
{
    xbound<-c(xbound,sub.tmp[[order[i]]]@coords[,1])
    ybound<-c(ybound,sub.tmp[[order[i]]]@coords[,2])
    labels<-c(labels,rep(i,length(sub.tmp[[order[i]]]@coords[,1])))
}

length(xbound)
length(ybound)
length(labels)

#Leggo anche i comuni (solo della provincia di venezia)
Comuni<-read.table(file="Coordinate.txt",header=T)



##### VENETO #####

#Ora devo studiare il caso del Veneto
#Il Veneto è nel poligono 1
#Studio tutti i dati in dipendenza dall'ascissa curvilinea


xVenezia_tot=xbound[labels==1]
yVenezia_tot=ybound[labels==1]
s<-0
previous<-0

for(i in 2:length(xVenezia_tot))
{
    #Aggiorno la distanza
    previous<-previous+sqrt((xVenezia_tot[i]-xVenezia_tot[i-1])^2+(yVenezia_tot[i]-yVenezia_tot[i-1])^2)
    s<-c(s,previous)  
}

#Ora devo fare smoothing di queste due funzioni
#Come si fa?

#Quanti punti fisso per la definizione del confine?
N<-100
#Nel range di ascissa curvilinea fisso N punti
Val<-seq(0,max(s),by=max(s)/N)
#Decido di usare delle splines cubiche
m<-4


##### REGRESSION SPLINES #####

nbasis<-90
basis = create.bspline.basis(c(0,max(s)), nbasis, m)
xsmooth = smooth.basis(argvals=s, y=xVenezia_tot, basis)
ysmooth = smooth.basis(argvals=s, y=yVenezia_tot, basis)
xVenezia  = eval.fd(Val, xsmooth$fd)
yVenezia  = eval.fd(Val, ysmooth$fd)
#Traccio ora i risultati
png(filename = "Miglior smoothing in x con regression splines.png")
plot(s,xVenezia_tot,main="Smoothing di xVenezia_tot")
points(Val,xVenezia,col="blue",type='l')
dev.off()
png(filename = "Miglior smoothing in y con regression splines.png")
plot(s,yVenezia_tot,main="Smoothing di yVenezia_tot")
points(Val,yVenezia,col="blue",type='l')
dev.off()
png(filename = "Miglior smoothing della laguna con regression splines.png")
plot(xVenezia_tot,yVenezia_tot,type='l',main="Nuova definizione della regione",xlim=c(12.1,12.6),ylim=c(45,45.6))
points(xVenezia,yVenezia,col="blue",type='l')
dev.off()
png(filename = "Miglior smoothing della regione con regression splines.png")
plot(xVenezia_tot,yVenezia_tot,type='l',main="Nuova definizione della regione")
points(xVenezia,yVenezia,col="blue",type='l')
dev.off()

#Controllo di intersezioni
Intersect<-Intersections(xVenezia,yVenezia)
Intersect

#Oggetti completi
x<-c(Comuni$Longitudine,xVenezia)
y<-c(Comuni$Latitudine,yVenezia)
#Creo i Boundaries
Boundaries<-NULL
for(i in (length(Comuni$Longitudine)+1):(length(x)-1))
{
    Boundaries<-rbind(Boundaries, c(i,i+1))
}
Boundaries<-rbind(Boundaries, c(length(x),length(Comuni$Longitudine)+1))
#Ora triangolazione
#Oggetto pslg
pslg_obj<-pslg(cbind(x,y),S=Boundaries)
#Creo la mesh
#Y dice di non aggiungere Steiner Points
#D dice di triangolare con Delaunay
mesh<-triangulate(pslg_obj,Y=TRUE,D=TRUE)
#Estrazione dei triangoli
Triang<-mesh$T
#Plot della triangolazione
plot(x,y,col="white")
for (ne in 1:dim(Triang)[1])
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
}
#Controllo se la triangolazione ha triangoli con solo punti di bordo
BorderTR<-BorderTriangles(mesh$T,Boundaries)
BorderTR
#Li coloro
for (ne in BorderTR)
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]),col="green")
}
#Ci sono comuni esterni alla frontiera?
PolyPoints<-cbind(xVenezia,yVenezia)
if(sum(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip)==length(Comuni$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}
points(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],Comuni$Latitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],col="red",pch=16)

##### COMUNI RIMASTI FUORI #####

#Quali comuni restano fuori?
Comuni$Comune[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0]

#Provo ad inserire Quartod'Altino
#Plot della triangolazione
plot(x,y,col="white",xlim=c(12.3,12.4),ylim=c(45.5,45.6))
for (ne in 1:dim(Triang)[1])
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
}
#Li coloro
for (ne in BorderTR)
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]),col="green")
}
points(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],Comuni$Latitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],col="red",pch=16)
points(xVenezia_tot,yVenezia_tot,type='l',col="red")
#Se allargo un pò la regione, ci sta dentro anche questo comune
#identify(xVenezia,yVenezia,pos=T,plot=T)
xVenezia[89]=12.38
yVenezia[89]=45.59
points(xVenezia[89],yVenezia[89])


#Riprovo tutto con queste modifiche
#Controllo di intersezioni
Intersect<-Intersections(xVenezia,yVenezia)
Intersect
#Oggetti completi
x<-c(Comuni$Longitudine,xVenezia)
y<-c(Comuni$Latitudine,yVenezia)
Codici<-Comuni$Codice
for (i in 1:length(xVenezia))
{
    Codici<-c(Codici,NA)
}
#Creo i Boundaries
Boundaries<-NULL
for(i in (length(Comuni$Longitudine)+1):(length(x)-1))
{
    Boundaries<-rbind(Boundaries, c(i,i+1))
}
Boundaries<-rbind(Boundaries, c(length(x),length(Comuni$Longitudine)+1))
#Ora triangolazione
#Oggetto pslg
pslg_obj<-pslg(cbind(x,y),S=Boundaries)
#Creo la mesh
#Y dice di non aggiungere Steiner Points
#D dice di triangolare con Delaunay
mesh<-triangulate(pslg_obj,Y=TRUE,D=TRUE)
#Estrazione dei triangoli
Triang<-mesh$T
#Plot della triangolazione
plot(x,y,col="white")
for (ne in 1:dim(Triang)[1])
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
}
#Controllo se la triangolazione ha triangoli con solo punti di bordo
BorderTR<-BorderTriangles(mesh$T,Boundaries)
BorderTR
#Li coloro
for (ne in BorderTR)
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]),col="green")
}
#Ci sono comuni esterni alla frontiera?
PolyPoints<-cbind(xVenezia,yVenezia)
points(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],Comuni$Latitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],col="red",pch=16)
dev.off()


##### SMOOTHING SPLINES DI VENEZIA #####

#Ora studio l'isola di venezia.
load("Veneto.RData")
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

#Venezia
xVenezia_tot=xbound[labels==6]
yVenezia_tot=ybound[labels==6]
s<-0
previous<-0

for(i in 2:length(xVenezia_tot))
{
    #Aggiorno la distanza
    previous<-previous+sqrt((xVenezia_tot[i]-xVenezia_tot[i-1])^2+(yVenezia_tot[i]-yVenezia_tot[i-1])^2)
    s<-c(s,previous)  
}
N<-15
#Nel range di ascissa curvilinea fisso N punti
Val<-seq(0,max(s),by=max(s)/N)
#Decido di usare delle splines cubiche
m<-3

#Troppi punti fanno in modo che si crei sempre un elevato numero di triangoli da eliminare
#Quindi devo fare qualcosa di più semmplice
#Scelgo nbasis 11
nbasis<-11
basis = create.bspline.basis(c(0,max(s)), nbasis, m)
xsmooth = smooth.basis(argvals=s, y=xVenezia_tot, basis)
ysmooth = smooth.basis(argvals=s, y=yVenezia_tot, basis)
xVenezia2  = eval.fd(Val, xsmooth$fd)
yVenezia2 = eval.fd(Val, ysmooth$fd)

#Devo modificare questa visione dell'isola
#identify(xVenezia,yVenezia,pos=T,plot=T)
IDDelete<-c(1,5:9,13,14,16)
#Alcun li devo modificare
xVenezia2[4]=12.345
yVenezia2[4]=45.443
xVenezia2[10]=12.342
yVenezia2[10]=45.433
xVenezia2[12]=12.312
yVenezia2[12]=45.434
xVenezia2=xVenezia2[-IDDelete]
yVenezia2=yVenezia2[-IDDelete]



#Unisco venezia alla regione totale
plot(xVenezia,yVenezia,type='l')
points(xVenezia2,yVenezia2,type='l')
#In pratica questi punti vanno inseriti tra il 54 e 55 punto del Venezia
#Si inseriscono i punti di venezia dal primo al settimo
#Troppo largo il settimo punto di Venezia
#Così va bene

#Ora devo unire le zone

xVenezia=c(xVenezia[1:57],12.3,xVenezia2,12.299,xVenezia[58:length(xVenezia)])
yVenezia=c(yVenezia[1:57],45.474,yVenezia2,45.4735,yVenezia[58:length(yVenezia)])

plot(xVenezia,yVenezia,type='l',xlim=c(12.26,12.36),ylim=c(45.42,45.49))
#Perfettamente unite

##### CONTROLLI FINALI DELLA REGIONE #####

x<-c(Comuni$Longitudine,xVenezia)
y<-c(Comuni$Latitudine,yVenezia)
#Creo i Boundaries
Boundaries<-NULL
for(i in (length(Comuni$Longitudine)+1):(length(x)-1))
{
    Boundaries<-rbind(Boundaries, c(i,i+1))
}
Boundaries<-rbind(Boundaries, c(length(x),length(Comuni$Longitudine)+1))
#Ora triangolazione
#Oggetto pslg
pslg_obj<-pslg(cbind(x,y),S=Boundaries)
#Creo la mesh
#Y dice di non aggiungere Steiner Points
#D dice di triangolare con Delaunay
mesh<-triangulate(pslg_obj,Y=TRUE,D=TRUE)
#Estrazione dei triangoli
Triang<-mesh$T
#Plot della triangolazione
plot(x,y,col="white",xlim=c(12.1,12.6),ylim=c(45,45.6))
for (ne in 1:dim(Triang)[1])
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
}
#Controllo se la triangolazione ha triangoli con solo punti di bordo
BorderTR<-BorderTriangles(mesh$T,Boundaries)
BorderTR
#Li coloro
for (ne in BorderTR)
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]),col="green")
}
#Ci sono comuni esterni alla frontiera?
PolyPoints<-cbind(xVenezia,yVenezia)
points(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],Comuni$Latitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],col="red",pch=16)

dev.off()


##### SALVATAGGIO DEGLI OGGETTI DEFINITIVI CON LA PENISOLA DI CT #####

#Creo i Boundaries
#Oggetti completi
x<-c(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip!=0],xVenezia)
y<-c(Comuni$Latitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip!=0],yVenezia)
Codici<-Comuni$Codice[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip!=0]
for (i in 1:length(xVenezia))
{
    Codici<-c(Codici,NA)
}
Boundaries<-NULL
for(i in (length(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip!=0])+1):(length(x)-1))
{
    Boundaries<-rbind(Boundaries, c(i,i+1))
}
Boundaries<-rbind(Boundaries, c(length(x),length(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip!=0])+1))
#Ora triangolazione
#Oggetto pslg
pslg_obj<-pslg(cbind(x,y),S=Boundaries)
#Creo la mesh
#Y dice di non aggiungere Steiner Points
#D dice di triangolare con Delaunay
mesh<-triangulate(pslg_obj,Y=TRUE,D=TRUE)
#Estrazione dei triangoli
Triang<-mesh$T


#Ora triangolazione, devo rifarla? No
plot(x,y,col="white")
for (ne in 1:dim(Triang)[1])
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
}

save(file="ProvinciaVenezia.RData",x,y,Codici,Triang)
