#Creazione della frontiera in poligono unico usamdo una tecnica di smoothing

#Funzioni di appoggio
source("Functions.R")
library(fda)
library(KernSmooth)
library(geosphere)
library(fda)
library(RTriangle)
library(SDMTools)

##### DOWNLOAD DELLA FRONTIERA #####

#Scarico la frontiera
#require(raster)
#veneto =  subset(getData('GADM', country='ITA', level=1), NAME_1=="Veneto")
#plot(veneto)
#save(file="Veneto.RData",veneto)
load("Veneto.RData")

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

#Leggo anche i comuni
Comuni<-read.table(file="Coordinate.txt",header=T)

##### VENETO #####

#Ora devo studiare il caso del Veneto
#Il Veneto è nel poligono 1
#Studio tutti i dati in dipendenza dall'ascissa curvilinea


xVeneto_tot=xbound[labels==1]
yVeneto_tot=ybound[labels==1]
s<-0
previous<-0

for(i in 2:length(xVeneto_tot))
{
    #Aggiorno la distanza
    previous<-previous+sqrt((xVeneto_tot[i]-xVeneto_tot[i-1])^2+(yVeneto_tot[i]-yVeneto_tot[i-1])^2)
    s<-c(s,previous)  
}

plot(s,xVeneto_tot,type='l',main="xVeneto_tot")
plot(s,yVeneto_tot,type='l',main="yVeneto_tot")

#Ora devo fare smoothing di queste due funzioni
#Come si fa?

#Quanti punti fisso per la definizione del confine?
N<-30
#Nel range di ascissa curvilinea fisso N punti
Val<-seq(0,max(s),by=max(s)/N)
#Decido di usare delle splines cubiche
m<-4


##### REGRESSION SPLINES #####

#Scelgo come upper bound il numero di punti della valutazione per il numero di basi con
#cui comporre la frontiera
#Come scelgo l'ottimo del numero di basi?
#Mi interessa che sia minimizzata la distanza totale quadratica euclidea
#poligono reale del Veneto
bestsum<-9999999999999
best<-0
sum_vect<-NULL
best_vect<-NULL
#Qui mi salvo un po' di risultati numerici
sink(file = "ROutput.txt", append = FALSE)
print("nbasis sum best Intersections Trinagoli Comuni")
for(nbasis in (m+1):N)
{
    basis = create.bspline.basis(c(0,max(s)), nbasis, m)
    xsmooth = smooth.basis(argvals=s, y=xVeneto_tot, basis)
    ysmooth = smooth.basis(argvals=s, y=yVeneto_tot, basis)
    xVeneto  = eval.fd(Val, xsmooth$fd)
    yVeneto  = eval.fd(Val, ysmooth$fd)
    png(filename = paste(nbasis,".png",sep=""))
    plot(xVeneto_tot,yVeneto_tot,type='l',main=paste(nbasis,sep=""))
    points(xVeneto,yVeneto,col="blue",type='l')
    dev.off()
    sum=0
    for(i in 1:N)
    {
        sum<-sum+dist2Line(c(xVeneto[i],yVeneto[i]), cbind(xVeneto_tot,yVeneto_tot), distfun=SquareEuclideanDistance)[1]
    }
    if(sum<bestsum)
    {
        best<-nbasis
        bestsum=sum
    }
    Intersect<-Intersections(xVeneto,yVeneto)
    
    #Oggetti completi
    x<-c(Comuni$Longitudine,xVeneto)
    y<-c(Comuni$Latitudine,yVeneto)
    #Creo i Boundaries
    Boundaries<-NULL
    for(i in (length(Comuni$Longitudine)+1):length(x))
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
    png(filename = paste("Triangolazione nella laguna ", nbasis, ".png", sep=""))
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
    PolyPoints<-cbind(xVeneto,yVeneto)
    points(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],Comuni$Latitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],col="red",pch=16)
    dev.off()
    
    length(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0)
    comuniout=580-length(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0])
    print(paste(
        nbasis,
        sum,
        best,
        dim(Intersect)[1],
        length(BorderTR),
        comuniout,
        sep=" "))
    sum_vect<-c(sum_vect,sum)
    best_vect<-c(best_vect,best)
}
sink()
png(filename = "Sum.png")
plot(sum_vect,type='l')
dev.off()
png(filename = "Best.png")
plot(best_vect,type='l')
dev.off()

#Si trova che il meglio si ha con nbasis 94
nbasis<-140
basis = create.bspline.basis(c(0,max(s)), nbasis, m)
xsmooth = smooth.basis(argvals=s, y=xVeneto_tot, basis)
ysmooth = smooth.basis(argvals=s, y=yVeneto_tot, basis)
xVeneto  = eval.fd(Val, xsmooth$fd)
yVeneto  = eval.fd(Val, ysmooth$fd)
#Traccio ora i risultati
png(filename = "Miglior smoothing in x con regression splines.png")
plot(s,xVeneto_tot,main="Smoothing di xVeneto_tot")
points(Val,xVeneto,col="blue",type='l')
dev.off()
png(filename = "Miglior smoothing in y con regression splines.png")
plot(s,yVeneto_tot,main="Smoothing di yVeneto_tot")
points(Val,yVeneto,col="blue",type='l')
dev.off()
png(filename = "Miglior smoothing della laguna con regression splines.png")
plot(xVeneto_tot,yVeneto_tot,type='l',main="Nuova definizione della regione",xlim=c(12.1,12.6),ylim=c(45,45.6))
points(xVeneto,yVeneto,col="blue",type='l')
dev.off()
png(filename = "Miglior smoothing della regione con regression splines.png")
plot(xVeneto_tot,yVeneto_tot,type='l',main="Nuova definizione della regione")
points(xVeneto,yVeneto,col="blue",type='l')
dev.off()

#Controllo di intersezioni
Intersect<-Intersections(xVeneto,yVeneto)
Intersect

#Oggetti completi
x<-c(Comuni$Longitudine,xVeneto)
y<-c(Comuni$Latitudine,yVeneto)
#Creo i Boundaries
Boundaries<-NULL
for(i in (length(Comuni$Longitudine)+1):length(x))
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
png(filename = "Triangolazione nella laguna.png")
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
PolyPoints<-cbind(xVeneto,yVeneto)
if(sum(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip)==length(Comuni$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}
points(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],Comuni$Latitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],col="red",pch=16)
dev.off()

##### SMOOTHING SPLINES CON PENALIZZAZIONE DELLA DERIVATA SECONDA #####
#Se penalizzo anche con la derivata seconda, si trovano risultati poco aderenti alla
#regione iniziale. Infatti, succede che la funzione viene resa più liscia possibile,
#e quindi molte zone di contorno della regione sono appianate e molte "valli"
#riempite.
#Questo provoca una regione che si discosta molto dalla realtà, e di conseguenza un aumento
#di triangoli da eliminare e una regione non ben definita
#Ma allora è meglio considerare lambda bassi, e si hanno buoni risultati, ma sempre di più
#tendenti al caso non penalizzato con gli stessi parametri
#Ma allora tanto vale tenere il caso non penalizzato
nbasis<-140
basis = create.bspline.basis(c(0,max(s)), nbasis=nbasis, norder=m)
functionalPar = fdPar(fdobj=basis, Lfdobj=2, lambda=10^(-5))
xsmooth = smooth.basis(s, xVeneto_tot, functionalPar)
ysmooth = smooth.basis(s, yVeneto_tot, functionalPar)
xVeneto  = eval.fd(Val, xsmooth$fd)
yVeneto  = eval.fd(Val, ysmooth$fd)
#Traccio ora i risultati
png(filename = "Miglior smoothing in x con regression splines penalizzate.png")
plot(s,xVeneto_tot,main="Smoothing di xVeneto_tot")
points(Val,xVeneto,col="blue",type='l')
dev.off()
png(filename = "Miglior smoothing in y con regression splines penalizzate.png")
plot(s,yVeneto_tot,main="Smoothing di yVeneto_tot")
points(Val,yVeneto,col="blue",type='l')
dev.off()
png(filename = "Miglior smoothing della laguna con regression splines penalizzate.png")
plot(xVeneto_tot,yVeneto_tot,type='l',main="Nuova definizione della regione",xlim=c(12.1,12.6),ylim=c(45,45.6))
points(xVeneto,yVeneto,col="blue",type='l')
dev.off()
png(filename = "Miglior smoothing della regione con regression splines penalizzate.png")
plot(xVeneto_tot,yVeneto_tot,type='l',main="Nuova definizione della regione")
points(xVeneto,yVeneto,col="blue",type='l')
dev.off()

#Controllo di intersezioni
Intersect<-Intersections(xVeneto,yVeneto)
Intersect

#Oggetti completi
x<-c(Comuni$Longitudine,xVeneto)
y<-c(Comuni$Latitudine,yVeneto)
#Creo i Boundaries
Boundaries<-NULL
for(i in (length(Comuni$Longitudine)+1):length(x))
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
png(filename = "Triangolazione nella laguna nel caso penalizzato.png")
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
PolyPoints<-cbind(xVeneto,yVeneto)
if(sum(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip)==length(Comuni$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}
points(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],Comuni$Latitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],col="red",pch=16)
dev.off()


##### KERNEL SMOOTHING #####
#Il kernel smoothing è una tecnica che scarto dopo alcuni tentativi
#Infatti genera più di altri triangoli con lati di bordo
#Questo probabilmente dovuto alla tecnica ma visto che li devo eliminare
#Non considero risultati con questa tecnica
bw=0.1
m
N
xVeneto <- locpoly(s,xVeneto_tot,degree=m-1,bandwidth=bw,gridsize=N, range.x=range(s))
yVeneto <- locpoly(s,yVeneto_tot,degree=m-1,bandwidth=bw,gridsize=N, range.x=range(s))

plot(s,xVeneto_tot,main="Smoothing di xVeneto_tot")
points(xVeneto$x,xVeneto$y,col="blue",type='l')
plot(s,yVeneto_tot,main="Smoothing di yVeneto_tot")
points(xVeneto$x,yVeneto$y,col="blue",type='l')
xVeneto<-xVeneto$y
yVeneto<-yVeneto$y
plot(xVeneto_tot,yVeneto_tot,type='l',main="Nuova definizione della regione",xlim=c(12.1,12.6),ylim=c(45,45.6))
points(xVeneto,yVeneto,col="blue",type='l')

#Controllo di intersezioni
Intersect<-Intersections(xVeneto,yVeneto)
Intersect

#Oggetti completi
x<-c(Comuni$Longitudine,xVeneto)
y<-c(Comuni$Latitudine,yVeneto)
#Creo i Boundaries
Boundaries<-NULL
for(i in (length(Comuni$Longitudine)+1):length(x))
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
#png(filename = paste("Triangolazione",year,".png",sep=""))
plot(x,y,type="n")
for (ne in 1:dim(Triang)[1])
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
}
#dev.off()
#Controllo se la triangolazione ha triangoli con solo punti di bordo
BorderTR<-BorderTriangles(mesh$T,Boundaries)
BorderTR
#Li coloro
for (ne in BorderTR)
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]),col="green")
}
#Ci sono comuni esterni alla frontiera?
PolyPoints<-cbind(xVeneto,yVeneto)
if(sum(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip)==length(Comuni$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}
points(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],Comuni$Latitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],col="red",pch=16)




























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

#Ci sono punti ripetuti?
D<-Duplicated(xBoundaries,yBoundaries)
D
#No

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
