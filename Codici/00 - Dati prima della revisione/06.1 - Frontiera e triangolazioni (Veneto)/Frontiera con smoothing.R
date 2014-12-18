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

#Ora devo fare smoothing di queste due funzioni
#Come si fa?

#Quanti punti fisso per la definizione del confine?
N<-150
#Nel range di ascissa curvilinea fisso N punti
Val<-seq(0,max(s),by=max(s)/N)
#Decido di usare delle splines cubiche
m<-4


##### REGRESSION SPLINES #####

do_cycles=FALSE

##### ANALISI DEL MIGLIOR nbasis #####

#Scelgo come upper bound il numero di punti della valutazione per il numero di basi con
#cui comporre la frontiera
#Come scelgo l'ottimo del numero di basi?
#Mi interessa che sia minimizzata la distanza totale quadratica euclidea
#poligono reale del Veneto
if(do_cycles==TRUE)
{
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
}

##### RISULTATI #####

#Se io scelgo di avere 150 punti, succede che il miglior nbasis è 136 con la distanza,
#ma a livello di quanti comuni restano fuori e di quanti triangoli vengono tolti, è meglio
#nbasis 145
nbasis<-145
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
PolyPoints<-cbind(xVeneto,yVeneto)
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

#Ingrandisco nelle zone interessate dai comuni rimasti fuori

#Inizio da lestebasse
#Plot della triangolazione
plot(x,y,col="white",xlim=c(11.2,11.3),ylim=c(45.85,45.95))
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
points(xVeneto_tot,yVeneto_tot,type='l',col="red")
#Se allargo un pò la regione, ci sta dentro anche questo comune
#identify(xVeneto,yVeneto,pos=T,plot=T)
xVeneto[127]=11.26
yVeneto[127]=45.92
points(xVeneto[127],yVeneto[127])


#Ora il comune di Melara
#Plot della triangolazione
plot(x,y,col="white",xlim=c(11.10,11.30),ylim=c(44.95,45.15))
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
points(xVeneto_tot,yVeneto_tot,type='l',col="red")
#Se allargo un pò la regione, ci sta dentro anche questo comune
#identify(xVeneto,yVeneto,pos=T,plot=T)
xVeneto[105]=11.19
yVeneto[105]=45.05
points(xVeneto[105],yVeneto[105])
xVeneto[106]=11.145
yVeneto[106]=45.10
points(xVeneto[106],yVeneto[106])

#Riprovo tutto con queste modifiche
#Controllo di intersezioni
Intersect<-Intersections(xVeneto,yVeneto)
Intersect
#Oggetti completi
x<-c(Comuni$Longitudine,xVeneto)
y<-c(Comuni$Latitudine,yVeneto)
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
points(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],Comuni$Latitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],col="red",pch=16)
dev.off()

#Quali comuni restano fuori?
Comuni$Comune[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0]
#Ottimo, solo chiogggia e venezia, ma sono le isole che devo ancora aggiungere

#Ingrandisco nelle zone interessate dai comuni rimasti fuori
#Inizio da lestebasse
#Plot della triangolazione
plot(x,y,col="white",xlim=c(11.2,11.3),ylim=c(45.85,45.95))
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
points(xVeneto_tot,yVeneto_tot,type='l',col="red")
#Ottima correzione
#Ora il comune di Melara
#Plot della triangolazione
plot(x,y,col="white",xlim=c(11.10,11.30),ylim=c(44.95,45.15))
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
points(xVeneto_tot,yVeneto_tot,type='l',col="red")
#Ottima correzione


##### SMOOTHING SPLINES DI VENEZIA #####

#Ora studio l'isola di venezia.
#Anche qui devo cercare di ridurre al minimo i triangoli con vertici tutti di frontiera
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

#Scelgo come upper bound il numero di punti della valutazione per il numero di basi con
#cui comporre la frontiera
#Come scelgo l'ottimo del numero di basi?
#Mi interessa che sia minimizzata la distanza totale quadratica euclidea
#poligono reale di Venezia
if(do_cycles==TRUE)
{
bestsum<-9999999999999
best<-0
sum_vect<-NULL
best_vect<-NULL
#Qui mi salvo un po' di risultati numerici
sink(file = "ROutputVenezia.txt", append = FALSE)
print("nbasis sum best Intersections Trinagoli")
for(nbasis in (m+1):N)
{
    basis = create.bspline.basis(c(0,max(s)), nbasis, m)
    xsmooth = smooth.basis(argvals=s, y=xVenezia_tot, basis)
    ysmooth = smooth.basis(argvals=s, y=yVenezia_tot, basis)
    xVenezia  = eval.fd(Val, xsmooth$fd)
    yVenezia  = eval.fd(Val, ysmooth$fd)
    png(filename = paste(nbasis,".png",sep=""))
    plot(xVenezia_tot,yVenezia_tot,type='l',main=paste(nbasis,sep=""))
    points(xVenezia,yVenezia,col="blue",type='l')
    dev.off()
    sum=0
    for(i in 1:N)
    {
        sum<-sum+dist2Line(c(xVenezia[i],yVenezia[i]), cbind(xVenezia_tot,yVenezia_tot), distfun=SquareEuclideanDistance)[1]
    }
    if(sum<bestsum)
    {
        best<-nbasis
        bestsum=sum
    }
    Intersect<-Intersections(xVenezia,yVenezia)
    
    #Oggetti completi
    x<-c(Comuni$Longitudine[Comuni$Codice==425],xVenezia)
    y<-c(Comuni$Latitudine[Comuni$Codice==425],yVenezia)
    #Creo i Boundaries
    Boundaries<-NULL
    for(i in 2:(length(x)-1))
    {
        Boundaries<-rbind(Boundaries, c(i,i+1))
    }
    Boundaries<-rbind(Boundaries, c(length(x),2))
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
    png(filename = paste("Triangolazione a Venezia", nbasis, ".png", sep=""))
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
    points(Comuni$Longitudine[Comuni$Codice==425],Comuni$Latitudine[Comuni$Codice==425],col="red",pch=16)
    
    print(paste(
        nbasis,
        sum,
        best,
        dim(Intersect)[1],
        length(BorderTR),
        sep=" "))
    sum_vect<-c(sum_vect,sum)
    best_vect<-c(best_vect,best)
    dev.off()
}
sink()
png(filename = "Sum.png")
plot(sum_vect,type='l')
dev.off()
png(filename = "Best.png")
plot(best_vect,type='l')
dev.off()
}
#Troppi punti fanno in modo che si crei sempre un elevato numero di triangoli da eliminare
#Quindi devo fare qualcosa di più semmplice
#Scelgo nbasis 11
nbasis<-11
basis = create.bspline.basis(c(0,max(s)), nbasis, m)
xsmooth = smooth.basis(argvals=s, y=xVenezia_tot, basis)
ysmooth = smooth.basis(argvals=s, y=yVenezia_tot, basis)
xVenezia  = eval.fd(Val, xsmooth$fd)
yVenezia  = eval.fd(Val, ysmooth$fd)
#Traccio ora i risultati
png(filename = "Miglior smoothing di Venezia in x con regression splines.png")
plot(s,xVenezia_tot,main="Smoothing di xVenezia_tot")
points(Val,xVenezia,col="blue",type='l')
dev.off()
png(filename = "Miglior smoothing di Venezia in y con regression splines.png")
plot(s,yVenezia_tot,main="Smoothing di yVenezia_tot")
points(Val,yVenezia,col="blue",type='l')
dev.off()
png(filename = "Miglior smoothing di Venezia della regione con regression splines.png")
plot(xVenezia_tot,yVenezia_tot,type='l',main="Nuova definizione della regione")
points(xVenezia,yVenezia,col="blue",type='l')
dev.off()

#Controllo di intersezioni
Intersect<-Intersections(xVenezia,yVenezia)
Intersect

#Oggetti completi
x<-c(Comuni$Longitudine[Comuni$Codice==425],xVenezia)
y<-c(Comuni$Latitudine[Comuni$Codice==425],yVenezia)
#Creo i Boundaries
Boundaries<-NULL
for(i in 2:(length(x)-1))
{
    Boundaries<-rbind(Boundaries, c(i,i+1))
}
Boundaries<-rbind(Boundaries, c(length(x),2))
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
#Venezia
points(Comuni$Longitudine[Comuni$Codice==425],Comuni$Latitudine[Comuni$Codice==425],col="red",pch=16)

#Devo modificare questa visione dell'isola
#identify(xVenezia,yVenezia,pos=T,plot=T)
IDDelete<-c(1,5:9,13,14,16)
#Alcun li devo modificare
xVenezia[4]=12.345
yVenezia[4]=45.443
points(xVenezia[4],yVenezia[4])
xVenezia[10]=12.342
yVenezia[10]=45.433
points(xVenezia[10],yVenezia[10])
xVenezia[12]=12.312
yVenezia[12]=45.434
points(xVenezia[12],yVenezia[12])
xVenezia=xVenezia[-IDDelete]
yVenezia=yVenezia[-IDDelete]
#Oggetti completi
x<-c(Comuni$Longitudine[Comuni$Codice==425],xVenezia)
y<-c(Comuni$Latitudine[Comuni$Codice==425],yVenezia)
#Creo i Boundaries
Boundaries<-NULL
for(i in 2:(length(x)-1))
{
    Boundaries<-rbind(Boundaries, c(i,i+1))
}
Boundaries<-rbind(Boundaries, c(length(x),2))
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
#Venezia
points(Comuni$Longitudine[Comuni$Codice==425],Comuni$Latitudine[Comuni$Codice==425],col="red",pch=16)

#Ora è bene unire venezia al resto della regione


#Unisco venezia alla regione totale
plot(xVeneto,yVeneto,type='l',xlim=c(12.29,12.33),ylim=c(45.42,45.49))
points(xVenezia,yVenezia,type='l')
#In pratica questi punti vanno inseriti tra il 54 e 55 punto del veneto
#Si inseriscono i punti di venezia dal primo al settimo
#Troppo largo il settimo punto di Venezia
xVenezia[7]=12.3139
yVenezia[7]=45.4465
points(xVenezia[7],yVenezia[7])
points(c(xVenezia[1],12.3),c(yVenezia[1],45.474),type='l')
points(c(xVenezia[7],12.299),c(yVenezia[7],45.4735),type='l')
#Così va bene

#Ora devo unire le zone

xVeneto=c(xVeneto[1:54],12.3,xVenezia,12.299,xVeneto[55:length(xVeneto)])
yVeneto=c(yVeneto[1:54],45.474,yVenezia,45.4735,yVeneto[55:length(yVeneto)])

plot(xVeneto,yVeneto,type='l',xlim=c(12.29,12.33),ylim=c(45.42,45.49))
#Perfettamente unite

##### CHIOGGIA #####

#Ora studio l'isola di venezia.
#Anche qui devo cercare di ridurre al minimo i triangoli con vertici tutti di frontiera

plot(xVeneto,yVeneto,type='l',xlim=c(12.24,12.32),ylim=c(45.1,45.3))

xChioggia1=xbound[labels==16]
yChioggia1=ybound[labels==16]

points(xChioggia1,yChioggia1,type='l')

xChioggia2=xbound[labels==4]
yChioggia2=ybound[labels==4]

points(xChioggia2,yChioggia2,type='l')

#In pratica occorre modificare il confine tra i punt 79 e 81 perchè
#il punto 80 è finito dentro chioggia

points(12.26,45.19)
points(12.27,45.23)
points(12.31,45.235)
points(12.31,45.16)

xVeneto=c(xVeneto[1:79],12.26,12.27,12.31,12.31,xVeneto[81:length(xVeneto)])
yVeneto=c(yVeneto[1:79],45.19,45.23,45.235,45.16,yVeneto[81:length(yVeneto)])


##### CONTROLLI FINALI DELLA REGIONE #####

x<-c(Comuni$Longitudine,xVeneto)
y<-c(Comuni$Latitudine,yVeneto)
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


#Devo eliminare i triangoli composti solo da punti di bordo
plot(x,y,col="white",xlim=c(12.1,12.6),ylim=c(45,45.6))
for (ne in 1:dim(Triang)[1])
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
}
#Controllo se la triangolazione ha triangoli con solo punti di bordo
#Li coloro
for (ne in BorderTR)
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]),col="green")
}
#Controllo i baricentri dei triangoli
xG<-NULL
yG<-NULL
for (i in 1:length(BorderTR))
{
    xG<-c(xG,(x[Triang[BorderTR[i],1]]+x[Triang[BorderTR[i],2]]+x[Triang[BorderTR[i],3]])/3)
    yG<-c(yG,(y[Triang[BorderTR[i],1]]+y[Triang[BorderTR[i],2]]+y[Triang[BorderTR[i],3]])/3)
}
points(xG,yG,pch=16,col="red")
#Ci sono in totale 39 triangoli che non vanno. Quali devo eliminare?
#identify(xG,yG)
IDtriang<-c(17,18,20,21,22,23,24,7,8,1,2,3,4,5,10,13,27,28)
IDkeep<-c(6,9,11,12)
#Siamo a 22,17 mancanti

#Devo eliminare i triangoli composti solo da punti di bordo
plot(x,y,col="white",xlim=c(12.4,12.8),ylim=c(45.4,45.7))
for (ne in 1:dim(Triang)[1])
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
}
#Controllo se la triangolazione ha triangoli con solo punti di bordo
#Li coloro
BorderTR<-BorderTriangles(mesh$T,Boundaries)
BorderTR
for (ne in BorderTR)
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]),col="green")
}
#Controllo i baricentri dei triangoli
xG<-NULL
yG<-NULL
for (i in 1:length(BorderTR))
{
    xG<-c(xG,(x[Triang[BorderTR[i],1]]+x[Triang[BorderTR[i],2]]+x[Triang[BorderTR[i],3]])/3)
    yG<-c(yG,(y[Triang[BorderTR[i],1]]+y[Triang[BorderTR[i],2]]+y[Triang[BorderTR[i],3]])/3)
}
points(xG,yG,pch=16,col="red")
#Ci sono in totale 39 triangoli che non vanno. Quali devo eliminare?
#identify(xG,yG)
IDtriang<-c(17,18,20,21,22,23,24,7,8,1,2,3,4,5,10,13,27,28,32,33)
IDCT<-c(31,34,35,26,30,25,29,37)
IDkeep<-c(6,9,11,12)
points(xG[IDtriang],yG[IDtriang],pch=16,col="blue")
#Siamo ad un totale di 32 punti, ne mancano 7

#Devo eliminare i triangoli composti solo da punti di bordo
plot(x,y,col="white",xlim=c(12.3,12.5),ylim=c(44.8,45.2))
for (ne in 1:dim(Triang)[1])
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
}
#Controllo se la triangolazione ha triangoli con solo punti di bordo
#Li coloro
for (ne in BorderTR)
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]),col="green")
}
#Controllo i baricentri dei triangoli
xG<-NULL
yG<-NULL
for (i in 1:length(BorderTR))
{
    xG<-c(xG,(x[Triang[BorderTR[i],1]]+x[Triang[BorderTR[i],2]]+x[Triang[BorderTR[i],3]])/3)
    yG<-c(yG,(y[Triang[BorderTR[i],1]]+y[Triang[BorderTR[i],2]]+y[Triang[BorderTR[i],3]])/3)
}
points(xG,yG,pch=16,col="red")
#Ci sono in totale 39 triangoli che non vanno. Quali devo eliminare?
#identify(xG,yG)
IDtriang<-c(1,2,3,4,5,7,8,10,13:24,27,28,32,33,36,38,39)
IDCT<-c(25,26,29,30,31,34,35,37)
IDkeep<-c(6,9,11,12)
#Ma questi sono di BorderTR
IDtriang<-BorderTR[IDtriang]
IDCT<-BorderTR[IDCT]
IDkeep<-BorderTR[IDkeep]
#Sono tutti qui
#Ora devo sicuramente togliere quelli di IDtriang, e poi salvarli in un RData


##### SALVATAGGIO DEGLI OGGETTI DEFINITIVI CON LA PENISOLA DI CT #####

#Creo i Boundaries
x<-c(Comuni$Longitudine,xVeneto)
y<-c(Comuni$Latitudine,yVeneto)
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

L<-CleanPoints(Triang,IDtriang,x,y)
x<-L[[1]]
y<-L[[2]]
Triang<-L[[3]]
#Controllo che il nuovo bordo non abbia intersezioni
Intersect<-Intersections(x[581:length(x)],y[581:length(x)])
Intersect
#Creo anche un vettore di codici di comune
Codici<-Comuni$Codice
for (i in (length(Codici)+1):length(x))
{
    Codici<-c(Codici,NA)
}
#Ora mi salvo il nuovo contorno
Boundaries<-NULL
for(i in (length(Comuni$Longitudine)+1):(length(x)-1))
{
    Boundaries<-rbind(Boundaries, c(i,i+1))
}
Boundaries<-rbind(Boundaries, c(length(x),length(Comuni$Longitudine)+1))
#Ora triangolazione, devo rifarla? No
plot(x,y,col="white")
for (ne in 1:dim(Triang)[1])
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
}

#Controllo i punti di bordo
points(x[is.na(Codici)],y[is.na(Codici)],pch=16,col="red")
points(x[!(is.na(Codici))],y[!(is.na(Codici))],pch=16,col="blue")
BorderTR<-BorderTriangles(Triang,Boundaries)
BorderTR
#Li coloro
for (ne in BorderTR)
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]),col="green")
}
#Come mi aspettavo
#Controllo che tutti i comuni siano dentro

PolyPoints<-cbind(xVeneto,yVeneto)
if(sum(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip)==length(Comuni$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}
#Salvo tutto
save(file="ConCT.RData",x,y,Codici,Triang)


##### SALVATAGGIO DEGLI OGGETTI DEFINITIVI SENZA LA PENISOLA DI CT #####

#Creo i Boundaries
x<-c(Comuni$Longitudine,xVeneto)
y<-c(Comuni$Latitudine,yVeneto)
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

L<-CleanPoints(Triang,c(IDtriang,IDCT),x,y)
x<-L[[1]]
y<-L[[2]]
Triang<-L[[3]]
#Controllo che il nuovo bordo non abbia intersezioni
Intersect<-Intersections(x[581:length(x)],y[581:length(x)])
Intersect
#Creo anche un vettore di codici di comune
Codici<-Comuni$Codice
for (i in (length(Codici)+1):length(x))
{
    Codici<-c(Codici,NA)
}
#Ora mi salvo il nuovo contorno
Boundaries<-NULL
for(i in (length(Comuni$Longitudine)+1):(length(x)-1))
{
    Boundaries<-rbind(Boundaries, c(i,i+1))
}
Boundaries<-rbind(Boundaries, c(length(x),length(Comuni$Longitudine)+1))
#Ora triangolazione, devo rifarla? No
plot(x,y,col="white")
for (ne in 1:dim(Triang)[1])
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
}

#Controllo i punti di bordo
points(x[is.na(Codici)],y[is.na(Codici)],pch=16,col="red")
points(x[!(is.na(Codici))],y[!(is.na(Codici))],pch=16,col="blue")
BorderTR<-BorderTriangles(Triang,Boundaries)
BorderTR
#Li coloro
for (ne in BorderTR)
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]),col="green")
}
#Come mi aspettavo
#Controllo che tutti i comuni siano dentro

PolyPoints<-cbind(xVeneto,yVeneto)
if(sum(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip)==length(Comuni$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}
#Salvo tutto
save(file="SenzaCT.RData",x,y,Codici,Triang)


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