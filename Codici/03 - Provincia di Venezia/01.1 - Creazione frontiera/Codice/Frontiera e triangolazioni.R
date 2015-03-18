#Creazione della frontiera in poligono unico usamdo una tecnica di smoothing
#regression splines


#Funzioni di appoggio
source("Functions.R")
library(fda)
library(KernSmooth)
library(geosphere)
library(RTriangle)
library(SDMTools)




##### DOWNLOAD DELLA FRONTIERA #####

#Scarico la frontiera
# library(raster)
# venezia =  subset(getData('GADM', country='ITA', level=2), NAME_2=="Venezia")
# plot(venezia)
# save(file="Venezia.RData",venezia)
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

#Leggo anche i comuni (solo della provincia di venezia)
Comuni<-read.table(file="Coordinate.txt",header=T)



##### PROVINCIA DI VENEZIA #####

#Ora devo studiare il caso della provincia di Venezia
#è nel poligono 1
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

do_cycles=FALSE






##### ANALISI DEL MIGLIOR nbasis #####

#Scelgo come upper bound il numero di punti della valutazione per il numero di basi con
#cui comporre la frontiera
#Come scelgo l'ottimo del numero di basi?
#Mi interessa che sia minimizzata la distanza totale quadratica euclidea
#poligono reale della provincia di Venezia
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
        PolyPoints<-cbind(xVenezia,yVenezia)
        points(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],Comuni$Latitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],col="red",pch=16)
        dev.off()
        
        length(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0)
        comuni=length(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0])
        print(paste(
            nbasis,
            sum,
            best,
            dim(Intersect)[1],
            length(BorderTR),
            comuni,
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





##### REGRESSION SPLINES #####

# Ho trovato come miglior nbasis 95
nbasis<-95
basis = create.bspline.basis(c(0,max(s)), nbasis, m)
xsmooth = smooth.basis(argvals=s, y=xVenezia_tot, basis)
ysmooth = smooth.basis(argvals=s, y=yVenezia_tot, basis)
xVenezia  = eval.fd(Val, xsmooth$fd)
yVenezia  = eval.fd(Val, ysmooth$fd)
# Traccio ora i risultati
png(filename = "Provincia di Venezia - Smoothing x (regression splines).png")
plot(s,xVenezia_tot,main="Smoothing Longitudine",xlab="Ascissa curvilinea",ylab="Longitudine")
points(Val,xVenezia,col="red",type='l',lwd=2)
legend("bottomleft", legend=c("Reale","Smoothing"), col=c("black","red"), lty=1,lwd=2)
dev.off()
png(filename = "Provincia di Venezia - Smoothing y (regression splines).png")
plot(s,yVenezia_tot,main="Smoothing Latitudine",xlab="Ascissa curvilinea",ylab="Latitudine")
points(Val,yVenezia,col="red",type='l',lwd=2)
legend("bottomleft", legend=c("Reale","Smoothing"), col=c("black","red"), lty=1,lwd=2)
dev.off()
png(filename = "Provincia di Venezia - Smoothing regione (regression splines).png")
plot(xVenezia_tot,yVenezia_tot,type='l',main="Entroterra Poligono 1",xlab="Longitudine",ylab="Latitudine",lwd=2)
points(xVenezia,yVenezia,col="red",type='l',lwd=2)
legend("bottomright", legend=c("Reale","Smoothing"), col=c("black","red"), lty=1,lwd=2)
dev.off()


#Controllo di intersezioni
Intersect<-Intersections(xVenezia,yVenezia)
Intersect

#Ci sono comuni esterni alla frontiera?
PolyPoints<-cbind(xVenezia,yVenezia)
if(sum(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip)==length(Comuni$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}

#Quali comuni restano fuori?
Comuni$Comune[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0]




##### POLIGONO CON CHIOGGIA #####

#Ora devo studiare il caso della provincia di Venezia
#è nel poligono 1
#Studio tutti i dati in dipendenza dall'ascissa curvilinea


xChioggia_tot=xbound[labels==2]
yChioggia_tot=ybound[labels==2]
s<-0
previous<-0

for(i in 2:length(xChioggia_tot))
{
    #Aggiorno la distanza
    previous<-previous+sqrt((xChioggia_tot[i]-xChioggia_tot[i-1])^2+(yChioggia_tot[i]-yChioggia_tot[i-1])^2)
    s<-c(s,previous)  
}

#Ora devo fare smoothing di queste due funzioni
#Come si fa?

#Quanti punti fisso per la definizione del confine?
N<-50
#Nel range di ascissa curvilinea fisso N punti
Val<-seq(0,max(s),by=max(s)/N)
#Decido di usare delle splines cubiche
m<-4

do_cycles=FALSE






##### ANALISI DEL MIGLIOR nbasis #####

#Scelgo come upper bound il numero di punti della valutazione per il numero di basi con
#cui comporre la frontiera
#Come scelgo l'ottimo del numero di basi?
#Mi interessa che sia minimizzata la distanza totale quadratica euclidea
#poligono reale della provincia di Venezia, PARTE DI CHIOGGIA
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
        xsmooth = smooth.basis(argvals=s, y=xChioggia_tot, basis)
        ysmooth = smooth.basis(argvals=s, y=yChioggia_tot, basis)
        xChioggia  = eval.fd(Val, xsmooth$fd)
        yChioggia  = eval.fd(Val, ysmooth$fd)
        png(filename = paste(nbasis,".png",sep=""))
        plot(xChioggia_tot,yChioggia_tot,type='l',main=paste(nbasis,sep=""))
        points(xChioggia,yChioggia,col="blue",type='l')
        dev.off()
        sum=0
        for(i in 1:N)
        {
            sum<-sum+dist2Line(c(xChioggia[i],yChioggia[i]), cbind(xChioggia_tot,yChioggia_tot), distfun=SquareEuclideanDistance)[1]
        }
        if(sum<bestsum)
        {
            best<-nbasis
            bestsum=sum
        }
        Intersect<-Intersections(xChioggia,yChioggia)
        
        #Oggetti completi
        x<-c(Comuni$Longitudine,xChioggia)
        y<-c(Comuni$Latitudine,yChioggia)
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
        PolyPoints<-cbind(xChioggia,yChioggia)
        points(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],Comuni$Latitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],col="red",pch=16)
        dev.off()
        
        length(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0)
        comuni=length(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0])
        print(paste(
            nbasis,
            sum,
            best,
            dim(Intersect)[1],
            length(BorderTR),
            comuni,
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





##### REGRESSION SPLINES #####

# Ho trovato come miglior nbasis 25
nbasis<-25
basis = create.bspline.basis(c(0,max(s)), nbasis, m)
xsmooth = smooth.basis(argvals=s, y=xChioggia_tot, basis)
ysmooth = smooth.basis(argvals=s, y=yChioggia_tot, basis)
xChioggia  = eval.fd(Val, xsmooth$fd)
yChioggia  = eval.fd(Val, ysmooth$fd)
# Traccio ora i risultati
png(filename = "Provincia di Venezia (Chioggia) - Smoothing x (regression splines).png")
plot(s,xChioggia_tot,main="Smoothing in x",xlab="Ascissa curvilinea",ylab="xChioggia")
points(Val,xChioggia,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smoothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Provincia di Venezia (Chioggia) - Smoothing y (regression splines).png")
plot(s,yChioggia_tot,main="Smoothing in y",xlab="Ascissa curvilinea",ylab="yChioggia")
points(Val,yChioggia,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smoothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Provincia di Venezia (Chioggia) - Smoothing regione (regression splines).png")
plot(xChioggia_tot,yChioggia_tot,type='l',main="Nuova definizione della regione",xlab="xChioggia",ylab="yChioggia")
points(xChioggia,yChioggia,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smoothing"), col=c("black","blue"), lty=1)
dev.off()


#Controllo di intersezioni
Intersect<-Intersections(xChioggia,yChioggia)
Intersect



##### COMUNI RIMASTI FUORI #####

#Provo ad inserire Quartod'Altino
plot(xVenezia_tot,yVenezia_tot,col="black",xlim=c(12.3,12.4),ylim=c(45.5,45.6),type='l',main="Comune di Quarto d'Altino",xlab="xVenezia",ylab="yVenezia")
points(xVenezia,yVenezia,type='l',col="blue")
points(xVenezia,yVenezia,pch=16,col="blue")
points(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],Comuni$Latitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],col="red",pch=16)
#Se allargo un pò la regione, ci sta dentro anche questo comune
#identify(xVenezia,yVenezia,pos=T,plot=T)
xVenezia[89]=12.35
yVenezia[89]=45.59
points(xVenezia[89],yVenezia[89],col="blue",pch=16)
points(xVenezia,yVenezia,type='l',col="blue",lwd=2)



##### COLLEGAMENTO CON LA LAGUNA GIA' ANALIZZATA #####

load("PerPV.RData")

plot(xVenezia,yVenezia,type='l')
points(xVeneto,yVeneto,type='l',col="red")

#Guardo come connetterli
plot(xVenezia,yVenezia,type='l',xlim=c(12.9,13.1),ylim=c(45.4,45.8))
points(xVenezia,yVenezia)
#identify(xVenezia,yVenezia)
points(xVeneto,yVeneto,type='l',col="red")
points(xVeneto,yVeneto,col="red")
#identify(xVeneto,yVeneto)

plot(xChioggia,yChioggia,type='l',xlim=c(12.25,12.35),ylim=c(45.1,45.2))
points(xChioggia,yChioggia)
#identify(xChioggia,yChioggia)
points(xVeneto,yVeneto,type='l',col="red")
points(xVeneto,yVeneto,col="red")
#identify(xVeneto,yVeneto)


plot(xChioggia,yChioggia,type='l',xlim=c(11.9,12.35),ylim=c(45.1,45.4))
points(xChioggia,yChioggia)
#identify(xChioggia,yChioggia)
points(xVenezia,yVenezia,type='l',col="red")
points(xVenezia,yVenezia,col="red")
#identify(xVeneto,yVeneto)
points(xVeneto,yVeneto,col="blue",type="l")
points(xVenezia_tot,yVenezia_tot,type='l',col="green")

# Unisco
xVenezia<-c(xVenezia[1:9],xVeneto[25:138],xChioggia[22:43],xChioggia[46],xVenezia[74:length(xVenezia)])
yVenezia<-c(yVenezia[1:9],yVeneto[25:138],yChioggia[22:43],yChioggia[46],yVenezia[74:length(yVenezia)])

plot(xVenezia,yVenezia,type='l')

#Ci sono comuni esterni alla frontiera?
PolyPoints<-cbind(xVenezia,yVenezia)
if(sum(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip)==length(Comuni$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}




##### PLOT DEI TRIANGOLI CON PUNTI DI BORDO #####

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

# Triangoli con soli punti di Bordo
BorderTR<-BorderTriangles(mesh$T,Boundaries)
BorderTR

#Controllo i baricentri dei triangoli
xG<-NULL
yG<-NULL
for (i in 1:length(BorderTR))
{
    xG<-c(xG,(x[Triang[BorderTR[i],1]]+x[Triang[BorderTR[i],2]]+x[Triang[BorderTR[i],3]])/3)
    yG<-c(yG,(y[Triang[BorderTR[i],1]]+y[Triang[BorderTR[i],2]]+y[Triang[BorderTR[i],3]])/3)
}

for(k in 1:length(xG))
{
    xlim1<-c(xG[k]-0.02,xG[k]+0.02)
    xlim2<-c(xG[k]-0.08,xG[k]+0.08)
    
    ylim1<-c(yG[k]-0.1,yG[k]+0.1)
    ylim2<-c(yG[k]-0.3,yG[k]+0.3)
    
    
    png(filename=paste("BorderTR",k," 1.png",sep=""))
    plot(x,y,col="white",xlim=xlim1,ylim=ylim1,main=paste(k))
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
    points(xG[k],yG[k],pch=16,col="red")
    dev.off()
    
    png(filename=paste("BorderTR",k," 2.png",sep=""))
    plot(x,y,col="white",xlim=xlim2,ylim=ylim2,main=paste(k))
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
    points(xG[k],yG[k],pch=16,col="red")
    dev.off()
}

# Dall'analisi dei grafici si ha che:
IDDelete=c(3,4,5)
IDChoose=c(33,34,37,91)
IDKeep=1:94
IDKeep=IDKeep[-c(IDDelete,IDChoose)]
#Ma questi sono secondo gli indici di BorderTR, ora devo trovare gli ID di Triang
IDDelete<-BorderTR[IDDelete]
IDChoose<-BorderTR[IDChoose]
IDKeep<-BorderTR[IDKeep]



##### SALVATAGGIO DEI RISULTATI SENZA ELIMINAZIONE DEI TRIANGOLI #####

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

#Ci sono comuni esterni alla frontiera?
PolyPoints<-cbind(xVenezia,yVenezia)
if(sum(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip)==length(Comuni$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}


Codici<-Comuni$Codice
#Creo anche un vettore di codici di comune
for (i in (length(Codici)+1):length(x))
{
    Codici<-c(Codici,NA)
}

save(file="ConTriangoli.RData",x,y,Codici,Triang)









##### SALVATAGGIO DEI RISULTATI TOGLIENDO SOLO ALCUNI TRIANGOLI #####

#Creo i Boundaries
x<-c(Comuni$Longitudine,xVenezia)
y<-c(Comuni$Latitudine,yVenezia)
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

L<-CleanPoints(Triang,IDDelete,x,y)
x<-L[[1]]
y<-L[[2]]
Triang<-L[[3]]

Codici<-Comuni$Codice
#Controllo che il nuovo bordo non abbia intersezioni
Intersect<-Intersections(x[(length(Codici)+1):length(x)],y[(length(Codici)+1):length(x)])
Intersect
#Creo anche un vettore di codici di comune
for (i in (length(Codici)+1):length(x))
{
    Codici<-c(Codici,NA)
}


#Controllo che tutti i comuni siano dentro

PolyPoints<-cbind(x[is.na(Codici)],y[is.na(Codici)])
if(sum(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip)==length(Comuni$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}
#Salvo tutto
save(file="ConMenoTriangoli.RData",x,y,Codici,Triang)




##### SALVATAGGIO DEI RISULTATI TOGLIENDO TUTTI TRIANGOLI POSSIBILI #####

#Creo i Boundaries
x<-c(Comuni$Longitudine,xVenezia)
y<-c(Comuni$Latitudine,yVenezia)
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

L<-CleanPoints(Triang,c(IDDelete,IDChoose),x,y)
x<-L[[1]]
y<-L[[2]]
Triang<-L[[3]]

Codici<-Comuni$Codice
#Controllo che il nuovo bordo non abbia intersezioni
Intersect<-Intersections(x[(length(Codici)+1):length(x)],y[(length(Codici)+1):length(x)])
Intersect
#Creo anche un vettore di codici di comune
for (i in (length(Codici)+1):length(x))
{
    Codici<-c(Codici,NA)
}

#Controllo che tutti i comuni siano dentro

PolyPoints<-cbind(x[is.na(Codici)],y[is.na(Codici)])
if(sum(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip)==length(Comuni$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}
#Salvo tutto
save(file="SenzaTriangoli.RData",x,y,Codici,Triang)