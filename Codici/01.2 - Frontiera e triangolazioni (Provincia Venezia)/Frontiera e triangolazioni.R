#Creazione della frontiera in poligono unico usamdo una tecnica di smoothing
#regression splines
#La provincia di Venezia risulterà così descritta:
#Ho scelto 100 punti (poi si aggiunge venezia)
#Ho ricavato che il miglior nbasis è 95


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

##### REGRESSION SPLINES #####

#Ho trovato come miglior nbasis 95
nbasis<-95
basis = create.bspline.basis(c(0,max(s)), nbasis, m)
xsmooth = smooth.basis(argvals=s, y=xVenezia_tot, basis)
ysmooth = smooth.basis(argvals=s, y=yVenezia_tot, basis)
xVenezia  = eval.fd(Val, xsmooth$fd)
yVenezia  = eval.fd(Val, ysmooth$fd)
#Traccio ora i risultati
png(filename = "Venezia - Smoothing x (regression splines).png")
plot(s,xVenezia_tot,main="Smoothing in x",xlab="Ascissa curvilinea",ylab="xVenezia")
points(Val,xVenezia,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Venezia - Smoothing y (regression splines).png")
plot(s,yVenezia_tot,main="Smoothing in y",xlab="Ascissa curvilinea",ylab="yVenezia")
points(Val,yVenezia,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Venezia - Smoothing regione (regression splines).png")
plot(xVenezia_tot,yVenezia_tot,type='l',main="Nuova definizione della regione",xlab="xVenezia",ylab="yVenezia")
points(xVenezia,yVenezia,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
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
xIsola  = eval.fd(Val, xsmooth$fd)
yIsola = eval.fd(Val, ysmooth$fd)

plot(xIsola,yIsola,type='l')
#Devo modificare questa visione dell'isola
#identify(xVenezia,yVenezia,pos=T,plot=T)
IDDelete<-c(1,5:9,13,14,16)
#Alcun li devo modificare
xIsola[4]=12.345
yIsola[4]=45.443
xIsola[10]=12.342
yIsola[10]=45.433
xIsola[12]=12.312
yIsola[12]=45.434
xIsola=xIsola[-IDDelete]
yIsola=yIsola[-IDDelete]



#Unisco venezia alla regione totale
plot(xVenezia,yVenezia,type='l',xlim=c(12.27,12.36),ylim=c(45.42,45.51))
points(xIsola,yIsola,type='l')
#identify(xVenezia,yVenezia)
#In pratica questi punti vanno inseriti tra il 54 e 55 punto del Venezia
#Si inseriscono i punti di venezia dal primo al settimo
#Troppo largo il settimo punto di Venezia
#Così va bene

#Ora devo unire le zone

xVenezia=c(xVenezia[1:57],12.3,xIsola,12.298,xVenezia[58:length(xVenezia)])
yVenezia=c(yVenezia[1:57],45.474,yIsola,45.4735,yVenezia[58:length(yVenezia)])
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
BorderTR<-BorderTriangles(mesh$T,Boundaries)


Codici<-Comuni$Codice
#Creo anche un vettore di codici di comune
for (i in (length(Codici)+1):length(x))
{
    Codici<-c(Codici,NA)
}

save(file="ConTriangoli.RData",x,y,Codici,Triang)




##### ELIMINAZIONE DEI TRIANGOLI COMPOSTI SOLO DA PUNTI DI BORDO #####

#Devo eliminare i triangoli composti solo da punti di bordo
plot(x,y,col="white",xlim=c(12.1,12.4),ylim=c(45.3,45.6))
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
#Ci sono in totale 45 triangoli che non vanno. Quali devo eliminare?

#Ce ne sono alcuni che potrei risparmiare, altri da tenere...
#Quelli da tenere
IDKeep<-NULL
#Quelli che posso togliere per una versione intermedia
IDChoose<-NULL
#Quelli assolutamente da togliere
IDDelete<-NULL

#identify(xG,yG)
#Metto la scelta su molti ma se creano troppa distanza tra comuni e bordo li tolgo
IDKeep<-c(IDKeep,8,9,10,11,12,13,17,22)
IDChoose<-c(IDChoose,1,3,6,7)
IDDelete<-c(IDDelete,2,4,5)
#Siamo a 15, 30 mancanti

#Devo eliminare i triangoli composti solo da punti di bordo
plot(x,y,col="white",xlim=c(12.4,12.7),ylim=c(45.4,45.6))
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
points(xG,yG,pch=16,col="red")
#identify(xG,yG)
IDKeep<-c(IDKeep,19,20,21,28,31,33,35,36,38,39)
IDDelete<-c(IDDelete,14,15,18,23,25,26,27,29,30,32,34,37)
IDChoose<-c(IDChoose,24,16)
#Siamo a 39, 6 mancanti

plot(x,y,col="white",xlim=c(12.9,13.1),ylim=c(45.6,45.9))
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
points(xG,yG,pch=16,col="red")
#identify(xG,yG)
IDChoose<-c(IDChoose,40,45)
IDDelete<-c(IDDelete,41,42,43,44)

#Ma questi sono secondo gli indici di BorderTR, ora devo trovare gli ID di Triang
IDDelete<-BorderTR[IDDelete]
IDChoose<-BorderTR[IDChoose]
IDkeep<-BorderTR[IDKeep]
#Sono tutti qui
#Ora devo sicuramente togliere quelli di IDtriang, e poi salvarli in un RData





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

PolyPoints<-cbind(xVenezia,yVenezia)
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

PolyPoints<-cbind(xVenezia,yVenezia)
if(sum(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip)==length(Comuni$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}
#Salvo tutto
save(file="SenzaTriangoli.RData",x,y,Codici,Triang)
