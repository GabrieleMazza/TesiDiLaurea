#Creazione della frontiera in poligono unico usamdo una tecnica di smoothing
#La tecnica scelta è regression splines
#Per il Veneto sono stati scelti 150 punti
#nbasis ottimo del Veneto è 145
#Per Venezia sono stati scelti 15 punti (poi ridotti a 7)
#nbasis ottimo di Venezia è 11
#Sono stati poi modificati un paio di vertici per fare in modo che contenesse tutti i comuni
#Dall'isola di venezia sono comunque tolti triangoli puramente di bordo, perchè ce ne sono
#molti

#Funzioni di appoggio
source("Functions.R")
load("MapVeneto.RData")
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




##### REGRESSION SPLINES VENETO #####

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

##### RISULTATI E SMOOTHING FINALE #####

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
png(filename = "Veneto - Smoothing x (regression splines).png")
plot(s,xVeneto_tot,main="Smoothing in x",xlab="Ascissa curvilinea",ylab="xVeneto")
points(Val,xVeneto,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Veneto - Smoothing y (regression splines).png")
plot(s,yVeneto_tot,main="Smoothing in y",xlab="Ascissa curvilinea",ylab="yVeneto")
points(Val,yVeneto,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Veneto - Smoothing laguna (regression splines).png")
plot(xVeneto_tot,yVeneto_tot,type='l',main="Nuova definizione della regione",xlab="xVeneto",ylab="yVeneto",xlim=c(12.1,12.6),ylim=c(45,45.6))
points(xVeneto,yVeneto,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Veneto - Smoothing regione (regression splines).png")
plot(xVeneto_tot,yVeneto_tot,type='l',main="Nuova definizione della regione",xlab="xVeneto",ylab="yVeneto")
points(xVeneto,yVeneto,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()

#Controllo di intersezioni
Intersect<-Intersections(xVeneto,yVeneto)
Intersect

#Oggetti completi
x<-c(Comuni$Longitudine,xVeneto)
y<-c(Comuni$Latitudine,yVeneto)
#Controllo se ci sono dei comuni che rimangono fuori dalla frontiera
PolyPoints<-cbind(xVeneto,yVeneto)
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

#Ingrandisco nelle zone interessate dai comuni rimasti fuori

#Inizio da Lestebasse
plot(xVeneto_tot,yVeneto_tot,col="black",xlim=c(11.2,11.3),ylim=c(45.85,45.95),type='l',main="Comune di Lestebasse",xlab="xVeneto",ylab="yVeneto")
points(xVeneto,yVeneto,type='l',col="blue")
points(xVeneto,yVeneto,pch=16,col="blue")
points(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],Comuni$Latitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],col="red",pch=16)
#Se allargo un pò la regione, ci sta dentro anche questo comune
#identify(xVeneto,yVeneto,pos=T,plot=T)
xVeneto[127]=11.26
yVeneto[127]=45.92
points(xVeneto[127],yVeneto[127],col="blue",pch=16)
points(xVeneto,yVeneto,type='l',col="blue",lwd=2)

#Ora il comune di Melara
plot(xVeneto_tot,yVeneto_tot,col="black",xlim=c(11.10,11.30),ylim=c(44.95,45.15),type='l',main="Comune di Melara",xlab="xVeneto",ylab="yVeneto")
points(xVeneto,yVeneto,type='l',col="blue")
points(xVeneto,yVeneto,pch=16,col="blue")
points(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],Comuni$Latitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],col="red",pch=16)
#Se allargo un pò la regione, ci sta dentro anche questo comune
#identify(xVeneto,yVeneto,pos=T,plot=T)
xVeneto[105]=11.19
yVeneto[105]=45.05
points(xVeneto[105],yVeneto[105],col="blue",pch=16)
xVeneto[106]=11.145
yVeneto[106]=45.10
points(xVeneto[106],yVeneto[106],col="blue",pch=16)
points(xVeneto,yVeneto,type='l',col="blue",lwd=2)

#Riprovo tutto con queste modifiche
#Controllo di intersezioni
Intersect<-Intersections(xVeneto,yVeneto)
Intersect




##### SMOOTHING SPLINES VENEZIA #####

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

#Ora ho un'isola di venezia molto più piccola... Ma senza troppi triangoli di bordo. 
#Infatti non è possibile pensare che siano lasciati tutti per il fatto che è presente un solo
#dato in questa zona a fronte di così pochi punti

#Perciò in questo caso devo toglierli per forza


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

#Ora studio l'isola di chioggia. Non la descrivo troppo bene però.
#In questo caso mi basta "inglobarla" alla regione
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




##### SALVATAGGIO DEI RISULTATI SENZA ELIMINAZIONE DEI TRIANGOLI #####

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
plot(x,y,col="white",xlim=c(12.0,12.7),ylim=c(45,45.7))
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


Codici<-Comuni$Codice
#Creo anche un vettore di codici di comune
for (i in (length(Codici)+1):length(x))
{
    Codici<-c(Codici,NA)
}

save(file="ConTriangoli.RData",x,y,Codici,Triang)





##### ELIMINAZIONE DEI TRIANGOLI COMPOSTI SOLO DA PUNTI DI BORDO #####

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
#Ci sono in totale 37 triangoli che non vanno. Quali devo eliminare?

#Ce ne sono alcuni che potrei risparmiare, altri da tenere...
#Quelli da tenere
IDKeep<-NULL
#Quelli che posso togliere per una versione intermedia
IDChoose<-NULL
#Quelli assolutamente da togliere
IDDelete<-NULL

#identify(xG,yG)
IDKeep<-c(IDKeep,6,7,11,12)
IDChoose<-c(IDChoose,1,2,3,4,5,10,25,26,36)
IDDelete<-c(IDDelete,8,9,21,27,28)
#Siamo a 18, 19 mancanti

#Devo eliminare i triangoli composti solo da punti di bordo
plot(x,y,col="white",xlim=c(12.4,12.8),ylim=c(45.4,45.7))
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
IDKeep<-c(IDKeep,29,32,33)
IDDelete<-c(IDDelete,30,31)
#Siamo a 23, 14 mancanti

plot(x,y,col="white",xlim=c(12.28,12.5),ylim=c(44.7,45.2))
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
IDChoose<-c(IDChoose,20,13,24,14,15)
IDDelete<-c(IDDelete,22,23,16,17,18,19)
#Li tolgo perchè creano un'area in cui il comune è troppo distante dal bordo e tutta coperta
#Siamo a 33, 3 mancanti

plot(x,y,col="white",xlim=c(12.5,12.8),ylim=c(46.4,46.7))
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
IDChoose<-c(IDChoose,37)
#Siamo a 34, 2 mancanti


plot(x,y,col="white",xlim=c(12.8,13.1),ylim=c(45.6,45.8))
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
IDDelete<-c(IDDelete,34,35)
#Li tolgo perchè creano un'area in cui il comune è troppo distante dal bordo e tutta coperta

#Sono finiti?
Valida<-1:37
Valida[-c(IDKeep,IDDelete,IDChoose)]
#Ottimo

#Ma questi sono secondo gli indici di BorderTR, ora devo trovare gli ID di Triang
IDDelete<-BorderTR[IDDelete]
IDChoose<-BorderTR[IDChoose]
IDkeep<-BorderTR[IDKeep]
#Sono tutti qui
#Ora devo sicuramente togliere quelli di IDtriang, e poi salvarli in un RData





##### SALVATAGGIO DEI RISULTATI TOGLIENDO SOLO ALCUNI TRIANGOLI #####

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

PolyPoints<-cbind(xVeneto,yVeneto)
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

PolyPoints<-cbind(xVeneto,yVeneto)
if(sum(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip)==length(Comuni$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}
#Salvo tutto
save(file="SenzaTriangoli.RData",x,y,Codici,Triang)







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
png(filename = "Veneto - Smoothing x (regression splines penalizzate).png")
plot(s,xVeneto_tot,main="Smoothing in x",xlab="Ascissa curvilinea",ylab="xVeneto")
points(Val,xVeneto,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Veneto - Smoothing y (regression splines penalizzate).png")
plot(s,yVeneto_tot,main="Smoothing in y",xlab="Ascissa curvilinea",ylab="yVeneto")
points(Val,yVeneto,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Veneto - Smoothing laguna (regression splines penalizzate).png")
plot(xVeneto_tot,yVeneto_tot,type='l',main="Nuova definizione della regione",xlab="xVeneto",ylab="yVeneto",xlim=c(12.1,12.6),ylim=c(45,45.6))
points(xVeneto,yVeneto,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Veneto - Smoothing regione (regression splines penalizzate).png")
plot(xVeneto_tot,yVeneto_tot,type='l',main="Nuova definizione della regione",xlab="xVeneto",ylab="yVeneto")
points(xVeneto,yVeneto,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
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

#Traccio ora i risultati
png(filename = "Veneto - Smoothing x (kernel smoothing).png")
plot(s,xVeneto_tot,main="Smoothing in x",xlab="Ascissa curvilinea",ylab="xVeneto")
points(Val,xVeneto,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Veneto - Smoothing y (kernel smoothing).png")
plot(s,yVeneto_tot,main="Smoothing in y",xlab="Ascissa curvilinea",ylab="yVeneto")
points(Val,yVeneto,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Veneto - Smoothing laguna (kernel smoothing).png")
plot(xVeneto_tot,yVeneto_tot,type='l',main="Nuova definizione della regione",xlab="xVeneto",ylab="yVeneto",xlim=c(12.1,12.6),ylim=c(45,45.6))
points(xVeneto,yVeneto,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Veneto - Smoothing regione (kernel smoothing).png")
plot(xVeneto_tot,yVeneto_tot,type='l',main="Nuova definizione della regione",xlab="xVeneto",ylab="yVeneto")
points(xVeneto,yVeneto,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
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