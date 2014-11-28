#PROVA ANNO PER ANNO
source("Functions.R")
source("2013_SSR_AllFunctions.R")
load("Frontiera.RData")
library(fda)
library(rgl)
library(RTriangle)
library(sp)
library(RgoogleMaps)

##### ESTRAZIONE DEI PUNTI DI COMUNE #####

#Carico i punti di comune
Coord<-read.table("Coordinate.txt",header=T)
names(Coord)
#Leggo la risposta (rifiuti prodotto).. anche questa va ridotta
Data<-read.table("Risposta.txt",header=T)
names(Data)
#Leggo il file con le covariate
PostiLetto<-read.table("PostiLetto.txt",header=T)
names(PostiLetto)

##### CICLO PER GLI ANNI #####

#Creo l'oggetto per la triangolazione
#Sebbene i comuni cambino con il corso degli anni, eseguo la triangolazione con tutti i
#581 comuni.
#Poi assegnerò i dati solo a chi deve averli

#Massimi e minimi (per i colori del plot)
m<-NULL
M<-NULL
#Ogni anno raccolgo il massimo e il minimo della funzione stimata
Years<-1997:2010
#Creo anche due matrici in cui salvo gli intervalli di confidenza, corretti e non
betahat<-NULL
CI<-matrix(0,nrow=length(Years),ncol=3)
CIsim<-matrix(0,nrow=length(Years),ncol=3)

for (yearindex in 1:length(Years))
{
    year<-Years[yearindex]
    print(paste("Sto calcolando l'anno ",year,sep=""))
    #Prendo solo i dati che mi servono
    #DataYear<-Data$TotalePC[Data$Anno==year]
    #CodiciYear<-Data$Codice[Data$Anno==year]
    #Nel caso dovessi eliminare il comune di Ariano nel Polesine
    DataYear<-Data$TotalePC[Data$Anno==year & Data$Codice!=533]
    CodiciYear<-Data$Codice[Data$Anno==year & Data$Codice!=533]
    
    if(year%%100<10)
    {
        name<-paste("PL0",year%%100,sep="")
    } else
    {
        name<-paste("PL",year%%100,sep="")
    }
    
    #Ora creo il vettore delle coordinate e della covariata corrispondenti
    ComuniLat<-NULL
    ComuniLon<-NULL
    PL<-NULL
    attach(PostiLetto)
    for(i in 1:length(CodiciYear))
    {
        ComuniLat<-c(ComuniLat,Coord$Latitudine[Coord$Codice==CodiciYear[i]])
        ComuniLon<-c(ComuniLon,Coord$Longitudine[Coord$Codice==CodiciYear[i]])
        #Devo estrarre anche la covariata
        PL<-c(PL,get(name)[Coord$Codice==CodiciYear[i]])
    }
    detach(PostiLetto)
    #Latitudine è la y, Longitudine la x
    #Ora creo la triangolazione
    #Mi servono i punti
    x<-c(ComuniLon,xBoundaries)
    y<-c(ComuniLat,yBoundaries)
    #La matrice di contorno
    Boundaries<-NULL
    for(i in (length(ComuniLon)+1):(length(x)-1))
    {
        Boundaries<-rbind(Boundaries,c(i,i+1))
    }
    Boundaries<-rbind(Boundaries,c(length(x),length(ComuniLon)+1))
    #Creo la triangolazione
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
    #plot(x,y,type="n")
    #for (ne in 1:dim(Triang)[1])
    #{
    #    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
    #}
    #dev.off()
    
    #Creo la matrice disegno
    X = matrix(1,nrow=length(PL),ncol=1)
    X[,1] = PL
    
    #Ora calcolo il modello
    #No edge matrix
    e = NULL
    #Ordine
    order=1
    basisobj = create.FEM.basis(cbind(x,y), e, Triang, order)
    #Creo l'oggetto fd
    Rifiutifd = fd(numeric(basisobj$nbasis),basisobj)
    #Risposta
    Z = matrix(0,nrow=length(DataYear),ncol=2)
    Z[,1] = 1:length(DataYear)[1]
    Z[,2] = DataYear
    #Applico il modello
    lambda = 10^0
    Rifiuti = smooth.FEM.fd.CovarCI(Z,X,Rifiutifd,lambda,c(0.95,1-0.05/14))
    
    #Ora salvo i valori che mi servono per il grafico degli intervalli di 
    #confidenza del regressore
    betahat<-c(betahat,Rifiuti$approxCIbeta[1,2])
    CI[yearindex,]=Rifiuti$approxCIbeta[1,]
    CIsim[yearindex,]<-Rifiuti$approxCIbeta[2,]
    
    #Massimi e minimi
    fhat=Rifiuti$felsplobj$coef
    m<-c(m,min(fhat))
    M<-c(M,max(fhat))
    #Ora ho bisogno del plot della soluzione
    png(filename = paste("Funzione",year,".png",sep=""))
    #Per lambda=10^2
    #plot.FEM.2D(Rifiuti$felsplobj,zlimits=c(350,427))
    #Per lambda=10^1
    #plot.FEM.2D(Rifiuti$felsplobj,zlimits=c(261,501))
    #Per lambda=10^0
    plot.FEM.2D(Rifiuti$felsplobj,zlimits=c(-175,640))
    dev.off()
    
    #Per il googlemaps plot
    #L <- list()
    #for(i in 1:dim(mesh$T)[1])
    #{
    #    Point<-NULL
    #Creo l'oggetto con i poligoni
    #    for (j in 1:3)
    #    {
    #        Point<-rbind(Point,mesh$P[mesh$T[i,j],])
    #    }
    #    Point<-rbind(Point,mesh$P[mesh$T[i,1],])
    #    str<-paste("TR",i,sep="")
    #    tmp<-Polygon(Point)
    #    Ps = Polygons(list(tmp), ID = str)
    #    L[[paste0("TR", i)]]<-Ps
    #}
    #Oggetto di tipo SpatialPolygons
    #SPs = SpatialPolygons(L)
    #Ora dovrei creare l'oggetto con la valutazione della funzione nei baricentri dei triangoli.
    #Ogni triangolo ha la funzione dipendente dal colore del suo baricentro.
    #Creo un oggetto contenente i punti dei baricentri dei triangoli
    #Barx<-NULL
    #Bary<-NULL
    #for(i in 1:dim(mesh$T)[1])
    #{
    #    xG=(mesh$P[mesh$T[i,1],1]+mesh$P[mesh$T[i,2],1]+mesh$P[mesh$T[i,3],1])/3
    #    yG=(mesh$P[mesh$T[i,1],2]+mesh$P[mesh$T[i,2],2]+mesh$P[mesh$T[i,3],2])/3
    #    Barx<-rbind(Barx,xG)
    #    Bary<-rbind(Bary,yG)
    #}
    #Ora che ho tutti i baricentri, valuto la funzione nei baricentri
    #zBar<-eval.FEM.fd(Barx,Bary,Rifiuti$felsplobj)
    #Devo normalizzare questo vettore per creare dei color
    #zCol<-(zBar-min(zBar))/(max(zBar)-min(zBar))
    #from_col <- "yellow"
    #to_col <- "red"
    #val_col <- colorRampPalette(c(from_col,to_col))(length(zCol))
    #Scarico la mappa
    #rawdata <- data.frame(as.numeric(ComuniLon), as.numeric(ComuniLat))
    #names(rawdata) <- c("lon", "lat")
    #comuni <- as.matrix(rawdata)
    #center <- rev(sapply(rawdata, mean))
    #map <- GetMap(center=center, zoom=8)
    #Plot dei poligoni
    #png(filename = paste("Map",year,".png",sep=""))
    #PlotOnStaticMap(map)
    #PlotPolysOnStaticMap(map, SPs, val_col)
    #dev.off()
}

#Plot degli intervalli di confidenza negli anni
png(filename = "IC.png")
plot(rep(Years,2),c(CIsim[,1],CIsim[,3]),type="n",main="Intervalli di confidenza di beta")
points(Years,betahat,pch=8)
for(i in 1:length(Years))
{
  plotx<-rep(Years[i],2)
  ploty1<-CIsim[i,c(1,3)]
  ploty2<-CI[i,c(1,3)]
  points(plotx,ploty1,pch=8,type='l',col="red")
  points(plotx,ploty2,pch=8,type='l',col="blue")
}
points(Years,CIsim[,1],pch=8,col="red")
points(Years,CIsim[,3],pch=8,col="red")
points(Years,CI[,1],pch=8,col="blue")
points(Years,CI[,3],pch=8,col="blue")
abline(h=0)
dev.off()
#Plot dei massimi e dei minimi
png(filename = "Massimi.png")
plot(Years,M,main="Massimi della funzione stimata",pch=8,type='l')
dev.off()
png(filename = "Minimi.png")
plot(Years,m,main="Minimi della funzione stimata",pch=8,type='l')
dev.off()
#Massimi e minimi per il grafico
print(paste("Lower Bound: ",min(m),sep=""))
print(paste("Upper Bound: ",max(M),sep=""))