#PROVA ANNO PER ANNO
source("Functions.R")
source("2013_SSR_AllFunctions.R")
load("Frontiera.RData")
load("MapVeneto.RData")
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
    
    #Ora calcolo il modello
    #No edge matrix
    e = NULL
    #Ordine
    order=1
    basisobj = create.FEM.basis(cbind(x,y), e, Triang, order)
    #Creo l'oggetto fd
    PLfd = fd(numeric(basisobj$nbasis),basisobj)
    #Risposta
    Z = matrix(0,nrow=length(PL),ncol=2)
    Z[,1] = 1:length(PL)[1]
    Z[,2] = PL
    #Applico il modello
    lambda = 10^2
    PLsmooth = smooth.FEM.fd(Z,PLfd,lambda)
    
    #Massimi e minimi
    fhat=PLsmooth$felsplobj$coef
    m<-c(m,min(fhat))
    M<-c(M,max(fhat))
    #Ora ho bisogno del plot della soluzione
    if (lambda==10^2)
    {
      limits=c(781,2430)
    } else
    {
      if (lambda==10^1)
      {
        limits=c(379,10947)
      } else
      {
        if (lambda==10^0)
        {
          limits=c(-58,51339)
        } else
        {
          limits=c(max(Z[,2]),min(Z[,2]))
        }
      }
    }
    png(filename = paste("Risposta",year,".png",sep=""))
    plot.FEM.2D(PLsmooth$felsplobj,zlimits=limits,title=paste("Funzione ",Years[yearindex],sep=""))
    dev.off()
    
    zBar<-eval.FEM.fd(ComuniLon,ComuniLat,PLsmooth$felsplobj)
    ScaledzBar<-NULL
    maximum<-max(zBar)
    minimum<-min(zBar)
    for (i in 1:length(zBar))
    {
      ScaledzBar<-c(ScaledzBar,(zBar[i]-minimum)/(maximum-minimum))
    }
    Colors = rgb(ScaledzBar,0,0)
    png(filename = paste("FunzioneMaps",year,".png",sep=""))
    PlotOnStaticMap(MapVeneto,lon=ComuniLon,lat=ComuniLat,fun="points",pch=16,col=Colors)
    dev.off()
    ScaledPL<-NULL
    maximum<-max(PL)
    minimum<-min(PL)
    for (i in 1:length(PL))
    {
      ScaledPL<-c(ScaledPL,(PL[i]-minimum)/(maximum-minimum))
    }
    Colors = rgb(ScaledPL,0,0)
    png(filename = paste("PLMaps",year,".png",sep=""))
    PlotOnStaticMap(MapVeneto,lon=ComuniLon,lat=ComuniLat,fun="points",pch=16,col=Colors)
    dev.off()
}

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