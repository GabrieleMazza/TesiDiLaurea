#PROVA CON DATASET RIDOTTO
source("Functions.R")
source("2013_SSR_AllFunctions.R")
load("Frontiera.RData")
library(fda)
library(rgl)
library(RTriangle)

##### ESTRAZIONE DEI PUNTI DI COMUNE #####

#Carico i punti di comune
Coord<-read.table("Coordinate.txt",header=T)
names(Coord)
#Leggo la risposta (rifiuti prodotto).. anche questa va ridotta
Data<-read.table("Risposta.txt",header=T)
names(Data)
#Leggo il file con le covariate
Covar<-read.table("Covariate.txt",header=T)
names(Covar)
#Leggo il file della popolazione
Pop<-read.table("Popolazione.txt",header=T)
names(Pop)

##### CICLO PER GLI ANNI #####

#Creo l'oggetto per la triangolazione
#Sebbene i comuni cambino con il corso degli anni, eseguo la triangolazione con tutti i
#581 comuni.
#Poi assegnerò i dati solo a chi deve averli

for (year in 1997:2010)
{
    print(year)
    #Prendo solo i dati che mi servono
    #DataYear<-Data$TotalePC[Data$Anno==year]
    #CodiciYear<-Data$Codice[Data$Anno==year]
    #Nel caso dovessi eliminare il comune di Ariano nel Polesine
    DataYear<-Data$TotalePC[Data$Anno==year & Data$Codice!=533]
    CodiciYear<-Data$Codice[Data$Anno==year & Data$Codice!=533]
    
    if(year%%100<10)
    {
        nameALB<-paste("ALB0",year%%100,sep="")
        nameCMP<-paste("CMP0",year%%100,sep="")
    } else
    {
        nameALB<-paste("ALB",year%%100,sep="")
        nameCMP<-paste("CMP",year%%100,sep="")
    }
    
    #Ora creo il vettore delle coordinate e delle covariate corrispondenti
    ComuniLat<-NULL
    ComuniLon<-NULL
    ALB<-NULL
    CMP<-NULL
    attach(Covar)
    for(i in 1:length(CodiciYear))
    {
        ComuniLat<-c(ComuniLat,Coord$Latitudine[Coord$Codice==CodiciYear[i]])
        ComuniLon<-c(ComuniLon,Coord$Longitudine[Coord$Codice==CodiciYear[i]])
        #Devo estrarre anche la variabile covariate
        ALB<-c(ALB,get(nameALB)[Coord$Codice==CodiciYear[i]])/Pop$Popolazione[Pop$Codice==CodiciYear[i] & Pop$Anno==year]
        CMP<-c(CMP,get(nameCMP)[Coord$Codice==CodiciYear[i]])/Pop$Popolazione[Pop$Codice==CodiciYear[i] & Pop$Anno==year]
    }
    detach(Covar)
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
    png(filename = paste("Triangolazione",year,".png",sep=""))
    plot(x,y,type="n")
    for (ne in 1:dim(Triang)[1])
    {
        polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
    }
    dev.off()
    
    #Creo la matrice disegno
    X = matrix(1,nrow=length(ALB),ncol=2)
    X[,1] = ALB
    X[,2] = CMP
    
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
    lambda = 10^(3)
    Rifiuti = smooth.FEM.fd.Covar(Z,X,Rifiutifd,lambda)
    fhat=Rifiuti$felsplobj$coef
    #Ora ho bisogno del plot della soluzione
    png(filename = paste("Risposta",year,".png",sep=""))
    plot.FEM.2D(Rifiuti$felsplobj,zlimits=c(350,500))
    dev.off()
}


