#PROVA ANNO PER ANNO
source("Functions.R")
source("2013_SSR_AllFunctions.R")
load("Territorio.RData")
load("MapVeneto.RData")
library(fda)
library(rgl)
library(RTriangle)
library(sp)
library(RgoogleMaps)


##### LETTURA DEI DATI #####

Comuni<-read.table("Comuni.txt",header=T)
Risposta<-read.table("Risposta.txt",header=T)
#Se eventualmente servisse, plot della triangolazione
#Plot della triangolazione
#png(filename = "Triangolazione.png")
#plot(x,y,type="n")
#for (ne in 1:dim(Triang)[1])
#{
#    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
#}
#dev.off()


##### CICLO PER GLI ANNI #####

#E' stato fissato negli anni l'oggetto triangolazione e contiene tutti i comuni del
#dataset

#Massimi e minimi (per i colori del plot)
m<-NULL
M<-NULL
#Ogni anno raccolgo il massimo e il minimo della funzione stimata
Years<-1997:2011

for (yearindex in 1:length(Years))
{
    year<-Years[yearindex]
    print(paste("Sto calcolando l'anno ",year,sep=""))
    #Prendo solo i dati che mi servono
    CodiciYear<-Comuni$Codice
  
    if(year%%100<10)
    {
        name<-paste("PL0",year%%100,sep="")
    } else
    {
        name<-paste("PL",year%%100,sep="")
    }
  
    #Ho già caricato, con la triangolazione, la variabile Codici che corrisponde ai codici
    #dei comuni
    #Devo solo controllare che corrisponda
    for (i in 1:length(Codici))
    {
        if (i<=580)
        {
            if(Codici[i]!=CodiciYear[i])
            {
                print("ERRORE! COntrolla l'assegnazione dei codici di comune!")
            }
        } else
        {
            if(!is.na(Codici[i]))
            {
                print("ERRORE! COntrolla l'assegnazione dei codici di comune!")
            }
        }
    }
    
    #Quindi se tutto combacia...
    attach(Comuni)
    PL<-get(name)
    detach(Comuni)
    
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
    
    zBar<-eval.FEM.fd(x[1:580],y[1:580],PLsmooth$felsplobj)
    ScaledzBar<-NULL
    maximum<-max(zBar)
    minimum<-min(zBar)
    for (i in 1:length(zBar))
    {
      ScaledzBar<-c(ScaledzBar,(zBar[i]-minimum)/(maximum-minimum))
    }
    Colors = rgb(1,ScaledzBar,0)
    png(filename = paste("FunzioneMaps",year,".png",sep=""))
    PlotOnStaticMap(MapVeneto,lon=x[1:580],lat=y[1:580],fun="points",pch=16,col=Colors)
    dev.off()
    ScaledPL<-NULL
    maximum<-max(PL)
    minimum<-min(PL)
    for (i in 1:length(PL))
    {
      ScaledPL<-c(ScaledPL,(PL[i]-minimum)/(maximum-minimum))
    }
    Colors = rgb(1,ScaledPL,0)
    png(filename = paste("PLMaps",year,".png",sep=""))
    PlotOnStaticMap(MapVeneto,lon=x[1:580],lat=y[1:580],fun="points",pch=16,col=Colors)
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