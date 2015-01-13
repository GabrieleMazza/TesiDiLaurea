#STIMA DEI VARIOGRAMMI PER OGNI ANNO
load("Territorio.RData")
library(gstat)
library(geoR)
library(spatstat)

##### LETTURA DEI DATI #####

Comuni<-read.table("Comuni.txt",header=T)
Risposta<-read.table("Risposta.txt",header=T)

##### VARIOGRAMMA #####

for(year in 1997:2011)
{
    #Scelgo la risposta
    if(year%%100<10)
    {
        name<-paste("PL0",year%%100,sep="")
    } else
    {
        name<-paste("PL",year%%100,sep="")
    }
    ANS<-NULL
    attach(Comuni)
    for (i in 1:length(Codici))
    {
        if(!is.na(Codici[i]))
        {
            ANS<-c(ANS,get(name)[Comuni$Codice==Codici[i]])
        }
    }
    detach(Comuni)
    
    #Assumo solo i punti dove ho i dati
    xsp<-NULL
    ysp<-NULL
    for(i in 1:length(x))
    {
        if(!is.na(Codici[i]))
        {
            xsp<-c(xsp,x[i])
            ysp<-c(ysp,y[i])
        }
    }
    C<-data.frame(xsp,ysp)
    A<-data.frame(ANS)
    data<-SpatialPointsDataFrame(C, A)
    v <- variogram(ANS ~ 1, data)
    png(filename = paste("Stima",year,".png",sep=""))
    plot(v,pch=16,main=paste("Variogramma stimato",year,".png",sep=""))
    dev.off()
    
    #Con gstat
    v.model <- fit.variogram(v, vgm(psill=0,"Lin",range=0,nugget=1))
    
    #Provo con un modello
    png(filename = paste("Modello",year,".png",sep=""))
    plot(v, v.model, pch = 19,main=paste("Modello di variogramma",year,".png",sep=""))
    dev.off()
}