#RICERCA DEL MIGLIOR LAMBDA PER SR
source("2013_SSR_AllFunctions.R")
load("Territorio.RData")
library(fda)
library(rgl)
library(sp)
library(fields)
library(rgdal)

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

#Ogni anno raccolgo il massimo e il minimo della funzione stimata
Years<-1997:2011

### IMPORTANTE ###
#Ora sto facendo GCV
GCVvect<-NULL
#Lambda
logLambda<-seq(-10,3,by=0.5)
BestlogLambda<-NULL

for (year in Years)
{
    print(paste("Sto calcolando l'anno ",year,sep=""))
    
    #Scelgo i dati che mi servono...
    ANS<-NULL
    for (i in 1:length(Codici))
    {
        if(!is.na(Codici[i]))
        {
            ANS<-c(ANS,Risposta$TotalePC[Risposta$Codice==Codici[i] & Risposta$Anno==year])
        }
    }
    
    #Ora calcolo il modello
    #No edge matrix
    e = NULL
    #Ordine
    order=1
    basisobj = create.FEM.basis(cbind(x,y), e, Triang, order)
    #Creo l'oggetto fd
    fdobj = fd(numeric(basisobj$nbasis),basisobj)
    #Risposta
    Z = matrix(0,nrow=length(ANS),ncol=2)
    Z[,1] = 1:length(ANS)[1]
    Z[,2] = ANS
    #Applico il modello
    GCV = smooth.FEM.fd.GCV(Z,fdobj,logLambda)
    
    GCVvect<-c(GCVvect,GCV)

    png(filename = paste("EDF",year,".png",sep=""))
    plot(logLambda,GCV$EDFseq,main=paste("EDF",year,sep=""),xlab="logLambda",ylab="EDF",type='l',lwd=2)
    points(logLambda,GCV$EDFseq,pch=8)
    dev.off()
    png(filename = paste("Sigmahat",year,".png",sep=""))
    plot(logLambda,GCV$sigmahatseq,main=paste("Sigmahat",year,sep=""),xlab="logLambda",ylab="Sigmahat",pch=8,type='l',lwd=2)
    points(logLambda,GCV$sigmahatseq,pch=8)
    dev.off()
    png(filename = paste("GCV",year,".png",sep=""))
    plot(logLambda,GCV$GCVseq,main=paste("GCV",year,sep=""),xlab="logLambda",ylab="GCV",pch=8,type='l',lwd=2)
    points(logLambda,GCV$GCVseq,pch=8)
    dev.off()

    #Mi salvo il logLambda che realizza il minimo GCV
    BestlogLambda<-c(BestlogLambda,logLambda[which.min(GCV$GCVseq)])
}

png(filename = "Miglior logLambda.png")
plot(Years,BestlogLambda,main="BestlogLambda",xlab="Years",ylab="BestlogLambda",pch=8,type='l',lwd=2)
points(Years,BestlogLambda,pch=8)
dev.off()

save(file="GCVvect.RData",GCVvect,BestlogLambda)
