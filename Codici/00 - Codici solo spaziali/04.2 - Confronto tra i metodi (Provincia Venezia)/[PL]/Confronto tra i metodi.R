#CONFRONTO DEI METODI
source("Functions.R")
source("2013_SSR_AllFunctions.R")
load("Territorio.RData")
library(fda)
library(rgl)
library(sp)
library(fields)
library(rgdal)
library(RgoogleMaps)
library(RTriangle)

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

#Ogni anno raccolgo il massimo e il minimo della funzione stimata
Years<-1997:2011
SRLambda = 10^-7
TPSLambda = 10^-3

#Qui raccolgo la somma degli errori per ogni anno
SRTot<-NULL
KRIGTot<-NULL
TPSTot<-NULL

#Voglio salvare tutto
StoricoSR<-NULL
StoricoKRIG<-NULL
StoricoTPS<-NULL

for (yearindex in 1:length(Years))
{
    year<-Years[yearindex]
    print(paste("Sto calcolando l'anno ",year,sep=""))
    
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
    
    #Ora calcolo il modello
    #No edge matrix
    e = NULL
    #Ordine
    order=1
    
    #Ora voglio valutare errore di previsione
    #Non traccio su tutto il dominio il modello ma mi interesso soltanto
    #Di raccogliere l'errore.
    #Uso il leave one out
    
    #Creo i vettori in cui raccolgo gli errori devo raccogliere  gli errori
    SRError<-NULL
    KRIGError<-NULL
    TPSError<-NULL
    
    for(i in 1:length(ANS))
    {
        
        ##############################
        ##### SPATIAL REGRESSION #####
        ##############################
        
        #Devo rifare tutto, triangolazione compresa
        xnew<-x[-i]
        ynew<-y[-i]
        Boundaries<-NULL
        #+1 -1 per il fatto che in realtà è più uno ma uso un dato in meno
        for(k in (length(ANS)+1-1):(length(xnew)-1))
        {
            Boundaries<-rbind(Boundaries, c(k,k+1))
        }
        Boundaries<-rbind(Boundaries, c(length(xnew),length(ANS)+1-1))
        #Ora triangolazione
        #Oggetto pslg
        pslg_obj<-pslg(cbind(xnew,ynew),S=Boundaries)
        #Creo la mesh
        #Y dice di non aggiungere Steiner Points
        #D dice di triangolare con Delaunay
        mesh<-triangulate(pslg_obj,Y=TRUE,D=TRUE)
        #Estrazione dei triangoli
        Triang<-mesh$T
        #Applico il modello
        #Devo rimuovere il dato i-simo
        #Risposta
        Znew = matrix(0,nrow=(length(ANS)-1),ncol=2)
        Znew[,1] = 1:(length(ANS)[1]-1)
        Znew[,2] = ANS[-i]
        #Oggetti utili
        basisobj = create.FEM.basis(cbind(xnew,ynew), e, Triang, order)
        #Creo l'oggetto fd
        fdobj = fd(numeric(basisobj$nbasis),basisobj)
        
        SRsmooth = smooth.FEM.fd(Znew,fdobj,SRLambda)
        SRz<-eval.FEM.fd(x[i],y[i],SRsmooth$felsplobj)
        
        ###################
        ##### KRIGING #####
        ###################
        
        xK<-x[1:length(ANS)]
        yK<-y[1:length(ANS)]
        KRIGsmooth=Krig(x=cbind(xK[-i],yK[-i]),Y=ANS[-i],m=2)
        KRIGz<-predict(KRIGsmooth,t(c(x[i],y[i])))
        
        ##############################
        ##### THIN PLATE SPLINES #####
        ##############################
        
        TPSsmooth=Tps(x=cbind(xK[-i],yK[-i]),Y=ANS[-i],lambda=TPSLambda)
        TPSz<-predict(TPSsmooth,t(c(x[i],y[i])))
        
        #Ora aggiungo i valori allo storico
        SRError<-c(SRError,(SRz-ANS[i])^2)
        KRIGError<-c(KRIGError,(KRIGz-ANS[i])^2)
        TPSError<-c(TPSError,(TPSz-ANS[i])^2)
        
    }
    
    #Stampo i risultati dell'anno
    png(filename = paste("Confronto",year,".png",sep=""))
    ylim<-c(min(min(SRError),min(KRIGError),min(TPSError)),max(max(SRError),max(KRIGError),max(TPSError)))
    plot(1:length(ANS),SRError,main=paste("Confronto",year,sep=""),pch=8,type='l',lwd=2,col="red",ylim=ylim,xlab="Dato",ylab="Prediction Error")
    points(1:length(ANS),KRIGError,pch=8,type='l',lwd=2,col="blue",ylim=ylim)
    points(1:length(ANS),TPSError,pch=8,type='l',lwd=2,col="darkgreen",ylim=ylim)
    legend("topleft", legend=c("SR","KRIG","TPS"), col=c("red","blue","darkgreen"), lty=1,lwd=2)
    dev.off() 
    
    #Ora aggiungo la somma
    SRTot<-c(SRTot,sum(SRError))
    KRIGTot<-c(KRIGTot,sum(KRIGError))
    TPSTot<-c(TPSTot,sum(TPSError))
    #Voglio salvare tutto
    StoricoSR<-cbind(StoricoSR,SRError)
    StoricoKRIG<-cbind(StoricoKRIG,KRIGError)
    StoricoTPS<-cbind(StoricoTPS,TPSError)
}

#Stampo la somma totale
png(filename = "Confronto tra i metodi.png")
ylim<-c(min(min(SRTot),min(KRIGTot),min(TPSTot)),max(max(SRTot),max(KRIGTot),max(TPSTot)))
plot(Years,SRTot,main="Confronto tra i metodi",pch=8,type='l',lwd=2,col="red",ylim=ylim,xlab="Dato",ylab="Prediction Error totale")
points(Years,KRIGTot,pch=8,type='l',lwd=2,col="blue",ylim=ylim)
points(Years,TPSTot,pch=8,type='l',lwd=2,col="darkgreen",ylim=ylim)
legend("topleft", legend=c("SR","KRIG","TPS"), col=c("red","blue","darkgreen"), lty=1,lwd=2)
dev.off() 

save(file="Risultati.RData",SRTot,KRIGTot,TPSTot,StoricoSR,StoricoKRIG,StoricoTPS)
