#MIGLIOR LAMBDA PER IL METODO THIN PLATE SPLINES
load("Territorio.RData")
library(fields)

##### LETTURA DEI DATI #####

Comuni<-read.table("Comuni.txt",header=T)
Risposta<-read.table("Risposta.txt",header=T)

##### CICLO PER GLI ANNI #####

#E' stato fissato negli anni l'oggetto triangolazione e contiene tutti i comuni del
#dataset

#Ogni anno raccolgo il massimo e il minimo della funzione stimata
Years<-1997:2011

#Lambda
logLambda<-seq(-10,3,by=0.5)
BestlogLambda<-NULL

for (yearindex in 1:length(Years))
{
    year<-Years[yearindex]
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
    
    #Ora voglio valutare errore di previsione
    #Non traccio su tutto il dominio il modello ma mi interesso soltanto
    #Di raccogliere l'errore.
    #Uso il leave one out
    
    #Creo i vettori in cui raccolgo gli errori devo raccogliere  gli errori
    ErrorYear<-NULL
    
    for(Lambda in logLambda)
    {
        Error<-0
        for(i in 1:length(ANS))
        {
            xK<-x[1:length(ANS)]
            yK<-y[1:length(ANS)]
            
            TPSsmooth=Tps(x=cbind(xK[-i],yK[-i]),Y=ANS[-i],lambda=10^Lambda)
            TPSz<-predict(TPSsmooth,t(c(x[i],y[i])))
            
            #Ora aggiungo i valori allo storico
            Error=Error+((TPSz-ANS[i])^2)
        }
        ErrorYear<-c(ErrorYear,Error)
    }
    
    
    #Stampo i risultati dell'anno
    png(filename = paste("Miglior logLambda",year,".png",sep=""))
    plot(logLambda,ErrorYear,main=paste("MigliorLambda",year,sep=""),pch=8,type='l',lwd=2,xlab="logLambda",ylab="Prediction Error")
    dev.off() 
    
    BestlogLambda<-c(BestlogLambda,logLambda[which.min(ErrorYear)])
}

png(filename = "Miglior logLambda.png")
plot(Years,BestlogLambda,main="BestlogLambda",xlab="Years",ylab="BestlogLambda",pch=8,type='l',lwd=2)
points(Years,BestlogLambda,pch=8)
dev.off()
