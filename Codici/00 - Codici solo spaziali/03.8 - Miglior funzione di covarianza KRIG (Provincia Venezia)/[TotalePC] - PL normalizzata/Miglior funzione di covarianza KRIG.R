#MIGLIOR FUNZIONE DI COVARIANZA DEL COMANDO KRIG
load("Territorio.RData")
library(fields)

##### LETTURA DEI DATI #####

Comuni<-read.table("Comuni.txt",header=T)
Risposta<-read.table("Risposta.txt",header=T)


##### CICLO PER GLI ANNI #####

#Somme anno per anno
Sum1<-NULL
Sum2<-NULL
Sum3<-NULL
Sum4<-NULL
Sum5<-NULL

#Tutti gli errori anno per anno
ErrorTot1<-NULL
ErrorTot2<-NULL
ErrorTot3<-NULL
ErrorTot4<-NULL
ErrorTot5<-NULL

Years<-1997:2011

for (year in Years)
{
    print(paste("Sto calcolando l'anno ",year,sep=""))
    
    #Studio la risposta
    ANS<-NULL
    for (i in 1:length(Codici))
    {
        if(!is.na(Codici[i]))
        {
            ANS<-c(ANS,Risposta$TotalePC[Risposta$Codice==Codici[i] & Risposta$Anno==year])
        }
    }
    #Scelgo la covariata
    if(year%%100<10)
    {
        name<-paste("PL0",year%%100,sep="")
    } else
    {
        name<-paste("PL",year%%100,sep="")
    }
    PL<-NULL
    attach(Comuni)
    for (i in 1:length(Codici))
    {
        if(!is.na(Codici[i]))
        {
            PL<-c(PL,(get(name)[Comuni$Codice==Codici[i]])/1000)
        }
    }
    detach(Comuni)
    
    #Ora voglio valutare errore di previsione
    #Non traccio su tutto il dominio il modello ma mi interesso soltanto
    #Di raccogliere l'errore.
    #Uso il leave one out
    
    #Creo i vettori in cui raccolgo gli errori devo raccogliere  gli errori
    Error1<-NULL
    Error2<-NULL
    Error3<-NULL
    Error4<-NULL
    Error5<-NULL
    
    for(i in 1:length(ANS))
    {
        xK<-x[1:length(ANS)]
        yK<-y[1:length(ANS)]
        #Creo la matrice disegno
        DesignMatrixnew = matrix(1,nrow=(length(PL)-1),ncol=1)
        DesignMatrixnew[,1] = PL[-i]
        
        KRIGsmooth=Krig(x=cbind(xK[-i],yK[-i]),cov.function="Exp.cov",Y=ANS[-i],Z=DesignMatrixnew,m=2)
        KRIGz<-predict(KRIGsmooth,t(c(x[i],y[i])),Z=PL[i])
        Error1<-c(Error1,(KRIGz-ANS[i])^2)
        
        KRIGsmooth=Krig(x=cbind(xK[-i],yK[-i]),cov.function="Exp.simple.cov",Y=ANS[-i],Z=DesignMatrixnew,m=2)
        KRIGz<-predict(KRIGsmooth,t(c(x[i],y[i])),Z=PL[i])
        Error2<-c(Error2,(KRIGz-ANS[i])^2)
        
        KRIGsmooth=Krig(x=cbind(xK[-i],yK[-i]),cov.function="stationary.cov",Y=ANS[-i],Z=DesignMatrixnew,m=2)
        KRIGz<-predict(KRIGsmooth,t(c(x[i],y[i])),Z=PL[i])
        Error3<-c(Error3,(KRIGz-ANS[i])^2)
        
        KRIGsmooth=Krig(x=cbind(xK[-i],yK[-i]),cov.function="stationary.taper.cov",Y=ANS[-i],Z=DesignMatrixnew,m=2)
        KRIGz<-predict(KRIGsmooth,t(c(x[i],y[i])),Z=PL[i])
        Error4<-c(Error4,(KRIGz-ANS[i])^2)

        KRIGsmooth=Krig(x=cbind(xK[-i],yK[-i]),cov.function="wendland.cov",Y=ANS[-i],Z=DesignMatrixnew,m=2)
        KRIGz<-predict(KRIGsmooth,t(c(x[i],y[i])),Z=PL[i])
        Error5<-c(Error5,(KRIGz-ANS[i])^2)
    }

    #Stampo i risultati dell'anno
    png(filename = paste("CovarianceFunctions",year,".png",sep=""))
    m<-min(min(Error1),min(Error2),min(Error3),min(Error4),min(Error5))
    M<-max(max(Error1),max(Error2),max(Error3),max(Error4),max(Error5))
    ylim<-c(m,M)
    plot(1:length(ANS),Error1,main=paste("CovarianceFunctions",year,sep=""),pch=8,type='l',lwd=2,col="red",ylim=ylim,xlab="Dato",ylab="Prediction Error")
    points(1:length(ANS),Error2,pch=8,type='l',lwd=2,col="blue")
    points(1:length(ANS),Error3,pch=8,type='l',lwd=2,col="green")
    points(1:length(ANS),Error4,pch=8,type='l',lwd=2,col="black")
    points(1:length(ANS),Error5,pch=8,type='l',lwd=2,col="orange")
    legend("topleft", legend=c("Exp","Exp.simple","Stationary","Stationary,taper","Wendland"), col=c("red","blue","green","black","orange"), lty=1,lwd=2)
    dev.off() 
    
    #Ora aggiungo la somma
    Sum1<-c(Sum1,sum(Error1))
    Sum2<-c(Sum2,sum(Error2))
    Sum3<-c(Sum3,sum(Error3))
    Sum4<-c(Sum4,sum(Error4))
    Sum5<-c(Sum5,sum(Error5))

    #Tutti gli errori
    ErrorTot1<-cbind(ErrorTot1,Error1)
    ErrorTot2<-cbind(ErrorTot2,Error2)
    ErrorTot3<-cbind(ErrorTot3,Error3)
    ErrorTot4<-cbind(ErrorTot4,Error4)
    ErrorTot5<-cbind(ErrorTot5,Error5)
}

#Stampo la somma totale
png(filename = "CovarianceTotale.png")
m<-min(min(Sum1),min(Sum2),min(Sum3),min(Sum4),min(Sum5))
M<-max(max(Sum1),max(Sum2),max(Sum3),max(Sum4),max(Sum5))
ylim<-c(m,M)
plot(Years,Sum1,main="Prediction Error Totale",pch=8,type='l',lwd=2,col="red",ylim=ylim,xlab="Dato",ylab="Prediction Error Totale")
points(Years,Sum2,pch=8,type='l',lwd=2,col="blue")
points(Years,Sum3,pch=8,type='l',lwd=2,col="green")
points(Years,Sum4,pch=8,type='l',lwd=2,col="black")
points(Years,Sum5,pch=8,type='l',lwd=2,col="orange")
legend("topleft", legend=c("Exp","Exp.simple","Stationary","Stationary,taper","Wendland"), col=c("red","blue","green","black","orange"), lty=1,lwd=2)
dev.off() 

save(file="Risultati.RData",Sum1,Sum2,Sum3,Sum4,Sum5,ErrorTot1,ErrorTot2,ErrorTot3,ErrorTot4,ErrorTot5)