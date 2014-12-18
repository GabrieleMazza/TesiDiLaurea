#PROVA ANNO PER ANNO
source("Functions.R")
source("2013_SSR_AllFunctions.R")
load("Territorio.RData")
load("MapVeneto.RData")
library(fda)
library(rgl)
library(sp)
library(gstat)
library(fields)
library(rgdal)
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
#Massimi e minimi del KRIGING (per i colori del plot)
KRIGm<-NULL
KRIGM<-NULL
#Ogni anno raccolgo il massimo e il minimo della funzione stimata
Years<-1997:2011
#Creo anche due matrici in cui salvo gli intervalli di confidenza, corretti e non
betahat<-NULL
CI<-matrix(0,nrow=length(Years),ncol=3)
CIsim<-matrix(0,nrow=length(Years),ncol=3)

### IMPORTANTE ###
#E' la prima volta che stai facendo questo ciclo?
FirstCycle<-FALSE
#Devo fare anche il Kriging?
DoKriging<-FALSE
Lambda = 10^0

SSEvect<-NULL
if(DoKriging)
{
    KRIGSSEvect<-NULL
}

for (yearindex in 1:length(Years))
{
    year<-Years[yearindex]
    print(paste("Sto calcolando l'anno ",year,sep=""))
    
    #Studio la risposta
    ANS<-NULL
    for (i in 1:length(Codici))
    {
        if(!is.na(Codici[i]))
        {
            ANS<-c(ANS,Risposta$ResiduoPC[Risposta$Codice==Codici[i] & Risposta$Anno==year])
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
    #Creo la matrice disegno
    DesignMatrix = matrix(1,nrow=length(PL),ncol=1)
    DesignMatrix[,1] = PL
    
    #Applico il modello
    smoothobj = smooth.FEM.fd.CovarCI(Z,DesignMatrix,fdobj,Lambda,c(0.95,1-0.05/14))
    
    #Ora salvo i valori che mi servono per il grafico degli intervalli di 
    #confidenza del regressore
    betahat<-c(betahat,smoothobj$approxCIbeta[1,2])
    CI[yearindex,]=smoothobj$approxCIbeta[1,]
    CIsim[yearindex,]<-smoothobj$approxCIbeta[2,]
    
    #Massimi e minimi
    fhat=smoothobj$felsplobj$coef
    m<-c(m,min(fhat))
    M<-c(M,max(fhat))
    
    if(!FirstCycle)
    {
        #Ora ho bisogno del plot della soluzione
        if (Lambda==10^2)
        {
            limits=c(122,332)
        } else
        {
            if (Lambda==10^1)
            {
                limits=c(89,380)
            } else
            {
                if (Lambda==10^0)
                {
                    limits=c(-134,445)
                } else
                {
                    limits=c(max(Z[,2]),min(Z[,2]))
                }
            }
        }
        png(filename = paste("Funzione",year,".png",sep=""))
        plot.FEM.2D(smoothobj$felsplobj,zlimits=limits,title=paste("Funzione ",Years[yearindex],sep=""))
        dev.off()
        
        #Funzione senza covariate
        zBar<-eval.FEM.fd(x[1:length(ANS)],y[1:length(ANS)],smoothobj$felsplobj)
        ScaledzBar<-NULL
        maximum<-max(zBar)
        minimum<-min(zBar)
        for (i in 1:length(zBar))
        {
            ScaledzBar<-c(ScaledzBar,(zBar[i]-minimum)/(maximum-minimum))
        }
        Colors = rgb(1,ScaledzBar,0)
        png(filename = paste("FunzionePoints",year,".png",sep=""))
        PlotOnStaticMap(MapVeneto,lon=x[1:length(ANS)],lat=y[1:length(ANS)], FUN = points ,pch=16,col=Colors)
        dev.off()
        
        #Funzione con le covariate
        zBar<-zBar+PL*smoothobj$approxCIbeta[1,2]
        ScaledzBar<-NULL
        maximum<-max(zBar)
        minimum<-min(zBar)
        for (i in 1:length(zBar))
        {
            ScaledzBar<-c(ScaledzBar,(zBar[i]-minimum)/(maximum-minimum))
        }
        Colors = rgb(1,ScaledzBar,0)
        png(filename = paste("FunzioneCovariatePoints",year,".png",sep=""))
        PlotOnStaticMap(MapVeneto,lon=x[1:length(ANS)],lat=y[1:length(ANS)], FUN = points ,pch=16,col=Colors)
        dev.off()
        
        #Con questo zBar ricavo una stima dell'errore
        SSE<-0
        for (i in 1:length(ANS))
        {
            SSE=SSE+(ANS[i]-(zBar[i]+smoothobj$approxCIbeta[1,2]*PL[i]))^2
        }
        SSEvect<-c(SSEvect,SSE)
        
        ScaledANS<-NULL
        maximum<-max(ANS)
        minimum<-min(ANS)
        for (i in 1:length(ANS))
        {
            ScaledANS<-c(ScaledANS,(ANS[i]-minimum)/(maximum-minimum))
        }
        Colors = rgb(1,ScaledANS,0)
        png(filename = paste("DatoInizialePoints",year,".png",sep=""))
        PlotOnStaticMap(MapVeneto,lon=x[1:length(ANS)],lat=y[1:length(ANS)], FUN = points, pch=16,col=Colors)
        dev.off()
    }
    
    if(DoKriging)
    {
        #Ora provo a stampare su RgoogleMaps
        xmin = min(x)
        xmax = max(x)
        nx   = 201
        X    = matrix(seq(xmin, xmax, len=nx),ncol=1)
        ymin = min(y)
        ymax = max(y)
        ny   = 201
        Y    = matrix(seq(ymin, ymax, len=ny),ncol=1)    
        
        Xmat = X %*% matrix(1,nrow=1,ncol=ny)
        Ymat = matrix(1,nrow=nx,ncol=1) %*% t(Y)
        Xvec = NULL
        for (numc in 1:nx)
        {Xvec=c(Xvec,Xmat[,numc])}
        Yvec = NULL
        for (numc in 1:ny)
        {Yvec=c(Yvec,Ymat[,numc])}
        
        evalmat = eval.FEM.fd(Xvec, Yvec, smoothobj$felsplobj)
        evalmat = matrix(evalmat[,1] ,nrow=nx, ncol=ny, byrow=F)
    }
    
    if(!FirstCycle)
    {
        #png(filename = paste("FunzioneMaps",year,".png",sep=""))
        #PlotOnStaticMap(MapVeneto, add = FALSE, TrueProj=F,  FUN = image, x=X,y=Y,z=evalmat,col=heat.colors(100, alpha=0.7))
        #image(X,Y,evalmat,col=heat.colors(100, alpha=0.7),add=TRUE)
        #contour(X,Y,evalmat,add=TRUE)
        #dev.off()
    }
    
    #CONFRONTO GRAFICO CON IL KRIGING
    if(DoKriging)
    {
        #Prima creo il modello
        fit=Krig(x=cbind(x[1:length(ANS)],y[1:length(ANS)]),Y=ANS,Z=DesignMatrix,m=2)
        #Ora devo valutarlo solo dove mi interessa, non solo nel convex hull
        #Devo prendere solo i punti interni al dominio
        #Sfrutto evalmat di prima
        evalmatkrig=evalmat
        valori<-NULL
        for (i in 1:dim(evalmat)[1])
        {
            for (j in 1:dim(evalmat)[2])
            {
                if(!is.na(evalmat[i,j]))
                {
                    evalmatkrig[i,j]=predict(fit, t(c(X[i],Y[j])),drop.Z=TRUE)
                    valori<-c(valori,evalmatkrig[i,j])
                }
            }
        }
        KRIGm<-c(KRIGm,min(valori))
        KRIGM<-c(KRIGM,max(valori))
        
        Klimits<-c(-320,1125)
        
        png(filename = paste("Kriging",year,".png",sep=""))
        image(X,Y,evalmatkrig,col=heat.colors(100), xlab="", ylab="",asp=1,main=paste("Kriging ",Years[yearindex],sep=""))
        contour(X,Y,evalmatkrig,add=T)
        dev.off()
        
        #Ora ricavo una sorta di errore per il kriging
        #Con questo zBar ricavo una stima dell'errore
        SSE<-0
        for (i in 1:length(ANS))
        {
            tmp<-predict(fit, t(c(x[i],y[i])),Z=PL[i])
            SSE=SSE+(ANS[i]-tmp)^2
        }
        KRIGSSEvect<-c(KRIGSSEvect,SSE)
    }
}

#Plot degli intervalli di confidenza negli anni
png(filename = "IC.png")
plot(rep(Years,2),c(CIsim[,1],CIsim[,3]),type="n",main="Intervalli di confidenza di beta",xlab="Anno",ylab="IC")
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
legend("topright", legend=c("Corretti","Non corretti"), col=c("red","blue"), lty=1,lwd=2)
dev.off()

if(!FirstCycle)
{
    #Plot degli SSE
    png(filename = "SSE.png")
    ylim<-c(min(min(SSEvect),min(KRIGSSEvect)),max(max(SSEvect),max(KRIGSSEvect)))
    plot(Years,SSEvect,main="SSE",pch=8,type='l',lwd=2,col="red",ylim=ylim,xlab="Anno",ylab="SSE")
    points(Years,KRIGSSEvect,pch=8,type='l',lwd=2,col="blue",ylim=ylim)
    legend("topleft", legend=c("SR","KRIG"), col=c("red","blue"), lty=1,lwd=2)
    dev.off()
    dev.off()
}

#Plot dei massimi e dei minimi
png(filename = "Massimi.png")
plot(Years,M,main="Massimi della funzione stimata per SR",xlab="Anno",ylab="Max SR",pch=8,type='l',lwd=2)
dev.off()
png(filename = "Minimi.png")
plot(Years,m,main="Minimi della funzione stimata per SR",xlab="Anno",ylab="Min SR",pch=8,type='l',lwd=2)
dev.off()

if(DoKriging)
{
    #Plot dei massimi e dei minimi del kriging
    png(filename = "MassimiKriging.png")
    plot(Years,KRIGM,main="Massimi della funzione stimata per KRIG",xlab="Anno",ylab="Max KRIG",pch=8,type='l',lwd=2)
    dev.off()
    png(filename = "MinimiKriging.png")
    plot(Years,KRIGm,main="Minimi della funzione stimata per KRIG",xlab="Anno",ylab="Min KRIG",pch=8,type='l',lwd=2)
    dev.off()
}

#Massimi e minimi per il grafico
print(paste("Lower Bound: ",min(m),sep=""))
print(paste("Upper Bound: ",max(M),sep=""))

if(DoKriging)
{
    #Massimi e minimi per il grafico del kriging
    print(paste("Lower Bound Kriging: ",min(KRIGm),sep=""))
    print(paste("Upper Bound Kriging: ",max(KRIGM),sep=""))
}