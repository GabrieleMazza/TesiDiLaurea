#ANALISI ANNO PER ANNO SR
source("Functions.R")
source("2013_SSR_AllFunctions.R")
load("Territorio.RData")
load("Maps.RData")
library(fda)
library(rgl)
library(sp)
library(RgoogleMaps)
library(loa)

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

### IMPORTANTE ###
#Devo fare tutti i plot?
PlotAll<-TRUE
Lambda = 10^-7

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
    basisobj = create.FEM.basis(cbind(x,y), e, Triang, order)
    #Creo l'oggetto fd
    fdobj = fd(numeric(basisobj$nbasis),basisobj)
    #Risposta
    Z = matrix(0,nrow=length(ANS),ncol=2)
    Z[,1] = 1:length(ANS)[1]
    Z[,2] = ANS
    #Applico il modello
    smoothobj = smooth.FEM.fd(Z,fdobj,Lambda)
    
    #Massimi e minimi
    fhat=smoothobj$felsplobj$coef
    m<-c(m,min(fhat))
    M<-c(M,max(fhat))
    
    if(PlotAll)
    {
        #Ora ho bisogno del plot della soluzione
        if (Lambda==10^-7)
        {
            limits=c(-15344,132846)
        } else
        {
            print("ERRORE! Non ho il bound del grafico")
        }
        png(filename = paste("FunzioneImage",year,".png",sep=""))
        plot.FEM.2D(smoothobj$felsplobj,zlimits=limits,title=paste("Funzione ",Years[yearindex],sep=""))
        dev.off()

        #Stampo la funzione su loa
#         xmin = min(x)
#         xmax = max(x)
#         nx   = 201
#         X    = matrix(seq(xmin, xmax, len=nx),ncol=1)
#         ymin = min(y)
#         ymax = max(y)
#         ny   = 201
#         Y    = matrix(seq(ymin, ymax, len=ny),ncol=1)    
#         
#         Xmat = X %*% matrix(1,nrow=1,ncol=ny)
#         Ymat = matrix(1,nrow=nx,ncol=1) %*% t(Y)
#         Xvec = NULL
#         for (numc in 1:nx)
#         {
#             Xvec=c(Xvec,Xmat[,numc])
#         }
#         Yvec = NULL
#         for (numc in 1:ny)
#         {
#             Yvec=c(Yvec,Ymat[,numc])
#         }
#         
#         evalmat = eval.FEM.fd(Xvec, Yvec, smoothobj$felsplobj)
#         evalmat = matrix(evalmat[,1] ,nrow=nx, ncol=ny, byrow=F)
#         #Costruisco i vettori
#         lat.loa<-NULL
#         lon.loa<-NULL
#         z.loa<-NULL
#         #Devo togliere gli na
#         for (i in 1:nx)
#         {
#             for (j in 1:ny)
#             {
#                 if(!is.na(evalmat[i,j]))
#                 {
#                     lon.loa<-c(lon.loa,X[i])
#                     lat.loa<-c(lat.loa,Y[j])
#                     z.loa<-c(z.loa,evalmat[i,j])
#                 }
#             }
#         }
#         png(filename = paste("FunzioneLoa",year,".png",sep=""))
#         GoogleMap(z.loa ~ lat.loa*lon.loa,col.regions=c("red","yellow"),contour=TRUE,alpha.regions=list(alpha=.5, alpha=.5),panel=panel.contourplot)
#         dev.off()
        
        zBar<-eval.FEM.fd(x[1:length(ANS)],y[1:length(ANS)],smoothobj$felsplobj)
        ScaledzBar<-NULL
        maximum<-max(zBar)
        minimum<-min(zBar)
        for (i in 1:length(zBar))
        {
            ScaledzBar<-c(ScaledzBar,(zBar[i]-minimum)/(maximum-minimum))
        }
        Colors = rgb(1,ScaledzBar,0)
        png(filename = paste("FunzioneMaps",year,".png",sep=""))
        PlotOnStaticMap(MapVeneto,lon=x[1:length(ANS)],lat=y[1:length(ANS)], FUN = points ,pch=16,col=Colors)
        dev.off()
        
        ScaledANS<-NULL
        maximum<-max(ANS)
        minimum<-min(ANS)
        for (i in 1:length(ANS))
        {
            ScaledANS<-c(ScaledANS,(ANS[i]-minimum)/(maximum-minimum))
        }
        Colors = rgb(1,ScaledANS,0)
        png(filename = paste("DatoInizialeMaps",year,".png",sep=""))
        PlotOnStaticMap(MapVeneto,lon=x[1:length(ANS)],lat=y[1:length(ANS)], FUN = points, pch=16,col=Colors)
        dev.off()
    }    
}

#Plot dei massimi e dei minimi
png(filename = "Massimi.png")
plot(Years,M,main="Massimi della funzione stimata per SR",xlab="Anno",ylab="Max SR",pch=8,type='l',lwd=2)
dev.off()
png(filename = "Minimi.png")
plot(Years,m,main="Minimi della funzione stimata per SR",xlab="Anno",ylab="Min SR",pch=8,type='l',lwd=2)
dev.off()

#Massimi e minimi per il grafico
print(paste("Lower Bound: ",min(m),sep=""))
print(paste("Upper Bound: ",max(M),sep=""))