#ANALISI SPAZIO TEMPO PER LA PROVINCIA DI VENEZIA

source("SpazioTempo.R")
load("Territorio.RData")
library(rgl)
library(mgcv)
library(SDMTools)
library(loa)
library(ggplot2)
library(ggmap)
# Leggo i dati
CoordinateCovariate<-read.table("CoordinateCovariate.txt",header=T)
Risposta<-read.table("Risposta.txt",header=T)

# Punti di tempo
TimePoints<-1997:2011
InternalPoints<-!is.na(Codici)
nint<-sum(InternalPoints)

# Basi
TimeBasisObj<-Create.Bspline.Time.Basis(TimePoints,TimeOrder=4,DerivativeOrder=2,PlotIt=F)
SpaceBasisObj<-Create.FEM.Space.Basis(cbind(x,y),Triang,InternalPoints,1)

# Matrici di dati
DataMatrix<-NULL
for(j in 1:length(TimePoints))
{
    Data<-numeric(length(Codici[1:nint]))
    for(i in 1:length(Codici[1:nint]))
    {
        Data[i]=Risposta$TotalePC[(Risposta$Codice==Codici[i]) & (Risposta$Anno==TimePoints[j])]
    }
    DataMatrix<-cbind(DataMatrix,Data)
}

#Matrici di Covariate
attach(CoordinateCovariate)
DesMat<-NULL
for(j in TimePoints)
{
    PL<-NULL
    for (i in Codici[!is.na(Codici)])
    {
        #Scelgo la covariata dell'anno j
        if(j%%100<10)
        {
            name<-paste("PL0",j%%100,sep="")
        } else
        {
            name<-paste("PL",j%%100,sep="")
        }
        PL<-c(PL,(get(name)[CoordinateCovariate$Codice==i]))
    }
    DesMat<-cbind(DesMat,PL)
}
detach(CoordinateCovariate)
DesMat=DesMat+10^-12

lon<-c(x[!is.na(Codici)],11,11)
lat<-c(y[!is.na(Codici)],0,1)
Map=get_map(location = c(lon = 12.50184, lat = 45.48403), zoom = 9, maptype = 'terrain', source="google", language = "it-IT")

Boundary=data.frame(xbound=xbound,yboun=ybound)
j=15

zfit=c(DataMatrix[,j],min(DataMatrix),max(DataMatrix))
Knot=data.frame(lon=lon,lat=lat,zfit=zfit)
ggmap(Map) + 
    geom_polygon(data=Boundary,aes(x=xbound,y=ybound),alpha=1,colour="black", fill=NA, size=0.7) +
    geom_point(data = Knot, aes(x = lon, y = lat,size=zfit),col=rgb(0/256,102/256,104/256,alpha=0.6)) + 
    scale_size_area(breaks=c(300,800,1300,1800), "Rifiuti pro capite", max_size=10) +
    labs(title=paste("Dati anno ",TimePoints[j],sep="")) +
    theme(plot.title = element_text(size = rel(1.5), face="bold")) 

zfit=c(DesMat[,j],min(DesMat),max(DesMat))
Knot=data.frame(lon=lon,lat=lat,zfit=zfit)
ggmap(Map) + 
    geom_polygon(data=Boundary,aes(x=xbound,y=ybound),alpha=1,colour="black", fill=NA, size=0.7) +
    geom_point(data = Knot, aes(x = lon, y = lat,size=zfit),col=rgb(0/256,102/256,0/256,alpha=0.6)) + 
    scale_size_area(breaks=c(1,3,5,7), "PL pro capite", max_size=10) +
    labs(title=paste("Posti letto anno ",TimePoints[j],sep="")) +
    theme(plot.title = element_text(size = rel(1.5), face="bold")) 