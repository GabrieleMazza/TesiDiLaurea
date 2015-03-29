#ANALISI SPAZIO TEMPO PER LA PROVINCIA DI VENEZIA

source("SpazioTempo.R")
load("Territorio.RData")
library(rgl)
library(mgcv)
library(SDMTools)
library(loa)

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

SolutionObj=ReadSolutionObj(SpaceBasisObj,TimeBasisObj,"VettoreC.txt")



##### GRAFICI ANNO PER ANNO #####

nx<-200
ny<-200
xvec<-seq(min(x),max(x),length.out=nx)
yvec<-seq(min(y),max(y),length.out=ny)
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))

ResultFitted<-NULL

for(j in TimePoints)
{
    Time<-rep(j,length(xx))
    Result<-ST.Eval(xx,yy,Time,SolutionObj)
    
    ResultFitted<-cbind(ResultFitted,Result)
}

save(file="ResultFitted.RData",ResultFitted)

load("ResultFitted.RData")

zlim=c(min(ResultFitted,na.rm=T),max(ResultFitted,na.rm=T))

for(j in 1:length(TimePoints))
{
    # Matrice
    Mat <- matrix(ResultFitted[,j],nrow=nx,ncol=ny,byrow=F)
    Mat<-cbind(Mat,numeric(dim(Mat)[1]))
    Mat[1,ny+1]=zlim[1]
    Mat[2,ny+1]=zlim[2]
    
    # Plot
    png(filename=paste("Anno",TimePoints[j],".png",sep=""))
    image(xvec,c(yvec,50),Mat,col=heat.colors(100),main=paste("Funzione stimata anno ",TimePoints[j],sep=""),xlab="",ylab="",ylim=c(yvec[1],yvec[ny]))
    lines(xbound,ybound,lwd=1)
    lines(c(xbound[1],xbound[length(xbound)]),c(ybound[1],ybound[length(ybound)]),lwd=1)
    contour(xvec,c(yvec,50),Mat,nlevels=10,add=TRUE)
    dev.off()
}




##### GRAFICO SU GOOGLE MAPS #####

for(j in 1:length(TimePoints))
{
    
    j<-j+1
    zfit<-NULL
    lon<-NULL
    lat<-NULL
    for(i in 1:length(xx))
    {
        if(!is.na(ResultFitted[i,j]))
        {
            lon<-c(lon,xx[i])
            lat<-c(lat,yy[i])
            zfit<-c(zfit,ResultFitted[i,j])
        }
    }
    GoogleMap(zfit ~ lat*lon,main=paste("Funzione stimata anno ",TimePoints[j],sep=""),panel= function(...) panel.contourplot(...,labels=T,label.style="align",at=c(trunc(min(zfit)),seq(400,2200,by=300),(trunc(max(zfit))+1)),col.regions=rgb(1,seq(0,1,length.out=1000),0,alpha=0.5),regions=FALSE,contour=TRUE))

}


##### GRAFICI NELLE CITTA' SCELTE #####

lim=c(410,1820)

# Plot per punto fissato
# Vediamo come va a venezia
png(filename="Venezia.png")
FixedPointPlot(12.327500,45.438056,SolutionObj,lwd=2,NameLocation="Venezia",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Venezia"],col="red",pch=16)
dev.off()

# Vediamo come va a San Michele al Tagliamento
png(filename="San Michele al Tagliamento.png")
FixedPointPlot(12.994722,45.767222,SolutionObj,lwd=2,NameLocation="San Michele al Tagliamento",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="SanMichelealTagliamento"],col="red",pch=16)
dev.off()


# Vediamo come va a Bibione
png(filename="Bibione.png")
FixedPointPlot(13.062642,45.648494,SolutionObj,lwd=2,NameLocation="Bibione",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="SanMichelealTagliamento"],col="red",pch=16)
dev.off()


# Vediamo come va a Lido di Venezia
png(filename="Lido(A).png")
FixedPointPlot(12.348115,45.384122,SolutionObj,lwd=2,NameLocation="Lido di Venezia",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Lido(A)"],col="red",pch=16)
dev.off()

# Vediamo come va a Pellestrina
png(filename="Pellestrina(A).png")
FixedPointPlot(12.30181,45.27324,SolutionObj,lwd=2,NameLocation="Pellestrina",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Pellestrina(A)"],col="red",pch=16)
dev.off()

# Vediamo come va a Murano
png(filename="Murano(A).png")
FixedPointPlot(12.35155,45.45810,SolutionObj,lwd=2,NameLocation="Murano",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Murano(A)"],col="red",pch=16)
dev.off()

# Vediamo come va a Cavallino-Treporti
png(filename="Cavallino-Treporti.png")
FixedPointPlot(12.51000,45.46500,SolutionObj,lwd=2,NameLocation="Cavallino-Treporti",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Cavallino-Treporti"],col="red",pch=16)
dev.off()

# Vediamo come va a Jesolo
png(filename="Jesolo.png")
FixedPointPlot(12.64139,45.54,SolutionObj,lwd=2,NameLocation="Jesolo",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Jesolo"],col="red",pch=16)
dev.off()

# Vediamo come va a Caorle
png(filename="Caorle.png")
FixedPointPlot(12.88833,45.60250,SolutionObj,lwd=2,NameLocation="Caorle",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Caorle"],col="red",pch=16)
dev.off()

# Vediamo come va a Chioggia
png(filename="Chioggia.png")
FixedPointPlot(12.27944,45.22056,SolutionObj,lwd=2,NameLocation="Chioggia",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Chioggia"],col="red",pch=16)
dev.off()

# Vediamo come va a Portogruaro
png(filename="Portogruaro.png")
FixedPointPlot(12.83722,45.77722,SolutionObj,lwd=2,NameLocation="Portogruaro",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Portogruaro"],col="red",pch=16)
dev.off()

# Vediamo come va a San Donà di Piave
png(filename="SanDonàdiPiave.png")
FixedPointPlot(12.56528,45.63389,SolutionObj,lwd=2,NameLocation="San Donà di Piave",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="SanDonàdiPiave"],col="red",pch=16)
dev.off()

# Vediamo come va a Cavarzere
png(filename="Cavarzere.png")
FixedPointPlot(12.08389,45.13667,SolutionObj,lwd=2,NameLocation="Cavarzere",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Cavarzere"],col="red",pch=16)
dev.off()



##### PLOT RESIDUI #####

# Matrici di Dati
Data<-NULL
for(i in 1:length(Codici[1:nint]))
{
    Datatmp<-numeric(length(TimePoints))
    for(j in 1:length(TimePoints))
    {
        Datatmp[j]=Risposta$TotalePC[(Risposta$Codice==Codici[i]) & (Risposta$Anno==TimePoints[j])]
    }
    Data<-c(Data,Datatmp)
}

xknot=x[!is.na(Codici)]
yknot=y[!is.na(Codici)]
#Ricavo il vettore con tutte le stime
zHat<-NULL
for(i in 1:length(xknot))
{
    zHat<-cbind(zHat,ST.Eval(rep(xknot[i],length(TimePoints)),rep(yknot[i],length(TimePoints)),TimePoints,SolutionObj))
}

Residuals=Data-zHat

# Faccio qualche plot
png(filename="Scatterplot1.png")
plot(Data,Residuals,main="Scatterplot Residui",xlab="Dati",ylab="Residui")
abline(h=0,col="blue",lwd=2)
dev.off()
# Faccio qualche plot
png(filename="Scatterplot2.png")
plot(zHat,Residuals,main="Scatterplot Residui",xlab="Valori stimati",ylab="Residui")
abline(h=0,col="blue",lwd=2)
dev.off()

png(filename="QQplot.png")
qqnorm(Residuals)
qqline(Residuals,lwd=2,col="blue")
dev.off()

sink(file="Shapiro Test.txt")
print(shapiro.test(Residuals))
sink()
