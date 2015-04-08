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

SolutionObj=ReadSolutionObjCovar(SpaceBasisObj,TimeBasisObj,"VettoreC.txt","BetaHat.txt")



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
    Mat <- matrix(ResultFitted[,j],nrow=nx,ncol=ny,byrow=F)
    Mat<-cbind(Mat,numeric(dim(Mat)[1]))
    Mat[1,ny+1]=zlim[1]
    Mat[2,ny+1]=zlim[2]
    
    # Plot
    png(filename=paste("Anno",TimePoints[j],".png",sep=""))
    image(xvec,c(yvec,50),Mat,col=heat.colors(100),main=paste("Funzione stimata anno ",TimePoints[j],sep=""),ylim=c(yvec[1],yvec[ny]),xlab="",ylab="")
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
    min(zfit)
    GoogleMap(zfit ~ lat*lon,main=paste("Funzione stimata anno ",TimePoints[j],sep=""),panel= function(...) panel.contourplot(...,labels=T,label.style="align",at=c((trunc(min(zfit))-1),seq(-700,600,by=100),(trunc(max(zfit))+1)),col.regions=rgb(1,seq(0,1,length.out=1000),0,alpha=0.5),regions=FALSE,contour=TRUE))
    
}



##### PLOT IN ALCUNE SPECIFICHE CITTA' #####

attach(CoordinateCovariate)

VeneziaCovar<-NULL
SanMicheleAlTagliamentoCovar<-NULL
LidoCovar<-NULL
PellestrinaCovar<-NULL
MuranoCovar<-NULL
CavallinoTreportiCovar<-NULL
JesoloCovar<-NULL
CaorleCovar<-NULL
ChioggiaCovar<-NULL
PortogruaroCovar<-NULL
SanDonàCovar<-NULL
CavarzereCovar<-NULL

for(j in TimePoints)
{
    #Scelgo la covariata dell'anno j
    if(j%%100<10)
    {
        name<-paste("PL0",j%%100,sep="")
    } else
    {
        name<-paste("PL",j%%100,sep="")
    }
    VeneziaCovar<-c(VeneziaCovar,(get(name)[CoordinateCovariate$Codice==425]))
    SanMicheleAlTagliamentoCovar<-c(SanMicheleAlTagliamentoCovar,(get(name)[CoordinateCovariate$Codice==417]))
    LidoCovar<-c(LidoCovar,(get(name)[CoordinateCovariate$Codice==582]))
    PellestrinaCovar<-c(PellestrinaCovar,(get(name)[CoordinateCovariate$Codice==583]))
    MuranoCovar<-c(MuranoCovar,(get(name)[CoordinateCovariate$Codice==584]))
    CavallinoTreportiCovar<-c(CavallinoTreportiCovar,(get(name)[CoordinateCovariate$Codice==581]))
    JesoloCovar<-c(JesoloCovar,(get(name)[CoordinateCovariate$Codice==402]))
    CaorleCovar<-c(CaorleCovar,(get(name)[CoordinateCovariate$Codice==388]))
    ChioggiaCovar<-c(ChioggiaCovar,(get(name)[CoordinateCovariate$Codice==391]))
    PortogruaroCovar<-c(PortogruaroCovar,(get(name)[CoordinateCovariate$Codice==412]))
    SanDonàCovar<-c(SanDonàCovar,(get(name)[CoordinateCovariate$Codice==416]))
    CavarzereCovar<-c(CavarzereCovar,(get(name)[CoordinateCovariate$Codice==389]))
    
}
detach(CoordinateCovariate)



## PLOT ##

##### GRAFICI NELLE CITTA' SCELTE #####

lim=c(400,1700)

# Plot per punto fissato
# Vediamo come va a venezia
png(filename="Venezia.png")
FixedPointPlot(12.327500,45.438056,SolutionObj,lwd=2,NameLocation="Venezia",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Venezia"]-SolutionObj$BetaHat*VeneziaCovar,col="red",pch=16)
dev.off()


# Vediamo come va a San Michele al Tagliamento
png(filename="San Michele al Tagliamento.png")
FixedPointPlot(12.994722,45.767222,SolutionObj,lwd=2,NameLocation="San Michele al Tagliamento",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="SanMichelealTagliamento"]-SolutionObj$BetaHat*SanMicheleAlTagliamentoCovar,col="red",pch=16)
dev.off()


# Vediamo come va a Bibione
png(filename="Bibione.png")
FixedPointPlot(13.062642,45.648494,SolutionObj,lwd=2,NameLocation="Bibione",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Bibione(A)"]-SolutionObj$BetaHat*SanMicheleAlTagliamentoCovar,col="red",pch=16)
dev.off()


# Vediamo come va a Lido di Venezia
png(filename="Lido(A).png")
FixedPointPlot(12.348115,45.384122,SolutionObj,lwd=2,NameLocation="Lido di Venezia",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Lido(A)"]-SolutionObj$BetaHat*LidoCovar,col="red",pch=16)
dev.off()

# Vediamo come va a Pellestrina
png(filename="Pellestrina(A).png")
FixedPointPlot(12.30181,45.27324,SolutionObj,lwd=2,NameLocation="Pellestrina",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Pellestrina(A)"]-SolutionObj$BetaHat*PellestrinaCovar,col="red",pch=16)
dev.off()

# Vediamo come va a Murano
png(filename="Murano(A).png")
FixedPointPlot(12.35155,45.45810,SolutionObj,lwd=2,NameLocation="Murano",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Murano(A)"]-SolutionObj$BetaHat*MuranoCovar,col="red",pch=16)
dev.off()

# Vediamo come va a Cavallino-Treporti
png(filename="Cavallino-Treporti.png")
FixedPointPlot(12.51000,45.46500,SolutionObj,lwd=2,NameLocation="Cavallino-Treporti",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Cavallino-Treporti"]-SolutionObj$BetaHat*CavallinoTreportiCovar,col="red",pch=16)
dev.off()

# Vediamo come va a Jesolo
png(filename="Jesolo.png")
FixedPointPlot(12.64139,45.54,SolutionObj,lwd=2,NameLocation="Jesolo",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Jesolo"]-SolutionObj$BetaHat*JesoloCovar,col="red",pch=16)
dev.off()

# Vediamo come va a Caorle
png(filename="Caorle.png")
FixedPointPlot(12.88833,45.60250,SolutionObj,lwd=2,NameLocation="Caorle",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Caorle"]-SolutionObj$BetaHat*CaorleCovar,col="red",pch=16)
dev.off()

# Vediamo come va a Chioggia
png(filename="Chioggia.png")
FixedPointPlot(12.27944,45.22056,SolutionObj,lwd=2,NameLocation="Chioggia",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Chioggia"]-SolutionObj$BetaHat*ChioggiaCovar,col="red",pch=16)
dev.off()

# Vediamo come va a Portogruaro
png(filename="Portogruaro.png")
FixedPointPlot(12.83722,45.77722,SolutionObj,lwd=2,NameLocation="Portogruaro",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Portogruaro"]-SolutionObj$BetaHat*PortogruaroCovar,col="red",pch=16)
dev.off()

# Vediamo come va a San Donà di Piave
png(filename="SanDonàdiPiave.png")
FixedPointPlot(12.56528,45.63389,SolutionObj,lwd=2,NameLocation="San Donà di Piave",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="SanDonàdiPiave"]-SolutionObj$BetaHat*SanDonàCovar,col="red",pch=16)
dev.off()

# Vediamo come va a Cavarzere
png(filename="Cavarzere.png")
FixedPointPlot(12.08389,45.13667,SolutionObj,lwd=2,NameLocation="Cavarzere",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Cavarzere"]-SolutionObj$BetaHat*CavarzereCovar,col="red",pch=16)
dev.off()




##### INTERVALLI DI CONFIDENZA #####

ICResult = ST.IC(DataMatrix,DesMat,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)
save(file="ICResult.RData",ICResult)