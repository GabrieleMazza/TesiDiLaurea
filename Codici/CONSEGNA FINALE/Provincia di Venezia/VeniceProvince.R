#ANALISI SPAZIO TEMPO PER LA PROVINCIA DI VENEZIA

source("ST-PDE.R")
library(rgl)
library(mgcv)
library(SDMTools)

# Data and covariate files
CoordinateCovariate<-read.table("CoordinateCovariate.txt",header=T)
Risposta<-read.table("Risposta.txt",header=T)

# Region Boundary (for plots)
Bound<-read.table("Boundary.txt",header=T)


# TIME B-SPLINES BASIS

TimePoints=1997:2011
TimeOrder=4                 # Cubic B-Splines (use order +1)
DerivativeOrder=2           # Penalizing second derivative in time
NumBasis=length(TimePoints)  # How many time basis? I use the same number as TimePoints

TimeBasisObj = create.bspline.basis(c(min(TimePoints),max(TimePoints)), NumBasis, TimeOrder)
# Use fda function!
plot(TimeBasisObj,main=paste(NumBasis,"Time B-Spline Basis"),xlab="Time",ylab="Basis")

# SPACE FE BASIS

# Data Space Points (municipalities)
SpacePoints<-CoordinateCovariate[,c(3,2)]
# Switched longitude and latitude...
# Municipalities codes (used to read data and covariate)
Codes<-CoordinateCovariate[,4]

TriangPoints=read.table(file="TriangulationPoints.txt",header=T)
Triang=read.table(file="Triangulation.txt",header=T)
FEOrder=1   #Linear Finite Elements

SpaceBasisObj<-create.FEM.basis(TriangPoints,e=NULL,Triang,FEOrder)



# DATA AND DESIGN MATRIX CREATION

# Data Vector
# If i is Space index (in SpacePoints), and j Time Index (in TimePoints)
# DATA MUST BE ORDERED FIRST VARYING TIME INDEX, THAN SPACE INDEX
# z11, z12, z13, ..., z1m, z21, z22, ..., z2m, ..., znm
# To read space points, I use municipalities Codes
Data<-NULL
for(i in Codes)
{
    for(j in TimePoints)
    {
        Data<-rbind(Data, Risposta$TotalePC[(Risposta$Codice==i) & (Risposta$Anno==j)])
    }
}
# THE SAME ORDER MUST BE USED FOR DESIGN MATRIX!
# Every row of design matrix is for (i,j) space-time point
# So Design Matrix is a nm x p matrix
# For the rows, use the same index of Data
DesMat<-NULL
for(i in 1:length(Codes))
{
    for (j in 1:length(TimePoints))
    {
        DesMat<-rbind(DesMat,CoordinateCovariate[i,4+j])
    }
}




#### MODEL WITHOUT COVARIATE ####

#Model without covariate
# GCV Optimum
LogLambdaS=-9
LogLambdaT=-0.125
SolutionObj<-ST.Smooth(Data=Data,
                       #DESMAT MUST BE A NA
                       DesMat=NA,
                       SpacePoints=SpacePoints,
                       SpaceBasisObj=SpaceBasisObj,
                       TimePoints=TimePoints,
                       TimeBasisObj=TimeBasisObj,
                       LogLambdaS=LogLambdaS,
                       LogLambdaT=LogLambdaT)
write.table(SolutionObj$cHat,file="cHat (without covariate).txt",row.names=FALSE,col.names=FALSE)

## FUNCTION PLOT AT GIVEN YEARS

nx<-200
ny<-200
xvec<-seq(min(Bound$xBound),max(Bound$xBound),length.out=nx)
yvec<-seq(min(Bound$yBound),max(Bound$yBound),length.out=ny)
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))

ResultFitted<-NULL

ntot=length(TimePoints)+1
k=1
pb <- txtProgressBar(min = 1, max = ntot, style = 3)

for(j in TimePoints)
{
    if(k==1)
    {
        setTxtProgressBar(pb, k)
    }
    Time<-rep(j,length(xx))
    Result<-ST.Eval(xx,yy,Time,SpaceBasisObj,TimeBasisObj,SolutionObj)
    
    ResultFitted<-cbind(ResultFitted,Result)
    k=k+1
    setTxtProgressBar(pb, k)
}

save(file="ResultFitted.RData",ResultFitted)
load(file="ResultFitted.RData")

zlim=c(min(ResultFitted,na.rm=T),max(ResultFitted,na.rm=T))

for(j in 1:length(TimePoints))
{
    # Matrice
    Mat <- matrix(ResultFitted[,j],nrow=nx,ncol=ny,byrow=F)
    Mat<-cbind(Mat,numeric(dim(Mat)[1]))
    Mat[1,ny+1]=zlim[1]
    Mat[2,ny+1]=zlim[2]
    
    # Plot
    png(filename=paste("Year",TimePoints[j],".png",sep=""))
    image(xvec,c(yvec,50),Mat,col=heat.colors(100),main=paste("Function in year ",TimePoints[j],sep=""),ylim=c(yvec[1],yvec[ny]))
    lines(Bound$xBound,Bound$yBound,lwd=1)
    lines(c(Bound$xBound[1],Bound$xBound[length(Bound$xBound)]),c(Bound$yBound[1],Bound$yBound[length(Bound$yBound)]),lwd=1)
    contour(xvec,c(yvec,50),Mat,nlevels=10,add=TRUE)
    dev.off()
}

## PLOT IN SOME MUNICIPALITIES

ylim=c(300,1900)
# venezia
png(filename="Venezia.png")
FixedPointPlot(12.327500,45.438056,SpaceBasisObj,TimeBasisObj,SolutionObj,lwd=2,NameLocation="Venezia",ylim=ylim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Venezia"],col="red",pch=16)
dev.off()

# San Michele al Tagliamento
png(filename="San Michele al Tagliamento.png")
FixedPointPlot(12.994722,45.767222,SpaceBasisObj,TimeBasisObj,SolutionObj,lwd=2,NameLocation="San Michele al Tagliamento",ylim=ylim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="SanMichelealTagliamento"],col="red",pch=16)
dev.off()

# Bibione
png(filename="Bibione.png")
FixedPointPlot(13.062642,45.648494,SpaceBasisObj,TimeBasisObj,SolutionObj,lwd=2,NameLocation="Bibione",ylim=ylim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="SanMichelealTagliamento"],col="red",pch=16)
dev.off()

# Lido di Venezia
png(filename="Lido(A).png")
FixedPointPlot(12.348115,45.384122,SpaceBasisObj,TimeBasisObj,SolutionObj,lwd=2,NameLocation="Lido di Venezia",ylim=ylim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Lido(A)"],col="red",pch=16)
dev.off()

# Pellestrina
png(filename="Pellestrina(A).png")
FixedPointPlot(12.30181,45.27324,SpaceBasisObj,TimeBasisObj,SolutionObj,lwd=2,NameLocation="Pellestrina",ylim=ylim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Pellestrina(A)"],col="red",pch=16)
dev.off()

# Murano
png(filename="Murano(A).png")
FixedPointPlot(12.35155,45.45810,SpaceBasisObj,TimeBasisObj,SolutionObj,lwd=2,NameLocation="Murano",ylim=ylim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Murano(A)"],col="red",pch=16)
dev.off()

# Cavallino-Treporti
png(filename="Cavallino-Treporti.png")
FixedPointPlot(12.51000,45.46500,SpaceBasisObj,TimeBasisObj,SolutionObj,lwd=2,NameLocation="Cavallino-Treporti",ylim=ylim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Cavallino-Treporti"],col="red",pch=16)
dev.off()

# Jesolo
png(filename="Jesolo.png")
FixedPointPlot(12.64139,45.54,SpaceBasisObj,TimeBasisObj,SolutionObj,lwd=2,NameLocation="Jesolo",ylim=ylim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Jesolo"],col="red",pch=16)
dev.off()

# Caorle
png(filename="Caorle.png")
FixedPointPlot(12.88833,45.60250,SpaceBasisObj,TimeBasisObj,SolutionObj,lwd=2,NameLocation="Caorle",ylim=ylim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Caorle"],col="red",pch=16)
dev.off()

# Chioggia
png(filename="Chioggia.png")
FixedPointPlot(12.27944,45.22056,SpaceBasisObj,TimeBasisObj,SolutionObj,lwd=2,NameLocation="Chioggia",ylim=ylim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Chioggia"],col="red",pch=16)
dev.off()

# Portogruaro
png(filename="Portogruaro.png")
FixedPointPlot(12.83722,45.77722,SpaceBasisObj,TimeBasisObj,SolutionObj,lwd=2,NameLocation="Portogruaro",ylim=ylim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Portogruaro"],col="red",pch=16)
dev.off()

# San Donà di Piave
png(filename="SanDonàdiPiave.png")
FixedPointPlot(12.56528,45.63389,SpaceBasisObj,TimeBasisObj,SolutionObj,lwd=2,NameLocation="San Donà di Piave",ylim=ylim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="SanDonàdiPiave"],col="red",pch=16)
dev.off()

# Cavarzere
png(filename="Cavarzere.png")
FixedPointPlot(12.08389,45.13667,SpaceBasisObj,TimeBasisObj,SolutionObj,lwd=2,NameLocation="Cavarzere",ylim=ylim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Cavarzere"],col="red",pch=16)
dev.off()













#### MODEL WITH COVARIATE ####

LogLambdaS=-9
LogLambdaT=-0.125
SolutionObjCovar<-ST.Smooth(Data=Data,
                       DesMat=DesMat,
                       SpacePoints=SpacePoints,
                       SpaceBasisObj=SpaceBasisObj,
                       TimePoints=TimePoints,
                       TimeBasisObj=TimeBasisObj,
                       LogLambdaS=LogLambdaS,
                       LogLambdaT=LogLambdaT)
write.table(SolutionObjCovar$cHat,file="cHat (with covariate).txt",row.names=FALSE,col.names=FALSE)
write.table(SolutionObjCovar$BetaHat,file="BetaHat (with covariate).txt",row.names=FALSE,col.names=FALSE)

## FUNCTION PLOT AT GIVEN YEARS

nx<-200
ny<-200
xvec<-seq(min(Bound$xBound),max(Bound$xBound),length.out=nx)
yvec<-seq(min(Bound$yBound),max(Bound$yBound),length.out=ny)
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))

ResultFitted<-NULL

ntot=length(TimePoints)+1
k=1
pb <- txtProgressBar(min = 1, max = ntot, style = 3)

for(j in TimePoints)
{
    if(k==1)
    {
        setTxtProgressBar(pb, k)
    }
    Time<-rep(j,length(xx))
    Result<-ST.Eval(xx,yy,Time,SpaceBasisObj,TimeBasisObj,SolutionObjCovar)
    
    ResultFitted<-cbind(ResultFitted,Result)
    k=k+1
    setTxtProgressBar(pb, k)
}

save(file="ResultFittedCovar.RData",ResultFitted)
load(file="ResultFittedCovar.RData")

zlim=c(min(ResultFitted,na.rm=T),max(ResultFitted,na.rm=T))

for(j in 1:length(TimePoints))
{
    # Matrice
    Mat <- matrix(ResultFitted[,j],nrow=nx,ncol=ny,byrow=F)
    Mat<-cbind(Mat,numeric(dim(Mat)[1]))
    Mat[1,ny+1]=zlim[1]
    Mat[2,ny+1]=zlim[2]
    
    # Plot
    png(filename=paste("Year",TimePoints[j],".png",sep=""))
    image(xvec,c(yvec,50),Mat,col=heat.colors(100),main=paste("Function in year ",TimePoints[j],sep=""),ylim=c(yvec[1],yvec[ny]))
    lines(Bound$xBound,Bound$yBound,lwd=1)
    lines(c(Bound$xBound[1],Bound$xBound[length(Bound$xBound)]),c(Bound$yBound[1],Bound$yBound[length(Bound$yBound)]),lwd=1)
    contour(xvec,c(yvec,50),Mat,nlevels=10,add=TRUE)
    dev.off()
}

## PLOT IN SOME MUNICIPALITIES

#First of all, I need covariates
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

lim=c(400,1700)

# Plot per punto fissato
# venezia
png(filename="Venezia.png")
FixedPointPlot(12.327500,45.438056,SpaceBasisObj,TimeBasisObj,SolutionObjCovar,lwd=2,NameLocation="Venezia",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Venezia"]-SolutionObjCovar$BetaHat*VeneziaCovar,col="red",pch=16)
dev.off()

# San Michele al Tagliamento
png(filename="San Michele al Tagliamento.png")
FixedPointPlot(12.994722,45.767222,SpaceBasisObj,TimeBasisObj,SolutionObjCovar,lwd=2,NameLocation="San Michele al Tagliamento",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="SanMichelealTagliamento"]-SolutionObjCovar$BetaHat*SanMicheleAlTagliamentoCovar,col="red",pch=16)
dev.off()

# Bibione
png(filename="Bibione.png")
FixedPointPlot(13.062642,45.648494,SpaceBasisObj,TimeBasisObj,SolutionObjCovar,lwd=2,NameLocation="Bibione",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Bibione(A)"]-SolutionObjCovar$BetaHat*SanMicheleAlTagliamentoCovar,col="red",pch=16)
dev.off()

# Lido di Venezia
png(filename="Lido(A).png")
FixedPointPlot(12.348115,45.384122,SpaceBasisObj,TimeBasisObj,SolutionObjCovar,lwd=2,NameLocation="Lido di Venezia",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Lido(A)"]-SolutionObjCovar$BetaHat*LidoCovar,col="red",pch=16)
dev.off()

# Pellestrina
png(filename="Pellestrina(A).png")
FixedPointPlot(12.30181,45.27324,SpaceBasisObj,TimeBasisObj,SolutionObjCovar,lwd=2,NameLocation="Pellestrina",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Pellestrina(A)"]-SolutionObjCovar$BetaHat*PellestrinaCovar,col="red",pch=16)
dev.off()

# Murano
png(filename="Murano(A).png")
FixedPointPlot(12.35155,45.45810,SpaceBasisObj,TimeBasisObj,SolutionObjCovar,lwd=2,NameLocation="Murano",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Murano(A)"]-SolutionObjCovar$BetaHat*MuranoCovar,col="red",pch=16)
dev.off()

# Cavallino-Treporti
png(filename="Cavallino-Treporti.png")
FixedPointPlot(12.51000,45.46500,SpaceBasisObj,TimeBasisObj,SolutionObjCovar,lwd=2,NameLocation="Cavallino-Treporti",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Cavallino-Treporti"]-SolutionObjCovar$BetaHat*CavallinoTreportiCovar,col="red",pch=16)
dev.off()

# Jesolo
png(filename="Jesolo.png")
FixedPointPlot(12.64139,45.54,SpaceBasisObj,TimeBasisObj,SolutionObjCovar,lwd=2,NameLocation="Jesolo",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Jesolo"]-SolutionObjCovar$BetaHat*JesoloCovar,col="red",pch=16)
dev.off()

# Caorle
png(filename="Caorle.png")
FixedPointPlot(12.88833,45.60250,SpaceBasisObj,TimeBasisObj,SolutionObjCovar,lwd=2,NameLocation="Caorle",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Caorle"]-SolutionObjCovar$BetaHat*CaorleCovar,col="red",pch=16)
dev.off()

# Chioggia
png(filename="Chioggia.png")
FixedPointPlot(12.27944,45.22056,SpaceBasisObj,TimeBasisObj,SolutionObjCovar,lwd=2,NameLocation="Chioggia",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Chioggia"]-SolutionObjCovar$BetaHat*ChioggiaCovar,col="red",pch=16)
dev.off()

# Portogruaro
png(filename="Portogruaro.png")
FixedPointPlot(12.83722,45.77722,SpaceBasisObj,TimeBasisObj,SolutionObjCovar,lwd=2,NameLocation="Portogruaro",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Portogruaro"]-SolutionObjCovar$BetaHat*PortogruaroCovar,col="red",pch=16)
dev.off()

# San Donà di Piave
png(filename="SanDonàdiPiave.png")
FixedPointPlot(12.56528,45.63389,SpaceBasisObj,TimeBasisObj,SolutionObjCovar,lwd=2,NameLocation="San Donà di Piave",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="SanDonàdiPiave"]-SolutionObjCovar$BetaHat*SanDonàCovar,col="red",pch=16)
dev.off()

# Cavarzere
png(filename="Cavarzere.png")
FixedPointPlot(12.08389,45.13667,SpaceBasisObj,TimeBasisObj,SolutionObjCovar,lwd=2,NameLocation="Cavarzere",ylim=lim)
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Cavarzere"]-SolutionObjCovar$BetaHat*CavarzereCovar,col="red",pch=16)
dev.off()

## CONFIDENCE INTERVAL
IC=ST.IC(Data,DesMat,SpaceBasisObj,TimeBasisObj,LogLambdaS,LogLambdaT,Alpha=0.05)

FixedTimePlot(1998,SpaceBasisObj,TimeBasisObj,SolutionObj)