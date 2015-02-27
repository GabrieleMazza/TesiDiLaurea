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

SolutionObj=ReadSolutionObjCovar(SpaceBasisObj,TimeBasisObj,"MatriceC.txt","BetaHat.txt")



##### GRAFICI ANNO PER ANNO #####

nx<-200
ny<-200
xvec<-seq(12.3,12.4,length.out=nx)
yvec<-seq(45.47,45.6,length.out=ny)
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))

j<-1997
Time<-rep(j,length(xx))
Result<-ST.Eval(xx,yy,Time,SolutionObj)

zlim=c(min(Result,na.rm=T),max(Result,na.rm=T))

j<-1
# Matrice
Mat <- matrix(Result,nrow=nx,ncol=ny,byrow=F)
# Plot
png(filename=paste("Zoom ",TimePoints[j],".png",sep=""))
image(xvec,yvec,Mat,zlim=zlim,main=paste("Funzione stimata zoom tempo ",TimePoints[j],sep=""))
lines(x[!InternalPoints],y[!InternalPoints],lwd=1)
contour(xvec,yvec,Mat,nlevels=10,add=TRUE)
dev.off()