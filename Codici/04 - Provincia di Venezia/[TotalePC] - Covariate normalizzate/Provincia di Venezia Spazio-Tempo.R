#ANALISI SPAZIO TEMPO PER LA PROVINCIA DI VENEZIA

source("SpazioTempo.R")
load("Territorio.RData")
load("Maps.RData")
library(rgl)
library(mgcv)
library(SDMTools)
library(fda)

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

# Matrici di Dati
DataMatrix<-NULL
for(j in TimePoints)
{
    Data<-numeric(nint)
    for(i in 1:length(Codici[1:nint]))
    {
        Data[i]=Risposta$TotalePC[(Risposta$Codice==Codici[i]) & (Risposta$Anno==j)]
    }
    DataMatrix<-cbind(DataMatrix,Data)
}

# Matrice Disegno

DesMat<-NULL
attach(CoordinateCovariate)
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
    PL<-NULL
    for (i in Codici[!is.na(Codici)])
    {
        PL<-c(PL,(get(name)[CoordinateCovariate$Codice==i]/(Risposta$TotalePC[(Risposta$Codice==i) & (Risposta$Anno==j)])))
        #PL<-c(PL,(get(name)[CoordinateCovariate$Codice==i]))
    }
    DesMat<-cbind(DesMat,PL)
}
detach(CoordinateCovariate)
dimnames(DesMat)[[2]]<-NULL


##### GCV #####

LogS<--10:0
LogT<--10:0
GCVResult = ST.GCV.Covar(DataMatrix,DesMat,SpaceBasisObj,TimeBasisObj,LogS,LogT)

png(filename="GCV Matrix.png")
image(LogS,LogT,GCVResult$GCVMatrix,col=heat.colors(100),main="GCV Matrix",xlab="logLambdaS",ylab="logLambdaT")
dev.off()

LambdaS=10^GCVResult$Best[1]
LambdaT=10^GCVResult$Best[2]

save.image()


##### RISOLUZIONE DEL SISTEMA #####

SolutionObj<-ST.Smooth.Covar(DataMatrix,DesMat,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)

# Ora salvo i risultati
write.table(SolutionObj$C,file="MatriceC.txt",row.names=FALSE,col.names=FALSE)
write.table(SolutionObj$BetaHat,file="BetaHat.txt",row.names=FALSE,col.names=FALSE)

nx<-100
ny<-100
xvec<-seq(min(x),max(x),length.out=nx)
yvec<-seq(min(y),max(y),length.out=ny)
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))

ResultFitted<-NULL

for(j in TimePoints)
{
    Time<-rep(j,length(xx))
    Result<-eval.ST.fd(xx,yy,Time,SolutionObj)
    
    ResultFitted<-cbind(ResultFitted,Result)
}

zlim=c(min(ResultFitted,na.rm=T),max(ResultFitted,na.rm=T))

for(j in 1:length(TimePoints))
{
    # Matrice con la reale
    Mat <- matrix(ResultFitted[,j],nrow=nx,ncol=ny,byrow=F)
    # Plot
    png(filename=paste("Anno",TimePoints[j],".png",sep=" "))
    image(xvec,yvec,Mat,zlim=zlim,main=paste("Funzione reale tempo ",TimePoints[j],sep=""))
    lines(x[!InternalPoints],y[!InternalPoints],lwd=1)
    contour(xvec,yvec,Mat,nlevels=10,add=TRUE)
    dev.off()
}

save.image()