#ANALISI SPAZIO TEMPO PER LA PROVINCIA DI VENEZIA

source("SpazioTempo.R")
load("Territorio.RData")
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





##### GCV #####
# 
# LogS<-seq(-11,-9,length.out=10)
# LogT<-seq(-1,1,length.out=10)
# GCVResult = ST.GCV(DataMatrix,SpaceBasisObj,TimeBasisObj,LogS,LogT)
# 
# png(filename="GCV Matrix.png")
# image(LogS,LogT,GCVResult$GCVMatrix,col=heat.colors(100),main="GCV Matrix",xlab="logLambdaS",ylab="logLambdaT")
# dev.off()
# 
# LambdaS=10^GCVResult$Best[1]
# LambdaT=10^GCVResult$Best[2]
# 
# save(file="GCVResult.RData",GCVResult,LogS,LogT)

LambdaS=10^-9.6666667
LambdaT=10^-0.1111111

##### RISOLUZIONE DEL SISTEMA #####

SolutionObj<-ST.Smooth(DataMatrix,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)

# Ora salvo i risultati
write.table(SolutionObj$C,file="MatriceC.txt",row.names=FALSE,col.names=FALSE)




##### GRAFICI ANNO PER ANNO #####


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
    Result<-ST.Eval(xx,yy,Time,SolutionObj)
    
    ResultFitted<-cbind(ResultFitted,Result)
}

zlim=c(min(ResultFitted,na.rm=T),max(ResultFitted,na.rm=T))

for(j in 1:length(TimePoints))
{
    # Matrice con la reale
    Mat <- matrix(ResultFitted[,j],nrow=nx,ncol=ny,byrow=F)
    # Plot
    png(filename=paste("Anno",TimePoints[j],".png",sep=" "))
    image(xvec,yvec,Mat,zlim=zlim,main=paste("Funzione stimata tempo ",TimePoints[j],sep=""))
    lines(x[!InternalPoints],y[!InternalPoints],lwd=1)
    contour(xvec,yvec,Mat,nlevels=10,add=TRUE)
    dev.off()
}



##### GRAFICI NELLE CITTA' SCELTE #####

# Plot per punto fissato
# Vediamo come va a venezia
png(filename="Venezia.png")
FixedPointPlot(12.327500,45.438056,SolutionObj,NameLocation="Venezia")
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Venezia"],col="red")
dev.off()

# Vediamo come va a San Michele al Tagliamento
png(filename="San Michele al Tagliamento.png")
FixedPointPlot(12.994722,45.767222,SolutionObj,NameLocation="San Michele al Tagliamento")
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="SanMichelealTagliamento"],col="red")
dev.off()

# Vediamo come va a Lido
png(filename="Lido(A).png")
FixedPointPlot(12.348115,45.384122,SolutionObj,NameLocation="Lido")
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Lido(A)"],col="red")
dev.off()
