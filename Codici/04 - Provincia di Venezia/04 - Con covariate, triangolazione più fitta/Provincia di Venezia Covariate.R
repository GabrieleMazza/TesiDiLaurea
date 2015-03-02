#ANALISI SPAZIO TEMPO PER LA PROVINCIA DI VENEZIA

source("SpazioTempo.R")
load("Territorio.RData")
library(rgl)
library(mgcv)
library(SDMTools)

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
for(i in 1:length(Codici[1:nint]))
{
    Data<-numeric(length(TimePoints))
    for(j in 1:length(TimePoints))
    {
        Data[j]=Risposta$TotalePC[(Risposta$Codice==Codici[i]) & (Risposta$Anno==TimePoints[j])]
    }
    DataMatrix<-c(DataMatrix,Data)
}

# Matrice Disegno

DesMat<-NULL
attach(CoordinateCovariate)
for(i in Codici[!is.na(Codici)])
{
    PL<-NULL
    for (j in TimePoints)
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
    DesMat<-c(DesMat,PL)
}
detach(CoordinateCovariate)
dimnames(DesMat)[[2]]<-NULL


##### GCV #####

LogS<-seq(-9,-8,length.out=5)
LogT<-seq(-1,+1,length.out=5)
GCVResult = ST.GCV.Covar(DataMatrix,DesMat,SpaceBasisObj,TimeBasisObj,LogS,LogT)

png(filename="GCV Matrix.png")
image(LogS,LogT,GCVResult$GCVMatrix,col=heat.colors(100),main="GCV Matrix",xlab="logLambdaS",ylab="logLambdaT")
dev.off()

LambdaS=10^GCVResult$Best[1]
LambdaT=10^GCVResult$Best[2]

save(file="GCVResult.RData",GCVResult,LogS,LogT)

# LambdaS=10^-7.888889
# LambdaT=10^-0.1111111

##### RISOLUZIONE DEL SISTEMA #####

SolutionObj<-ST.Smooth.Covar(DataMatrix,DesMat,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)




##### RISOLUZIONE DEL SISTEMA #####

SolutionObj<-ST.Smooth(Data,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)


# Ora salvo i risultati
write.table(SolutionObj$C,file="VettoreC.txt",row.names=FALSE,col.names=FALSE)


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

zlim=c(min(ResultFitted,na.rm=T),max(ResultFitted,na.rm=T))

for(j in 1:length(TimePoints))
{
    # Matrice con la reale
    Mat <- matrix(ResultFitted[,j],nrow=nx,ncol=ny,byrow=F)
    # Plot
    png(filename=paste("Anno",TimePoints[j],".png",sep=" "))
    image(xvec,yvec,Mat,zlim=zlim,main=paste("Funzione stimata tempo ",TimePoints[j],sep=""))
    lines(xbound,ybound,lwd=1)
    contour(xvec,yvec,Mat,nlevels=10,add=TRUE)
    dev.off()
}



##### PLOT IN ALCUNE SPECIFICHE CITTA' #####

attach(CoordinateCovariate)

VeneziaCovar<-NULL
SanMicheleAlTagliamentoCovar<-NULL
LidoCovar<-NULL

VeneziaCode<-425
SanMicheleAlTagliamentoCode<-417
LidoCode<-582

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
    VeneziaCovar<-c(VeneziaCovar,(get(name)[CoordinateCovariate$Codice==VeneziaCode]))
    SanMicheleAlTagliamentoCovar<-c(SanMicheleAlTagliamentoCovar,(get(name)[CoordinateCovariate$Codice==SanMicheleAlTagliamentoCode]))
    LidoCovar<-c(LidoCovar,(get(name)[CoordinateCovariate$Codice==LidoCode]))
    
}
detach(CoordinateCovariate)

# Plot per punto fissato
# Vediamo come va a venezia
png(filename="Venezia.png")
FixedPointPlot(12.327500,45.438056,SolutionObj,NameLocation = "Venezia")
points(1997:2011,(Risposta$TotalePC[Risposta$Comune=="Venezia"]-SolutionObj$BetaHat*VeneziaCovar),col="red",pch=16)
dev.off()

# Vediamo come va a San Michele al Tagliamento
png(filename="San Michele al Tagliamento.png")
FixedPointPlot(12.994722,45.767222,SolutionObj,NameLocation = "San Michele al Tagliamento")
points(1997:2011,(Risposta$TotalePC[Risposta$Comune=="SanMichelealTagliamento"]-SolutionObj$BetaHat*SanMicheleAlTagliamentoCovar),col="red",pch=16)
dev.off()

# Vediamo come va a Lido
png(filename="Lido(A).png")
FixedPointPlot(12.348115,45.384122,SolutionObj,NameLocation = "Lido")
points(1997:2011,(Risposta$TotalePC[Risposta$Comune=="Lido(A)"]-SolutionObj$BetaHat*LidoCovar),col="red",pch=16)
dev.off()




##### INTERVALLI DI CONFIDENZA #####

ICResult = ST.IC(DataMatrix,DesMat,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)
save(file="ICResult.RData",ICResult)
