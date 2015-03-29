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
Data<-NULL
for(i in 1:length(Codici[1:nint]))
{
    DataTmp<-numeric(length(TimePoints))
    for(j in 1:length(TimePoints))
    {
        DataTmp[j]=Risposta$TotalePC[(Risposta$Codice==Codici[i]) & (Risposta$Anno==TimePoints[j])]
    }
    Data<-c(Data,DataTmp)
}


##### GCV #####

LogS = seq(-9.125,-8.875,by=0.125)
LogT = seq(-0.25,0,by=0.125)
GCVResult = ST.GCV(Data,SpaceBasisObj,TimeBasisObj,LogS,LogT)

png(filename="GCV Matrix.png")
image(LogS,LogT,GCVResult$GCVMatrix,col=heat.colors(100),main="GCV Matrix",xlab="logLambdaS",ylab="logLambdaT")
dev.off()

LambdaS=10^GCVResult$Best[1]
LambdaT=10^GCVResult$Best[2]

save(file="GCVResult.RData",GCVResult,LogS,LogT)

LambdaS=10^-9.625
LambdaT=10^0





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
