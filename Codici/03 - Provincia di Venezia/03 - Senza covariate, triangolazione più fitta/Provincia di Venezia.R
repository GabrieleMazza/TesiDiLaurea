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

# LogS = seq(-10,-9,by=0.125)
# LogT = seq(-0.5,0.5,by=0.125)
# GCVResult = ST.GCV(Data,SpaceBasisObj,TimeBasisObj,LogS,LogT)
# 
# png(filename="GCV Matrix.png")
# image(LogS,LogT,GCVResult$GCVMatrix,col=heat.colors(100),main="GCV Matrix",xlab="logLambdaS",ylab="logLambdaT")
# dev.off()
# 
# LambdaS=10^GCVResult$Best[1]
# LambdaT=10^GCVResult$Best[2]
# 
# save(file="GCVResult.RData",GCVResult,LogS,LogT)

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
    # Matrice con la reale
    Mat <- matrix(ResultFitted[,j],nrow=nx,ncol=ny,byrow=F)
    # Plot
    png(filename=paste("Anno",TimePoints[j],".png",sep=" "))
    image(xvec,yvec,Mat,zlim=zlim,main=paste("Funzione stimata tempo ",TimePoints[j],sep=""))
    lines(xbound,ybound,lwd=1)
    contour(xvec,yvec,Mat,nlevels=10,add=TRUE)
    dev.off()
}



##### GRAFICI NELLE CITTA' SCELTE #####


##### GRAFICI NELLE CITTA' SCELTE #####

# Plot per punto fissato
# Vediamo come va a venezia
png(filename="Venezia.png")
FixedPointPlot(12.327500,45.438056,SolutionObj,NameLocation="Venezia")
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Venezia"],col="red",pch=16)
dev.off()

# Vediamo come va a San Michele al Tagliamento
png(filename="San Michele al Tagliamento.png")
FixedPointPlot(12.994722,45.767222,SolutionObj,NameLocation="San Michele al Tagliamento")
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="SanMichelealTagliamento"],col="red",pch=16)
dev.off()

# Vediamo come va a Lido di Venezia
png(filename="Lido(A).png")
FixedPointPlot(12.348115,45.384122,SolutionObj,NameLocation="Lido di Venezia")
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Lido(A)"],col="red",pch=16)
dev.off()

# Vediamo come va a Pellestrina
png(filename="Pellestrina(A).png")
FixedPointPlot(12.30181,45.27324,SolutionObj,NameLocation="Pellestrina")
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Pellestrina(A)"],col="red",pch=16)
dev.off()

# Vediamo come va a Murano
png(filename="Murano(A).png")
FixedPointPlot(12.35155,45.45810,SolutionObj,NameLocation="Murano")
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Murano(A)"],col="red",pch=16)
dev.off()

# Vediamo come va a Cavallino-Treporti
png(filename="Cavallino-Treporti.png")
FixedPointPlot(12.51000,45.46500,SolutionObj,NameLocation="Cavallino-Treporti")
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Cavallino-Treporti"],col="red",pch=16)
dev.off()

# Vediamo come va a Jesolo
png(filename="Jesolo.png")
FixedPointPlot(12.64139,45.54,SolutionObj,NameLocation="Jesolo")
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Jesolo"],col="red",pch=16)
dev.off()

# Vediamo come va a Caorle
png(filename="Caorle.png")
FixedPointPlot(12.88833,45.60250,SolutionObj,NameLocation="Caorle")
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Caorle"],col="red",pch=16)
dev.off()

# Vediamo come va a Chioggia
png(filename="Chioggia.png")
FixedPointPlot(12.27944,45.22056,SolutionObj,NameLocation="Chioggia")
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Chioggia"],col="red",pch=16)
dev.off()

# Vediamo come va a Portogruaro
png(filename="Portogruaro.png")
FixedPointPlot(12.83722,45.77722,SolutionObj,NameLocation="Portogruaro")
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Portogruaro"],col="red",pch=16)
dev.off()

# Vediamo come va a San Donà di Piave
png(filename="SanDonàdiPiave.png")
FixedPointPlot(12.56528,45.63389,SolutionObj,NameLocation="San Donà di Piave")
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="SanDonàdiPiave"],col="red",pch=16)
dev.off()

# Vediamo come va a Cavarzere
png(filename="Cavarzere.png")
FixedPointPlot(12.08389,45.13667,SolutionObj,NameLocation="Cavarzere")
points(1997:2011,Risposta$TotalePC[Risposta$Comune=="Cavarzere"],col="red",pch=16)
dev.off()