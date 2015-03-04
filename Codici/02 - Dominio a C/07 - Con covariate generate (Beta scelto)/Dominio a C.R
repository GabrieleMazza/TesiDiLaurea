library(rgl)
library(mgcv)
library(SDMTools)
library(ggplot2)
source("SpazioTempo.R")
load("FrontieraC.RData")
load("Validazione.RData")

# Devo aggiungere del rumore al dato?
noise<-TRUE

# Per quanti istanti di tempo eseguire la stima?
TimePoints<-0:5


# DEVIAZIONE STANDARD DEL RUMORE
SigmaNoise=0.05

# Per la covariata
Mu<-1
Sigma<-0.5
Beta<-1


# Distinguo i punti interni dai punti di bordo
xknot<-x[TypePoint]
yknot<-y[TypePoint]
xbound<-x[!TypePoint]
ybound<-y[!TypePoint]

# Function for perturbation
fun=function(x,y,t)
{
    y=fs.test(x,y)*cos(t)
    return(y)
}

# Tutti i punti della validazione sono contenuti all'interno del bordo?
Bound<-cbind(xbound,ybound)
sum(pnt.in.poly(cbind(xValid,yValid),Bound)$pip)==length(xValid)




##### PERTURBAZIONE TEMPORALE DEI DATI #####

#Ricavo i dati
DataMatrix<-NULL
for(t in TimePoints)
{
    if(noise)
    {
        DataMatrix<-cbind(DataMatrix,fun(xknot,yknot,rep(t,length(xknot)))+rnorm(length(xknot),0,SigmaNoise))
    } else
    {
        DataMatrix<-cbind(DataMatrix,fun(xknot,yknot,j))
    }
}
DataFUN<-NULL
for(i in 1:(dim(DataMatrix)[1]))
{
    DataFUN<-c(DataFUN,DataMatrix[i,])
}


##### MATRICE DISEGNO #####

# Voglio la covariata non fissa nel tempo
# Genero la covariata da una normale
DesMat<-rnorm(length(DataFUN),Mu,Sigma)

# Li aggiungo poi ai dati secondo il Beta reale
Data=DataFUN+Beta*DesMat




##### BUBBLEPLOT #####

# Ho bisogno dei un dataframe
Knot=data.frame(xknot=c(xknot,10,10),yknot=c(yknot,1,-1))
Boundary<-data.frame(xbound,ybound)

for(j in 1:length(TimePoints))
{
    Gen<-NULL
    for(i in 0:(dim(DataMatrix)[1]-1))
    {
        Gen<-c(Gen,Data[i*length(TimePoints)+j])
    }
    Gen<-c(Gen,min(DataMatrix),max(DataMatrix))
    
    ggplot(data=Knot, aes(x=xknot, y=yknot)) +
        geom_point(aes(size=(Gen),color=Gen)) +
        scale_color_gradient(low="yellow", high="red") +
        geom_polygon(data=Boundary,aes(x=xbound,y=ybound),alpha=1,colour="black", fill=NA, size=1.1) +
        theme_bw() +
        labs(color="Dato")+
        scale_x_continuous(name="x", limits=c(-1,3.5)) +
        scale_y_continuous(name="y",limits=c(-1,1)) +
        labs(title=paste("Dati Generati al tempo ",trunc(TimePoints[j],2),sep="")) + 
        theme(plot.title = element_text(size = rel(1.5), face="bold")) +
        guides(size=FALSE)
    ggsave(file=paste("Dati Generati al tempo ",TimePoints[j],".png",sep=""))
}

# NOTA BENE
# I warning sono previsti... Ho aggiunto due punti, uno con il massimo e uno con il
# minimo di tutti i grafici, per uniformare la scala
rm(Knot,Boundary)

##### GCV #####

#Creo le basi in spazio e tempo
TimeBasisObj<-Create.Bspline.Time.Basis(TimePoints,TimeOrder=4,DerivativeOrder=2,PlotIt=F)
SpaceBasisObj<-Create.FEM.Space.Basis(cbind(x,y),Triang,TypePoint,1)

LogS<--8:+1
LogT<--8:+1
GCVResult<-ST.GCV.Covar(Data,DesMat,SpaceBasisObj,TimeBasisObj,LogS,LogT)

png(filename="GCV Matrix.png")
image(LogS,LogT,GCVResult$GCVMatrix,col=heat.colors(100),main="GCV Matrix",xlab="logLambdaS",ylab="logLambdaT")
dev.off()

LambdaS=10^GCVResult$Best[1]
LambdaT=10^GCVResult$Best[2]

save(file="GCVResult.RData",GCVResult,LogS,LogT)

# Con Normale(1,0.25), SigmaNoise=0.05, Beta=1
# LambdaS=10^-2
# LambdaT=10^-8

# Con Normale(1,1), SigmaNoise=0.05, Beta=1
# LambdaS=10^-3
# LambdaT=10^-5

# Con Normale(5,1), SigmaNoise=0.05, Beta=1
# LambdaS=10^-3
# LambdaT=10^-5

# Con Normale(0,1), SigmaNoise=1, Beta=1
# LambdaS=10^1
# LambdaT=10^-3


##### RISOLUZIONE DEL SISTEMA #####

SolutionObj<-ST.Smooth.Covar(Data,DesMat,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)

# Ora salvo i risultati
write.table(SolutionObj$C,file="MatriceC.txt",row.names=FALSE,col.names=FALSE)
write.table(SolutionObj$BetaHat,file="BetaHat.txt",row.names=FALSE,col.names=FALSE)



##### GRAFICI ANNO PER ANNO SOLO DELLA FUNZIONE SENZA COVARIATE #####

PlotMatrix=matrix(ncol=length(yvec),nrow=length(xvec))
           
# Voglio ricreare perfettamente le posizioni della image matrix tru

# RICERCA DEL MASSIMO
# PER AVERE UNA SCALA DI COLORI PER IL GRAFICO UNIFORME IN TUTTI GLI ISTANTI DI TEMPO

ResultFitted<-NULL
        
for(j in TimePoints)
{
    Time<-rep(j,length(xValid))
    Result<-ST.Eval(xValid,yValid,Time,SolutionObj)
            
    ResultFitted<-cbind(ResultFitted,Result)
}

ResultReal<-NULL

for(j in TimePoints)
{
    Time<-rep(j,length(xValid))
    Result<-fun(xValid,yValid,Time)
    ResultReal<-cbind(ResultReal,Result)
}

Max=max(ResultFitted,ResultReal,na.rm=TRUE)
Min=min(ResultFitted,ResultReal,na.rm=TRUE)
zlim<-c(Min,Max)
        

# GRAFICI

# Ora salvo tutti i grafici
# Voglio ricreare perfettamente le posizioni della image matrix tru

nx<-length(xvec)
ny<-length(yvec)
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))

for(j in 1:length(TimePoints))
{
    # Funzione Reale
    for(i in 1:(dim(PosMatrix)[1]))
    {
        PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=ResultReal[i,j]
    }
    # Plot
    png(filename=paste("Tempo ",trunc(TimePoints[j],2)," reale.png",sep=""))
    image(xvec,yvec,PlotMatrix,zlim=zlim,main=paste("Funzione reale tempo ",trunc(TimePoints[j],2),sep=""))
    lines(x[TypePoint==FALSE],y[TypePoint==FALSE],lwd=3)
    contour(xvec,yvec,PlotMatrix,nlevels=10,add=TRUE)
    dev.off()
    
    # Matrice con la stimata
    for(i in 1:(dim(PosMatrix)[1]))
    {
        PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=ResultFitted[i,j]
    }
    # Plot
    png(filename=paste("Tempo ",trunc(TimePoints[j],2)," stimata.png",sep=""))
    image(xvec,yvec,PlotMatrix,zlim=zlim,main=paste("Funzione stimata tempo ",trunc(TimePoints[j],2),sep=""))
    lines(x[TypePoint==FALSE],y[TypePoint==FALSE],lwd=3)
    contour(xvec,yvec,PlotMatrix,nlevels=10,add=TRUE)
    dev.off()
}

png(filename=paste("Plot per un punto fissato.png",sep=" "))
FixedPointPlot(xknot[1],yknot[2],SolutionObj)
points(TimePoints,fun(rep(xknot[1],length(TimePoints)),rep(yknot[1],length(TimePoints)),TimePoints),pch=16,col="red")
dev.off()


##### INTERVALLO DI CONFIDENZA #####
ICResult<-ST.IC(Data,DesMat,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)
save(file="ICResult.RData",ICResult,Beta)
