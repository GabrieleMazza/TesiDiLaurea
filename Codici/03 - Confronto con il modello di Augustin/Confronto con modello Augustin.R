load("FrontieraC.RData")
load("Validazione.RData")

source("SpazioTempo.R")
library(mgcv)
library(SDMTools)
library(MASS)
library(fda)

# Devo aggiungere del rumore al dato?
noise<-TRUE

# Per quanti istanti di tempo eseguire la stima?
TimePoints<-0:5

# Distinguo i punti interni dai punti di bordo
xknot<-x[TypePoint]
yknot<-y[TypePoint]
xbound<-x[!TypePoint]
ybound<-y[!TypePoint]

# FUNZIONE PER LA PERTURBAZIONE
fun=function(x,y,t,xbound,ybound)
{
    Bound<-cbind(xbound,ybound)
    Result<-numeric(length(x))
    
    for(i in 1:length(x))
    {
        if(pnt.in.poly(cbind(x[i],y[i]),Bound)$pip==0)
        {
            Result[i]=NA
        } else
        {
            Result[i]=(x[i]+y[i]+t[i])^2
        }
    }
    
    return(Result)
}

# A che tempo devo farre la previsione con i modelli?
TimePrevision<-2.5

# Tutti i punti della validazione sono contenuti all'interno del bordo?
Bound<-cbind(xbound,ybound)
sum(pnt.in.poly(cbind(xValid,yValid),Bound)$pip)==length(xValid)




##### APPLICO IL GAM DI AUGUSTIN #####

Data<-NULL
TimeData<-NULL
xData<-NULL
yData<-NULL

for(i in TimePoints)
{
    Time<-rep(i,length(x[TypePoint]))
    if(noise)
    {
        Data<-c(Data,(fun(xknot,yknot,Time,xbound,ybound)+rnorm(length(xknot),0,0.005)))
    } else
    {
        Data<-c(Data,fun(xknot,yknot,Time,xbound,ybound))
    }
    TimeData<-c(TimeData,Time)
    xData<-c(xData,xknot)
    yData<-c(yData,yknot)
}


#Definisco il bordo
# knots<-data.frame(xData=xknot,yData=yknot)
# fsb<-list(xData=xbound,yData=ybound)
# fsb<-list(fsb)
# xt=list(bnd=fsb)
# knots=knots

mod <-gam(Data~
              te(xData,yData,TimeData,d=c(2,1),bs=c("tp","cr")),
          method="GCV.Cp")

summary(mod)

# Ora previsione
NewPoints<-data.frame(xData=xValid,yData=yValid,TimeData=rep(TimePrevision,length(xValid)))
PredictionGAM<-as.vector(predict(mod,NewPoints))



##### APPLICO IL MODELLO SPAZIOTEMPO #####

# Studio il valore reale della funzione visto che posso conoscerlo...
TimeBasisObj<-Create.Bspline.Time.Basis(TimePoints,TimeOrder=4,DerivativeOrder=2,PlotIt=F)
SpaceBasisObj<-Create.FEM.Space.Basis(cbind(x,y),Triang,TypePoint,1)

# Il valore ottimo dei lambda era -6 e -6
DataMatrix<-matrix(data=Data,nrow=length(xknot),ncol=length(TimePoints),byrow=F)

LogS<-seq(-7,-5,length.out=5)
LogT<-seq(-7,-5,length.out=5)
GCVResult<-ST.GCV(DataMatrix,SpaceBasisObj,TimeBasisObj,LogS,LogT)

LambdaS=10^GCVResult$Best[1]
LambdaT=10^GCVResult$Best[2]

save(file="GCVResult.RData",GCVResult)

# Migliore (finora) LogLambdaS =
# Migliore (finora) LogLambdaT =

SolutionObj<-ST.Smooth(DataMatrix,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)

PredictionST<-ST.Eval(xValid,yValid,rep(TimePrevision,length(xValid)),SolutionObj)




##### VALORE REALE DELLA FUNZIONE #####

# Studio il valore reale della funzione visto che posso conoscerlo...
RealValues<-fun(xValid,yValid,rep(TimePrevision,length(xValid)),xbound,ybound)





##### CONFRONTI TRA I METODI #####

# Iniziamo da un plot dei valori assunti...
plot(RealValues,type='l')
points(PredictionGAM,type='l',col="red")
points(PredictionST,type='l',col="blue")

# Si sovrappongono
# Traccio le funzioni...
Max<-max(RealValues,PredictionGAM,PredictionST)
Min<-min(RealValues,PredictionGAM,PredictionST)
zlim=c(Min,Max)

PlotMatrix=matrix(ncol=length(yvec),nrow=length(xvec))

# REALE
for(i in 1:(dim(PosMatrix)[1]))
{
    PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=RealValues[i]
}
# Plot
png(filename=paste("Funzione Reale al tempo ",TimePrevision,".png",sep=""))
image(xvec,yvec,PlotMatrix,zlim=zlim,main=paste("Funzione Reale al tempo",TimePrevision,sep=" "))
lines(x[TypePoint==FALSE],y[TypePoint==FALSE],lwd=3)
contour(xvec,yvec,PlotMatrix,nlevels=10,add=TRUE)
dev.off()


# GAM
for(i in 1:(dim(PosMatrix)[1]))
{
    PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=PredictionGAM[i]
}
# Plot
png(filename=paste("Funzione GAM al tempo ",TimePrevision,".png",sep=""))
image(xvec,yvec,PlotMatrix,zlim=zlim,main=paste("Funzione stimata da GAM al tempo",TimePrevision,sep=" "))
lines(x[TypePoint==FALSE],y[TypePoint==FALSE],lwd=3)
contour(xvec,yvec,PlotMatrix,nlevels=10,add=TRUE)
dev.off()


# ST
for(i in 1:(dim(PosMatrix)[1]))
{
    PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=PredictionST[i]
}
# Plot
png(filename=paste("Funzione ST al tempo ",TimePrevision,".png",sep=""))
image(xvec,yvec,PlotMatrix,zlim=zlim,main=paste("Funzione stimata da ST al tempo",TimePrevision,sep=" "))
lines(x[TypePoint==FALSE],y[TypePoint==FALSE],lwd=3)
contour(xvec,yvec,PlotMatrix,nlevels=10,add=TRUE)
dev.off()


# Studio degli errori
ErrorGAM=(PredictionGAM-RealValues)^2
ErrorST=(PredictionST-RealValues)^2

sum(ErrorGAM)
sum(ErrorST)
