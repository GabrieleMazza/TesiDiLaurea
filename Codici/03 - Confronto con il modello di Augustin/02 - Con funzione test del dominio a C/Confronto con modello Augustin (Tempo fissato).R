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

# Function for perturbation
fun=function(x)
{
    y=cos(x)
    return(y)
}

# A che tempo devo farre la previsione con i modelli?
TimePrevision<-2.5

# Tutti i punti della validazione sono contenuti all'interno del bordo?
Bound<-cbind(xbound,ybound)
sum(pnt.in.poly(cbind(xValid,yValid),Bound)$pip)==length(xValid)

# Quante volte eseguire il ciclo?
NTIMES=10


##### GCV #####

# ORA GCV, LO APPLICO AI DATI NON PERTURBATI, NON VOGLIO APPLICARLO TUTTE LE VOLTE NEL FOR

TimeBasisObj<-Create.Bspline.Time.Basis(TimePoints,TimeOrder=4,DerivativeOrder=2,PlotIt=F)
SpaceBasisObj<-Create.FEM.Space.Basis(cbind(x,y),Triang,TypePoint,1)

DataMatrix<-NULL
TimeData<-NULL
xData<-NULL
yData<-NULL
for(j in TimePoints)
{
    DataMatrix<-cbind(DataMatrix,Data*fun(j))
    
    #Per dopo..
    TimeData<-c(TimeData,rep(j,length(x[TypePoint])))
    xData<-c(xData,xknot)
    yData<-c(yData,yknot)
}

LogS<-seq(-20,-10,length.out=20)
LogT<-seq(-20,-10,length.out=20)
GCVResult<-ST.GCV(DataMatrix,SpaceBasisObj,TimeBasisObj,LogS,LogT)
png(filename="GCV Matrix.png")
image(LogS,LogT,GCVResult$GCVMatrix,col=heat.colors(100),main="GCV Matrix",xlab="logLambdaS",ylab="logLambdaT")
dev.off()

save(file="GCVResult.RData",GCVResult,LogS,LogT)



##### CICLO PER IL CONFRONTO DELL'ALGORITMO #####

PredictionSTMat<-NULL
PredictionGAMMat<-NULL
ErrorGAM<-NULL
ErrorST<-NULL

RealValues<-fs.test(xValid,yValid)*fun(TimePrevision)

for(ntimes in 1:NTIMES)
{
    #Ricavo i dati
    DataVec<-NULL
    for(i in TimePoints)
    {
        if(noise)
        {
            DataVec<-c(DataVec,Data*fun(i)+rnorm(length(xknot),0,0.05))
        } else
        {
            DataVec<-c(DataVec,Data*fun(i))
        }
    }
    
    ## MODELLO AUGUSTIN ##
    mod <-gam(DataVec~
                  te(xData,yData,TimeData,d=c(2,1),bs=c("tp","cr")),
              method="GCV.Cp")
    
    summary(mod)
    # Ora previsione
    NewPoints<-data.frame(xData=xValid,yData=yValid,TimeData=rep(TimePrevision,length(xValid)))
    PredictionGAM<-as.vector(predict(mod,NewPoints))
    PredictionGAMMat<-cbind(PredictionGAMMat,PredictionGAM)
    
    ## MODELLO SPAZIOTEMPO ##
    DataMatrix<-matrix(data=DataVec,nrow=length(xknot),ncol=length(TimePoints),byrow=F)
    LambdaS=10^GCVResult$Best[1]
    LambdaT=10^GCVResult$Best[2]
    # Risolvo
    SolutionObj<-ST.Smooth(DataMatrix,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)
    # Previsione
    PredictionST<-ST.Eval(xValid,yValid,rep(TimePrevision,length(xValid)),SolutionObj)
    PredictionSTMat<-cbind(PredictionSTMat,PredictionST)
    
    # Studio degli errori
    ErrorGAM<-c(ErrorGAM,sum((PredictionGAM-RealValues)^2))
    ErrorST<-c(ErrorST,sum((PredictionST-RealValues)^2))
    
}

png(filename="Confronto tra i metodi.png")
plot(ErrorGAM,ylim=c(min(ErrorGAM,ErrorST),max(ErrorGAM,ErrorST)),main=paste("Confronto tra metodi al tempo fissato ",TimePrevision,sep=""),,type='l')
points(ErrorST,type='l',col="red")
legend("bottomright",col=c("black","red"),legend=c("GAM","ST"),lty=1)
dev.off()

save(file="Risultati.RData",PredictionSTMat,PredictionGAMMat,RealValues,ErrorGAM,ErrorST)




##### PLOT DELL'ULTIMO CASO #####

#Provo un plot, ultimo istante provato
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
image(xvec,yvec,PlotMatrix,zlim=zlim,main=paste("Funzione Reale al tempo",TimePrevision,"(ultimo caso)",sep=" "))
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
image(xvec,yvec,PlotMatrix,zlim=zlim,main=paste("Funzione stimata da GAM al tempo",TimePrevision,"(ultimo caso)",sep=" "))
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
image(xvec,yvec,PlotMatrix,zlim=zlim,main=paste("Funzione stimata da ST al tempo",TimePrevision,"(ultimo caso)",sep=" "))
lines(x[TypePoint==FALSE],y[TypePoint==FALSE],lwd=3)
contour(xvec,yvec,PlotMatrix,nlevels=10,add=TRUE)
dev.off()

