load("FrontieraC.RData")
load("Validazione.RData")

source("SpazioTempo.R")
library(mgcv)
library(SDMTools)
library(MASS)
library(ggplot2)
library(gstat)
library(spacetime)
library(sp)

set.seed(1)

# Devo aggiungere del rumore al dato?
noise<-TRUE

# Per quanti istanti di tempo eseguire la stima?
TimePoints<-seq(0,2*pi,length.out=9)

# DEVIAZIONE STANDARD DEL RUMORE
SigmaNoise=0.5

# Per la covariata
Mu<-0
Sigma<-1
Beta<-1

# Distinguo i punti interni dai punti di bordo
xknot<-x[TypePoint]
yknot<-y[TypePoint]
xbound<-x[!TypePoint]
ybound<-y[!TypePoint]

# Per il kriginig, preparazione degli oggetti
punti <- SpatialPoints(cbind(x[TypePoint],y[TypePoint]))  # object of class SpatialPoints
istanti <- as.POSIXct("2010-08-05")+3600*TimePoints


# Function for perturbation
fun=function(x,y,t)
{
    y=fs.test(x,y)*cos(t)
    return(y)
}

# Quante volte eseguire il ciclo?
NTIMES=50

# DIMENSIONI DELLA GRIGLIA
NX=50
NY=30
NT=20


##### OGGETTI UTILI #####

TimeBasisObj<-Create.Bspline.Time.Basis(TimePoints,TimeOrder=4,DerivativeOrder=2,PlotIt=F)
SpaceBasisObj<-Create.FEM.Space.Basis(cbind(x,y),Triang,TypePoint,1)

DataMatrix<-NULL
TimeDataMatrix<-NULL
xDataMatrix<-NULL
yDataMatrix<-NULL
for(t in TimePoints)
{
    if(noise)
    {
        DataMatrix<-cbind(DataMatrix,fun(xknot,yknot,rep(t,length(xknot)))+rnorm(length(xknot),0,SigmaNoise))
    } else
    {
        DataMatrix<-cbind(DataMatrix,fun(xknot,yknot,rep(t,length(xknot))))
    }
    
    #Per dopo..
    TimeDataMatrix<-cbind(TimeDataMatrix,rep(t,length(xknot)))
    xDataMatrix<-cbind(xDataMatrix,xknot)
    yDataMatrix<-cbind(yDataMatrix,yknot)
}


TimeData<-NULL
xData<-NULL
yData<-NULL
Data<-NULL
for(i in 1:(dim(TimeDataMatrix)[1]))
{
    Data<-c(Data,DataMatrix[i,])
    TimeData<-c(TimeData,TimeDataMatrix[i,])
    xData<-c(xData,xDataMatrix[i,])
    yData<-c(yData,yDataMatrix[i,])
}

rm(DataMatrix,xDataMatrix,yDataMatrix,TimeDataMatrix)




##### ESEMPIO DI BUBBLEPLOT DEI DATI #####

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

DesMatMat<-NULL
for(t in TimePoints)
{
    DesMatMat<-cbind(DesMatMat,rnorm(length(xknot),Mu,Sigma)) 
}

DataMatrix<-DataMatrix+Beta*DesMatMat

#Per Dopo
Data<-NULL
for(i in 1:(dim(DataMatrix)[1]))
{
    Data<-c(Data,DataMatrix[i,])
}

DesMat<-NULL
for(i in 1:(dim(DesMatMat)[1]))
{
    DesMat<-c(DesMat,DesMatMat[i,])
}

# Ho bisogno dei un dataframe
Knot=data.frame(xknot=c(xknot,10,10),yknot=c(yknot,1,-1))
Boundary<-data.frame(xbound,ybound)

for(j in 1: length(TimePoints))
{
    Gen<-c(DataMatrix[,j],max(DataMatrix),min(DataMatrix))
    
    ggplot(data=Knot, aes(x=xknot, y=yknot)) +
        geom_point(aes(size=(Gen),color=Gen)) +
        scale_color_gradient(low="yellow", high="red") +
        geom_polygon(data=Boundary,aes(x=xbound,y=ybound),alpha=1,colour="black", fill=NA, size=1.1) +
        theme_bw() +
        labs(color="Dato")+
        scale_x_continuous(name="x", limits=c(-1,3.5)) +
        scale_y_continuous(name="y",limits=c(-1,1)) +
        labs(title=paste("Dati Generati al tempo ",round(TimePoints[j],2),sep="")) + 
        theme(plot.title = element_text(size = rel(1.5), face="bold")) +
        guides(size=FALSE)
    ggsave(file=paste("Dati Generati al tempo ",round(TimePoints[j],2),".png",sep=""))
}

# NOTA BENE
# I warning sono previsti... Ho aggiunto due punti, uno con il massimo e uno con il
# minimo di tutti i grafici, per uniformare la scala
rm(Knot,Boundary)




##### GCV #####

# LogS<--8:1
# LogT<--8:1
# GCVResult<-ST.GCV.Covar(Data,DesMat,SpaceBasisObj,TimeBasisObj,LogS,LogT)
# png(filename="GCV Matrix.png")
# image(LogS,LogT,GCVResult$GCVMatrix,col=heat.colors(100),main="GCV Matrix",xlab="logLambdaS",ylab="logLambdaT")
# dev.off()
# 
# save(file="GCVResult.RData",GCVResult,LogS,LogT)
# 
# LambdaS=10^GCVResult$Best[1]
# LambdaT=10^GCVResult$Best[2]

LambdaS=10^-0.5
LambdaT=10^-3.25



##### PUNTI DI VALIDAZIONE, PER LA GRIGLIA SPAZIOTEMPO #####

fsb <- list(fs.boundary())
Bound<-cbind(xbound,ybound)


# Genero i dati per la previsione
xGen<-seq(-1,3.5,length.out=NX)
yGen<-seq(-1,1,length.out=NY)
tGen<-seq(min(TimePoints),max(TimePoints),length.out = NT)

#Ho bisogno di un vettore con i punti di validazione
xValid<-NULL
yValid<-NULL
for(i in 1:length(xGen))
{
    for(j in 1:length(yGen))
    {
        x=xGen[i]
        y=yGen[j]
        if((pnt.in.poly(cbind(x,y),Bound)$pip==1) & inSide(fsb,x=x,y=y))
        {
            xValid<-c(xValid,x)
            yValid<-c(yValid,y)
        }
    }
}

tValid<-NULL
for (t in tGen)
{
    tValid<-c(tValid,rep(t,length(xValid)))    
}
xValid<-rep(xValid,length(tGen))
yValid<-rep(yValid,length(tGen))

#Valori reali sulla griglia
RealValues<-fun(xValid,yValid,tValid)


# PER IL KRIGING
# definisco griglia spaziotempo pred.grid_C solo sul dominio a C su cui calcolare errore

# t_grid<-rep(seq(from=min(TimePoints),to=max(TimePoints),length.out=NT),each=NY*NX)
# n_t_grid <- length(t_grid)
# x_grid<-rep(rep(seq(from=-1,to=3.5,length=NX),times=NY),NT)
# y_grid<-rep(rep(seq(from=-1,to=1,length=NY),each=NX),NT)
# grd <- SpatialPoints(
#     matrix(data = c(x_grid,y_grid),
#            nrow = NX*NY, ncol = 2, byrow = FALSE)
# )
# tgrd <- as.POSIXct("2010-08-05")+3600*seq(from=min(TimePoints),to=max(TimePoints),length.out=NT)
# pred.grid_C <- STF(sp=grd, time=tgrd)
# 
# 


##### CICLO DEI MODELLI #####

RMSE<-NULL
LABEL<-NULL
BETA<-NULL

#Ricavo la frontiera da usare nel caso in cui ci sia Soap Film Smoothing nel GAM
fsb=list(list(xbound,ybound))
names(fsb[[1]]) <- c("xData","yData")
knots <- data.frame(xData=xknot,yData=yknot)
fsbReal <- list(fs.boundary())

n=length(RealValues)

# CICLO
for(ntimes in 1:NTIMES)
{
    print(paste("Simulazione ",ntimes,"/",NTIMES,sep=""))
    #Ricavo i dati
    DataMatrix<-NULL
    for(t in TimePoints)
    {
        if(noise)
        {
            DataMatrix<-cbind(DataMatrix,fun(xknot,yknot,rep(t,length(xknot)))+rnorm(length(xknot),0,SigmaNoise))
        } else
        {
            DataMatrix<-cbind(DataMatrix,fun(xknot,yknot,rep(t,length(xknot))))
        }
    }
    
    DesMatMat<-NULL
    for(t in TimePoints)
    {
        DesMatMat<-cbind(DesMatMat,rnorm(length(xknot),Mu,Sigma)) 
    }
    
    DataMatrix<-DataMatrix+Beta*DesMatMat
    
    #Per Dopo
    Data<-NULL
    for(i in 1:(dim(DataMatrix)[1]))
    {
        Data<-c(Data,DataMatrix[i,])
    }
    
    DesMat<-NULL
    for(i in 1:(dim(DesMatMat)[1]))
    {
        DesMat<-c(DesMat,DesMatMat[i,])
    }
    
    CovarValid<-rep(0,length(xValid))
    NewPoints<-data.frame(xData=xValid,yData=yValid,DesMat=CovarValid,TimeData=tValid)
        
    ## MODELLI TIPO AUGUSTIN ##
    
    # Primo caso
    #   Prodotto tensoriale di tps e cr
    mod <-gam(Data~
                  DesMat+
                  te(xData,yData,TimeData,k=c(NA,length(TimePoints)),d=c(2,1),bs=c("tp","cr")),
              method="GCV.Cp",knots=knots)
    Prediction<-as.vector(predict(mod,NewPoints))
    RMSE<-c(RMSE,sqrt(sum((Prediction-RealValues)^2)/n))
    LABEL<-c(LABEL,"TPS")
    BETA<-c(BETA,mod$coefficients[2])
    
    # Secondo caso
    #   Prodotto tensoriale di sf e cr
    mod <-gam(Data~
                  DesMat+
                  te(xData,yData,TimeData,k=c(NA,length(TimePoints)),d=c(2,1),bs=c("sf","cr"),xt=list(bnd=fsb)),
              method="GCV.Cp",knots=knots)
    Prediction<-as.vector(predict(mod,NewPoints))
    RMSE<-c(RMSE,sqrt(sum((Prediction-RealValues)^2)/n))
    LABEL<-c(LABEL,"SOAP")
    BETA<-c(BETA,mod$coefficients[2])
    
    ## MODELLO SPAZIOTEMPO ##
    # Risolvo
    SolutionObj<-ST.Smooth.Covar(Data,DesMat,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)
    Prediction<-ST.Eval(xValid,yValid,tValid,SolutionObj)
    RMSE<-c(RMSE,sqrt(sum((Prediction-RealValues)^2)/n))
    LABEL<-c(LABEL,"STR-PDE") 
    BETA<-c(BETA,SolutionObj$BetaHat)
    
    ## KRIGING ##
#     STFDF.data <- stConstruct(data.frame(DataK=t(DataMatrix),DesMatK=t(DesMatMat)), 
#                               space = list(values = 1:nrow(DataMatrix)),
#                               time = istanti, 
#                               SpatialObj = punti, 
#                               interval = TRUE)         # object of class "STFDF"
#     # definisco la griglia di predizione pred.grid
#     sequenza_x <- seq(from=min(x), to=max(x), length=120)
#     sequenza_y <- seq(from=min(y), to=max(y), length=40)
#     grd <- SpatialPoints(
#         matrix(data = c(rep(sequenza_x,length(sequenza_y)),rep(sequenza_y,each=length(sequenza_x))),
#                nrow = length(sequenza_x)*length(sequenza_y), ncol = 2, byrow = FALSE)
#     )
#     tgrd <- istanti
#     pred.grid  <- STF(sp=grd, time=tgrd)
#     # interpolate with a separable exponential covariance model
#     vv = variogramST(DataK ~ DesMatK, STFDF.data)  # the spatio-temporal sample variogram
#     Vgm_model <- vgmST("separable", space = vgm(1, "Exp", 1, 1), time = vgm(1, "Exp", 1, 1), sill=1)  # the desired spatio-temporal model (separable exponential covariance model)
#     v <- fit.StVariogram(vv, Vgm_model, method = "L-BFGS-B", lower = c(0,0,0,0,0), upper = c(2,1,6,1,200)) # fit a spatio-temporal variogram of the given type (separable exponential covariance model) to spatio-temporal sample variogram (vv).
#     
#     # Valori stimati col kriging
#     valore_stimato <- krigeST(values ~ DesMat, STFDF.data, pred.grid_C, v)[[1]]
#     
#     # Valori veri della funzione 
#     RealValuesKrig<-fs.test(x_grid,y_grid)*cos(t_grid)
#     
#     
#     # RMSE
#     somma_quadrato_errori <- sum(( RealValuesKrig - valore_stimato)^2,na.rm=T)
#     media_quadrato_errori <- somma_quadrato_errori/sum(!is.na(RealValuesKrig))
#     RMSE1 <- sqrt(media_quadrato_errori)
#     RMSE<-c(RMSE,RMSE1)
#     LABEL<-c(LABEL,"KRIG")
    
}

save(file="Risultati.RData",RMSE,LABEL,BETA)

LABEL<-as.factor(LABEL)
png(filename="Confronto tra i metodi.png")
boxplot(RMSE ~ LABEL, main="Confronto tra metodi", xlab="Metodo",ylab="RMSE")
dev.off()

png(filename="Confronto tra Beta.png")
boxplot(BETA ~ LABEL, main="Confronto tra Beta", xlab="Metodo",ylab="Beta")
dev.off()
