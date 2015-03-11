rm(list=ls())
graphics.off()

set.seed(1)

###########################################

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


# Devo aggiungere del rumore al dato?
noise<-TRUE

# DEVIAZIONE STANDARD DELL'ERRORE
SigmaNoise=0.5

# Per quanti istanti di tempo eseguire la stima?
TimePoints<-seq(0,2*pi,length.out=9)

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
NX=80
NY=40
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



##### GCV #####

# LogS<-seq(-1,+1,by=0.25)
# LogT<-seq(-4,-2,by=0.25)
# GCVResult<-ST.GCV(Data,SpaceBasisObj,TimeBasisObj,LogS,LogT)
# png(filename="GCV Matrix.png")
# image(LogS,LogT,GCVResult$GCVMatrix,col=heat.colors(100),main="GCV Matrix",xlab="logLambdaS",ylab="logLambdaT")
# dev.off()
# 
# save(file="GCVResult.RData",GCVResult,LogS,LogT)
# 
# LambdaS=10^GCVResult$Best[1]
# LambdaT=10^GCVResult$Best[2]

LambdaS=10^-0.25
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

t_grid<-rep(seq(from=min(TimePoints),to=max(TimePoints),length.out=NT),each=NY*NX)
n_t_grid <- length(t_grid)
x_grid<-rep(rep(seq(from=-1,to=3.5,length=NX),times=NY),NT)
y_grid<-rep(rep(seq(from=-1,to=1,length=NY),each=NX),NT)
grd <- SpatialPoints(
    matrix(data = c(rep(seq(from=-1,to=3.5,length=NX),times=NY),rep(seq(from=-1,to=1,length=NY),each=NX)),
           nrow = NX*NY, ncol = 2, byrow = FALSE)
)
tgrd <- as.POSIXct("2010-08-05")+3600*seq(from=min(TimePoints),to=max(TimePoints),length.out=NT)
pred.grid_C <- STF(sp=grd, time=tgrd)


# ##### ESEMPIO DI BUBBLEPLOT DEI DATI #####
# 
# # Genero una volta i dati allo stesso modo di come farò dopo nel ciclo,
# # per vederne i bubbleplot...
# 
# # Ho bisogno dei un dataframe
# Knot=data.frame(xknot=c(xknot,10,10),yknot=c(yknot,1,-1))
# Boundary<-data.frame(xbound,ybound)
# # Genero i dati
# DataMatrix<-NULL
# for(t in TimePoints)
# {
#     if(noise)
#     {
#         DataMatrix<-cbind(DataMatrix,fun(xknot,yknot,rep(t,length(xknot)))+rnorm(length(xknot),0,SigmaNoise))
#     } else
#     {
#         DataMatrix<-cbind(DataMatrix,fun(xknot,yknot,rep(t,length(xknot))))
#     }
# }
# 
# for(j in 1:length(TimePoints))
# {
#     Gen<-c(DataMatrix[,j],min(DataMatrix),max(DataMatrix))
#     
#     ggplot(data=Knot, aes(x=xknot, y=yknot)) +
#         geom_point(aes(size=(Gen),color=Gen)) +
#         scale_color_gradient(low="yellow", high="red") +
#         geom_polygon(data=Boundary,aes(x=xbound,y=ybound),alpha=1,colour="black", fill=NA, size=1.1) +
#         theme_bw() +
#         labs(color="Dato")+
#         scale_x_continuous(name="x", limits=c(-1,3.5)) +
#         scale_y_continuous(name="y",limits=c(-1,1)) +
#         labs(title=paste("Dati Generati al tempo ",round(TimePoints[j],3),sep="")) + 
#         theme(plot.title = element_text(size = rel(1.5), face="bold")) +
#         guides(size=FALSE)
#     ggsave(file=paste("Dati_tempo",j,".png",sep=""))
# }
# 
# rm(Knot,Boundary,DataMatrix)
# 
# dev.off()



# ##### CICLO DEI MODELLI #####
# 
# RMSE<-NULL
# LABEL<-NULL
# 
# #Ricavo la frontiera da usare nel caso in cui ci sia Soap Film Smoothing nel GAM
# fsb=list(list(xbound,ybound))
# names(fsb[[1]]) <- c("xData","yData")
# knots <- data.frame(xData=xknot,yData=yknot)
# 
# n=length(RealValues)
# 
# # Creo la matrice di rumore, per la replicabilità del dato
# NoiseMatrix<-NULL
# for(ntimes in 1:NTIMES)
# {
#     NoiseMatrix<-cbind(NoiseMatrix,rnorm(length(xknot),0,SigmaNoise))
# }
# 
# # CICLO
# for(ntimes in 1:NTIMES)
# {
#     print(paste("Simulazione ",ntimes,"/",NTIMES,sep=""))
#     #Ricavo i dati
#     DataMatrix<-NULL
#     for(t in TimePoints)
#     {
#         if(noise)
#         {
#             DataMatrix<-cbind(DataMatrix,fun(xknot,yknot,rep(t,length(xknot)))+NoiseMatrix[,ntimes])
#         } else
#         {
#             DataMatrix<-cbind(DataMatrix,fun(xknot,yknot,rep(t,length(xknot))))
#         }
#     }
#     
#     # Per tutti, eccetto il kriging
#     Data<-NULL
#     for(i in 1:(dim(DataMatrix)[1]))
#     {
#         Data<-c(Data,DataMatrix[i,])
#     }
#     NewPoints<-data.frame(xData=xValid,yData=yValid,TimeData=tValid)
#         
#     ## MODELLI TIPO AUGUSTIN ##
#     
#     # Primo caso
#     #   Prodotto tensoriale di tps e cr
#     mod <-gam(Data~
#                   te(xData,yData,TimeData,k=c(NA,length(TimePoints)),d=c(2,1),bs=c("tp","cr")),
#                   method="REML",knots=knots)
#     Prediction<-as.vector(predict(mod,NewPoints))
#     RMSE<-c(RMSE,sqrt(sum((Prediction-RealValues)^2)/n))
#     LABEL<-c(LABEL,"TPS")
#     
#     # Secondo caso
#     #   Prodotto tensoriale di sf e cr
#     mod <-gam(Data~
#                   te(xData,yData,TimeData,k=c(NA,length(TimePoints)),d=c(2,1),bs=c("sf","cr"),xt=list(bnd=fsb)),
#               method="REML",knots=knots)
#     Prediction<-as.vector(predict(mod,NewPoints))
#     RMSE<-c(RMSE,sqrt(sum((Prediction-RealValues)^2)/n))
#     LABEL<-c(LABEL,"SOAP")
# 
#     
#     ## MODELLO SPAZIOTEMPO ##
#     # Risolvo
#     SolutionObj<-ST.Smooth(Data,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)
#     Prediction<-ST.Eval(xValid,yValid,tValid,SolutionObj)
#     RMSE<-c(RMSE,sqrt(sum((Prediction-RealValues)^2)/n))
#     LABEL<-c(LABEL,"STR-PDE")
#     
#     
#     ## KRIGING ##
#     STFDF.data <- stConstruct(t(DataMatrix), 
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
#     vv = variogramST(values ~ 1, STFDF.data)  # the spatio-temporal sample variogram
#     Vgm_model <- vgmST("separable", space = vgm(1, "Exp", 1, 1), time = vgm(1, "Exp", 1, 1), sill=1)  # the desired spatio-temporal model (separable exponential covariance model)
#     v <- fit.StVariogram(vv, Vgm_model, method = "L-BFGS-B", lower = c(0,0,0,0,0), upper = c(2,1,6,1,200)) # fit a spatio-temporal variogram of the given type (separable exponential covariance model) to spatio-temporal sample variogram (vv).
#     
#     # Valori stimati col kriging
#     valore_stimato <- krigeST(values ~ 1, STFDF.data, pred.grid_C, v)[[1]]
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
#     
# }
# 
# save(file="Risultati.RData",RMSE,LABEL)
# 
# 
# LABEL[which(LABEL=="STR-PDE")]<-"STSR"
# LABEL<-as.factor(LABEL)
# png(filename="Confronto_metodi.png")
# boxplot(RMSE ~ LABEL, main="Confronto tra metodi", xlab="Metodo",ylab="RMSE")
# dev.off()





##### PLOT PER I CONFRONTI DEI METODI #####

# Creo dei nuovi dati
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

# Per tutti, eccetto il kriging
Data<-NULL
for(i in 1:(dim(DataMatrix)[1]))
{
    Data<-c(Data,DataMatrix[i,])
}

fsb=list(list(xbound,ybound))
names(fsb[[1]]) <- c("xData","yData")
knots <- data.frame(xData=xknot,yData=yknot)


#PRIMO CICLO: CALCOLA TUTTI I VALORI
PredictionTPSVect<-NULL
PredictionSOAPVect<-NULL
PredictionKRIGVect<-NULL
PredictionSTSRVect<-NULL
RealValuesVect<-NULL
NXplot=300
NYplot=100
#Voglio una griglia fine
xGen=seq(from=-1,to=3.5,length.out=NXplot)
yGen=seq(from=-1,to=1,length.out=NYplot)
xPlot<-NULL
yPlot<-NULL
PosMatrix<-NULL
for(i in 1:NXplot)
{
    for(j in 1:NYplot)
    {
        if(!is.na(fs.test(xGen[i],yGen[j])))
        {
            xPlot<-c(xPlot,xGen[i])
            yPlot<-c(yPlot,yGen[j])
            PosMatrix<-rbind(PosMatrix,c(i,j))
        }
    }
}

for(jj in 1:length(TimePoints))
{
    TimeSelected<-TimePoints[jj]
    tPlot<-NULL
    RealValues<-NULL
    
    for(i in 1:NXplot)
    {
        for(j in 1:NYplot)
        {
            if(!is.na(fs.test(xGen[i],yGen[j])))
            {
                tPlot<-c(tPlot,TimeSelected)
                RealValues<-c(RealValues,fun(xGen[i],yGen[j],TimeSelected))
            }
        }
    }
    
    NewPoints<-data.frame(xData=xPlot,yData=yPlot,TimeData=tPlot)
    
    t_grid<-rep(TimeSelected,each=NYplot*NXplot)
    x_grid<-rep(seq(from=-1,to=3.5,length=NXplot),times=NYplot)
    y_grid<-rep(seq(from=-1,to=1,length=NYplot),each=NXplot)
    grd <- SpatialPoints(
        matrix(data = c(x_grid,y_grid),
               nrow = NXplot*NYplot, ncol = 2, byrow = FALSE)
    )
    tgrd <- as.POSIXct("2010-08-05")+3600*c(TimeSelected,7)
    pred.grid_C <- STF(sp=grd, time=tgrd)
    
    # Primo caso
    #   Prodotto tensoriale di tps e cr
    mod <-gam(Data~
                  te(xData,yData,TimeData,k=c(NA,length(TimePoints)),d=c(2,1),bs=c("tp","cr")),
              method="REML",knots=knots)
    PredictionTPS<-as.vector(predict(mod,NewPoints))
    
    # Secondo caso
    #   Prodotto tensoriale di sf e cr
    mod <-gam(Data~
                  te(xData,yData,TimeData,k=c(NA,length(TimePoints)),d=c(2,1),bs=c("sf","cr"),xt=list(bnd=fsb)),
              method="REML",knots=knots)
    PredictionSOAP<-as.vector(predict(mod,NewPoints))
    
    
    ## MODELLO SPAZIOTEMPO ##
    # Risolvo
    SolutionObj<-ST.Smooth(Data,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)
    PredictionSTSR<-ST.Eval(xPlot,yPlot,tPlot,SolutionObj)
    
    ## KRIGING ##
    STFDF.data <- stConstruct(t(DataMatrix), 
                              space = list(values = 1:nrow(DataMatrix)),
                              time = istanti, 
                              SpatialObj = punti, 
                              interval = TRUE)         # object of class "STFDF"
    # definisco la griglia di predizione pred.grid
    sequenza_x <- seq(from=min(x), to=max(x), length=120)
    sequenza_y <- seq(from=min(y), to=max(y), length=40)
    grd <- SpatialPoints(
        matrix(data = c(rep(sequenza_x,length(sequenza_y)),rep(sequenza_y,each=length(sequenza_x))),
               nrow = length(sequenza_x)*length(sequenza_y), ncol = 2, byrow = FALSE)
    )
    tgrd <- istanti
    pred.grid  <- STF(sp=grd, time=tgrd)
    # interpolate with a separable exponential covariance model
    vv = variogramST(values ~ 1, STFDF.data)  # the spatio-temporal sample variogram
    Vgm_model <- vgmST("separable", space = vgm(1, "Exp", 1, 1), time = vgm(1, "Exp", 1, 1), sill=1)  # the desired spatio-temporal model (separable exponential covariance model)
    v <- fit.StVariogram(vv, Vgm_model, method = "L-BFGS-B", lower = c(0,0,0,0,0), upper = c(2,1,6,1,200)) # fit a spatio-temporal variogram of the given type (separable exponential covariance model) to spatio-temporal sample variogram (vv).
    # Valori stimati col kriging
    PredictionKRIG <- (krigeST(values ~ 1, STFDF.data, pred.grid_C, v)[[1]])[1:(NXplot*NYplot)]

    PredictionTPSVect<-cbind(PredictionTPSVect,PredictionTPS)
    PredictionSOAPVect<-cbind(PredictionSOAPVect,PredictionSOAP)
    PredictionKRIGVect<-cbind(PredictionKRIGVect,PredictionKRIG)
    PredictionSTSRVect<-cbind(PredictionSTSRVect,PredictionSTSR)
    RealValuesVect<-cbind(RealValuesVect,RealValues)
}

    
Max<-max(DataMatrix,PredictionTPSVect,PredictionSOAPVect,PredictionKRIGVect,PredictionSTSRVect,RealValuesVect,na.rm=FALSE)
Min<-min(DataMatrix,PredictionTPSVect,PredictionSOAPVect,PredictionKRIGVect,PredictionSTSRVect,RealValuesVect,na.rm=FALSE)
## PLOT ##


## PLOT DATI GENERATI
Knot=data.frame(xknot=c(xknot,10,10),yknot=c(yknot,1,-1))
Boundary<-data.frame(xbound,ybound)

for(j in 1:length(TimePoints))
{
    Gen<-c(DataMatrix[,j],Min,Max)
    
    ggplot(data=Knot, aes(x=xknot, y=yknot)) +
        geom_point(aes(size=(Gen),color=Gen)) +
        scale_color_gradient(low="yellow", high="red") +
        geom_polygon(data=Boundary,aes(x=xbound,y=ybound),alpha=1,colour="black", fill=NA, size=1.1) +
        theme_bw() +
        labs(color="Dato")+
        scale_x_continuous(name="x", limits=c(-1,3.5)) +
        scale_y_continuous(name="y",limits=c(-1,1)) +
        labs(title=paste("Dati Generati al tempo ",round(TimePoints[j],3),sep="")) + 
        theme(plot.title = element_text(size = rel(1.5), face="bold")) +
        guides(size=FALSE)
    ggsave(file=paste("Dati_tempo",j,".png",sep=""),height=5,width=6)
}

rm(Knot,Boundary)

dev.off()


## PLOT FUNZIONI STIMATE E REALI
PlotMatrix=matrix(nrow=NXplot,ncol=(NYplot+1))
#Agiungo una colonna che sarà sfruttata per inserire i massimi e i minimi
# e avere una scala di colori uniforme

PlotMatrix[1,(NYplot+1)]=Max
PlotMatrix[2,(NYplot+1)]=Min

#Ora Taglierò questa colonna con ylab
for(jj in 1:length(TimePoints))
{
    
    TimeSelected<-TimePoints[jj]
    # TPS
    for(i in 1:(dim(PosMatrix)[1]))
    {
        PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=PredictionTPSVect[i,jj]
    }
    # Plot
    png(filename=paste("TPStempo",jj,".png",sep=""))
    image(xGen,c(yGen,2),PlotMatrix,col=heat.colors(100),main=paste("Estimated function at time ",round(TimeSelected,3)," (TPS)",sep=""),xlab="",ylab="",ylim=c(-1,1))
    lines(xbound,ybound,lwd=3)
    lines(c(xbound[1],xbound[length(xbound)]),c(ybound[1],ybound[length(ybound)]),lwd=3)
    contour(xGen,c(yGen,2),PlotMatrix,nlevels=10,add=TRUE,lwd=2,labcex=1.1)
    dev.off()
    
    # SOAP
    for(i in 1:(dim(PosMatrix)[1]))
    {
        PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=PredictionSOAPVect[i,jj]
    }
    # Plot
    png(filename=paste("SOAPtempo",jj,".png",sep=""))
    image(xGen,c(yGen,2),PlotMatrix,col=heat.colors(100),main=paste("Estimated function at time ",round(TimeSelected,3)," (SOAP)",sep=""),xlab="",ylab="",ylim=c(-1,1))
    lines(xbound,ybound,lwd=3)
    lines(c(xbound[1],xbound[length(xbound)]),c(ybound[1],ybound[length(ybound)]),lwd=3)
    contour(xGen,c(yGen,2),PlotMatrix,nlevels=10,add=TRUE,lwd=2,labcex=1.1)
    dev.off()
    
    # ST
    for(i in 1:(dim(PosMatrix)[1]))
    {
        PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=PredictionSTSRVect[i,jj]
    }
    # Plot
    png(filename=paste("STSRtempo",jj,".png",sep=""))
    image(xGen,c(yGen,2),PlotMatrix,col=heat.colors(100),main=paste("Estimated function at time ",round(TimeSelected,3)," (STSR)",sep=""),xlab="",ylab="",ylim=c(-1,1))
    lines(xbound,ybound,lwd=3)
    lines(c(xbound[1],xbound[length(xbound)]),c(ybound[1],ybound[length(ybound)]),lwd=3)
    contour(xGen,c(yGen,2),PlotMatrix,nlevels=10,add=TRUE,lwd=2,labcex=1.1)
    dev.off()
    
    # REALE
    for(i in 1:(dim(PosMatrix)[1]))
    {
        PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=RealValuesVect[i,jj]
    }
    # Plot
    png(filename=paste("REALEtempo",jj,".png",sep=""))
    image(xGen,c(yGen,2),PlotMatrix,col=heat.colors(100),main=paste("Real function at time ",round(TimeSelected,3),sep=""),xlab="",ylab="",ylim=c(-1,1))
    lines(xbound,ybound,lwd=3)
    lines(c(xbound[1],xbound[length(xbound)]),c(ybound[1],ybound[length(ybound)]),lwd=3)
    contour(xGen,c(yGen,2),PlotMatrix,nlevels=10,add=TRUE,lwd=2,labcex=1.1)
    dev.off()
    
    # KRIG
    for(i in 1:(dim(PosMatrix)[1]))
    {
        ind<-(PosMatrix[i,2]-1)*NXplot+PosMatrix[i,1]
        PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=PredictionKRIGVect[ind,jj]
    }
    # Plot
    png(filename=paste("KRIGtempo",jj,".png",sep=""))
    image(xGen,c(yGen,2),PlotMatrix,col=heat.colors(100),main=paste("Estimated function at time ",round(TimeSelected,3)," (KRIG)",sep=""),xlab="",ylab="",ylim=c(-1,1))
    lines(xbound,ybound,lwd=3)
    lines(c(xbound[1],xbound[length(xbound)]),c(ybound[1],ybound[length(ybound)]),lwd=3)
    contour(xGen,c(yGen,2),PlotMatrix,nlevels=10,add=TRUE,lwd=2,labcex=1.1)
    dev.off()
    
}