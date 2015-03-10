load("FrontieraC.RData")
load("Validazione.RData")

source("SpazioTempo.R")
library(mgcv)
library(SDMTools)
library(MASS)
library(ggplot2)

# Devo aggiungere del rumore al dato?
noise<-TRUE

# DEVIAZIONE STANDARD DELL'ERRORE
SigmaNoise=0.1

# Per quanti istanti di tempo eseguire la stima?
TimePoints<-0:5

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

# Quante volte eseguire il ciclo?
NTIMES=10



##### OGGETTI UTILI #####

TimeBasisObj<-Create.Bspline.Time.Basis(TimePoints,TimeOrder=4,DerivativeOrder=2,PlotIt=F)
SpaceBasisObj<-Create.FEM.Space.Basis(cbind(x,y),Triang,TypePoint,1)

DataMatrix<-NULL
TimeDataMatrix<-NULL
xDataMatrix<-NULL
yDataMatrix<-NULL
for(t in TimePoints)
{
    DataMatrix<-cbind(DataMatrix,fun(xknot,yknot,rep(t,length(xknot))))
    
    #Per dopo..
    TimeDataMatrix<-cbind(TimeDataMatrix,rep(t,length(x[TypePoint])))
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

# LogS<--8:0
# LogT<--8:0
# GCVResult<-ST.GCV(Data,SpaceBasisObj,TimeBasisObj,LogS,LogT)
# png(filename="GCV Matrix.png")
# image(LogS,LogT,GCVResult$GCVMatrix,col=heat.colors(100),main="GCV Matrix",xlab="logLambdaS",ylab="logLambdaT")
# dev.off()
# 
# save(file="GCVResult.RData",GCVResult,LogS,LogT)
# 
# LambdaS=10^GCVResult$Best[1]
# LambdaT=10^GCVResult$Best[2]

LambdaS=10^-8
LambdaT=10^-8



##### PUNTI DI VALIDAZIONE, PER LA GRIGLIA SPAZIOTEMPO #####

fsb <- list(fs.boundary())

# Genero i dati per la previsione
xGen<-seq(-1,3.5,length.out=20)
yGen<-seq(-1,1,length.out=10)
tGen<-seq(min(TimePoints),max(TimePoints),length.out = 4)

# Trasformo i dati in x e y in una griglia...
xx<-NULL
for(i in 1:length(yGen))
{
    xx<-cbind(xx,xGen)
}
yy<-NULL
for(j in 1:length(xGen))
{
    yy<-rbind(yy,yGen)
}
xValid<-NULL
yValid<-NULL
for(i in 1:(dim(xx)[1]))
{
    for(j in 1:(dim(xx)[2]))
    {
        x=xx[i,j]
        y=yy[i,j]
        if((pnt.in.poly(cbind(x,y),Bound)$pip==1) & inSide(fsb,x=x,y=y))
        {
            xValid<-c(xValid,xx[i,j])
            yValid<-c(yValid,yy[i,j])
        }
    }
}

tValid<-NULL
for (i in tGen)
{
    tValid<-c(tValid,rep(i,length(xValid)))    
}
xValid<-rep(xValid,length(tGen))
yValid<-rep(yValid,length(tGen))

#Valori reali sulla griglia
RealValues<-fun(xValid,yValid,tValid)




##### ESEMPIO DI BUBBLEPLOT DEI DATI #####

# Genero una volta i dati allo stesso modo di come farò dopo nel ciclo,
# per vederne i bubbleplot...

# Ho bisogno dei un dataframe
Knot=data.frame(xknot=c(xknot,10,10),yknot=c(yknot,1,-1))
Boundary<-data.frame(xbound,ybound)
# Genero i dati
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

for(j in 1: length(TimePoints))
{
    Gen<-c(DataMatrix[,j],min(DataMatrix),max(DataMatrix))
    
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

rm(Knot,Boundary,DataMatrix)




##### CICLO DEI MODELLI #####

RMSE<-NULL
LABEL<-NULL

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
            DataMatrix<-cbind(DataMatrix,fun(xknot,yknot,j))
        }
    }
    Data<-NULL
    for(i in 1:(dim(DataMatrix)[1]))
    {
        Data<-c(Data,DataMatrix[i,])
    }
    
    NewPoints<-data.frame(xData=xValid,yData=yValid,TimeData=tValid)
    
    ## MODELLI TIPO AUGUSTIN ##
    
    # Primo caso
    #   Prodotto tensoriale di tps e cr
    mod <-gam(Data~
                  te(xData,yData,TimeData,k=c(NA,length(TimePoints)),d=c(2,1),bs=c("tp","cr")),
                  method="REML",knots=knots)
    Prediction<-as.vector(predict(mod,NewPoints))
    RMSE<-c(RMSE,sqrt(sum((Prediction-RealValues)^2)/n))
    LABEL<-c(LABEL,"AUG1")
    
    # Secondo caso
    #   Prodotto tensoriale di sf e cr
    mod <-gam(Data~
                  te(xData,yData,TimeData,k=c(NA,length(TimePoints)),d=c(2,1),bs=c("sf","cr"),xt=list(bnd=fsb)),
              method="REML",knots=knots)
    Prediction<-as.vector(predict(mod,NewPoints))
    RMSE<-c(RMSE,sqrt(sum((Prediction-RealValues)^2)/n))
    LABEL<-c(LABEL,"AUG2")
    
    # Terzo caso
    #   Marginale tps
    #   Marginale cr
    mod <-gam(Data~
                  s(xData,yData,bs="tp")+
                  s(TimeData,k=length(TimePoints),bs=c("cr")),
              method="REML",knots=knots)
    Prediction<-as.vector(predict(mod,NewPoints))
    RMSE<-c(RMSE,sqrt(sum((Prediction-RealValues)^2)/n))
    LABEL<-c(LABEL,"AUG3")
    
    # Quarto caso
    #   Marginale sf
    #   Marginale cr
    mod <-gam(Data~
                  s(xData,yData,bs="so",xt=list(bnd=fsb))+
                  s(TimeData,k=length(TimePoints),bs=c("cr")),
              method="REML",knots=knots)
    Prediction<-as.vector(predict(mod,NewPoints))
    RMSE<-c(RMSE,sqrt(sum((Prediction-RealValues)^2)/n))
    LABEL<-c(LABEL,"AUG4")
    
    # Quinto caso
    #   Prodotto tensoriale di tps e cr
    #   Marginale tps
    #   Marginale cr
    mod <-gam(Data~
                  te(xData,yData,TimeData,k=c(NA,length(TimePoints)),d=c(2,1),bs=c("tp","cr"))+
                  s(xData,yData,bs="tp")+
                  s(TimeData,k=length(TimePoints),bs=c("cr")),
              method="REML",knots=knots)
    Prediction<-as.vector(predict(mod,NewPoints))
    RMSE<-c(RMSE,sqrt(sum((Prediction-RealValues)^2)/n))
    LABEL<-c(LABEL,"AUG5")
    
    # Sesto caso
    #   Prodotto tensoriale di sf e cr
    #   Marginale sf
    #   Marginale cr
    mod <-gam(Data~
                  te(xData,yData,TimeData,k=c(NA,length(TimePoints)),d=c(2,1),bs=c("sf","cr"),xt=list(bnd=fsb))+
                  s(xData,yData,bs="so",xt=list(bnd=fsb))+
                  s(TimeData,k=length(TimePoints),bs=c("cr")),
              method="REML",knots=knots)
    Prediction<-as.vector(predict(mod,NewPoints))
    RMSE<-c(RMSE,sqrt(sum((Prediction-RealValues)^2)/n))
    LABEL<-c(LABEL,"AUG6")
    
    # Settimo caso
    #   Prodotto tensoriale di sw e cr
    #     mod <-gam(Data~
    #                   te(xData,yData,TimeData,k=c(NA,length(TimePoints)),d=c(2,1),bs=c("sw","cr"),xt=list(bnd=fsb)),
    #               method="REML",knots=knots)
    #     Prediction<-as.vector(predict(mod,NewPoints))
    #     RMSE<-c(RMSE,sqrt(sum((Prediction-RealValues)^2)/n))
    #     LABEL<-c(LABEL,"AUG7")

    
    ## MODELLO SPAZIOTEMPO ##
    # Risolvo
    SolutionObj<-ST.Smooth(Data,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)
    Prediction<-ST.Eval(xValid,yValid,tValid,SolutionObj)
    RMSE<-c(RMSE,sqrt(sum((Prediction-RealValues)^2)/n))
    LABEL<-c(LABEL,"ST") 
        
    
}

save(file="Risultati.RData",RMSE,LABEL)

LABEL<-as.factor(LABEL)
png(filename="Confronto tra i metodi.png")
boxplot(RMSE ~ LABEL, main="Confronto tra metodi", xlab="Metodo",ylab="RMSE")
dev.off()



# Eventuale zoom


RMSE2<-RMSE[(LABEL!="AUG3" & LABEL!="AUG4" & LABEL!="AUG5" & LABEL!="AUG6")]
LABEL2<-LABEL[(LABEL!="AUG3" & LABEL!="AUG4" & LABEL!="AUG5" & LABEL!="AUG6")]

png(filename="Zoom.png")
boxplot(RMSE2 ~ LABEL2, main="Confronto tra metodi zoom", xlab="Metodo",ylab="RMSE")
dev.off()
