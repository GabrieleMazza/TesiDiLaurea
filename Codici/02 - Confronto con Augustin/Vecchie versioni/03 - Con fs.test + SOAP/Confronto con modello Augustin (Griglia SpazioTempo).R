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
NGEN=1000


##### GCV #####

# ORA GCV, LO APPLICO AI DATI NON PERTURBATI, NON VOGLIO APPLICARLO TUTTE LE VOLTE NEL FOR

TimeBasisObj<-Create.Bspline.Time.Basis(TimePoints,TimeOrder=4,DerivativeOrder=2,PlotIt=F)
SpaceBasisObj<-Create.FEM.Space.Basis(cbind(x,y),Triang,TypePoint,1)

DataMatrix<-NULL
TimeDataMatrix<-NULL
xDataMatrix<-NULL
yDataMatrix<-NULL
for(j in TimePoints)
{
    DataMatrix<-cbind(DataMatrix,Data*fun(j))
    
    #Per dopo..
    TimeDataMatrix<-cbind(TimeDataMatrix,rep(j,length(x[TypePoint])))
    xDataMatrix<-cbind(xDataMatrix,xknot)
    yDataMatrix<-cbind(yDataMatrix,yknot)
}
TimeData<-NULL
xData<-NULL
yData<-NULL
for(i in 1:(dim(TimeDataMatrix)[1]))
{
    TimeData<-c(TimeData,TimeDataMatrix[i,])
    xData<-c(xData,xDataMatrix[i,])
    yData<-c(yData,yDataMatrix[i,])
}
# LogS<-seq(-13,-12,length.out=10)
# LogT<-seq(-15,-14,length.out=10)
# GCVResult<-ST.GCV(DataMatrix,SpaceBasisObj,TimeBasisObj,LogS,LogT)
# png(filename="GCV Matrix.png")
# image(LogS,LogT,GCVResult$GCVMatrix,col=heat.colors(100),main="GCV Matrix",xlab="logLambdaS",ylab="logLambdaT")
# dev.off()
# 
# save(file="GCVResult.RData",GCVResult,LogS,LogT)

# LambdaS=10^GCVResult$Best[1]
# LambdaT=10^GCVResult$Best[2]

LambdaS=10^-12.88889
LambdaT=10^-14.22222



##### CICLO PER IL CONFRONTO DELL'ALGORITMO #####

ErrorGAM<-NULL
ErrorST<-NULL

fsb <- list(fs.boundary())


# Genero i dati per la previsione
xGen<-seq(-1,3.5,length.out=20)
yGen<-seq(-1,1,length.out=10)
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
xValidaz<-NULL
yValidaz<-NULL
for(i in 1:(dim(xx)[1]))
{
    for(j in 1:(dim(xx)[2]))
    {
        x=xx[i,j]
        y=yy[i,j]
        if((pnt.in.poly(cbind(x,y),Bound)$pip==1) & inSide(fsb,x=x,y=y))
        {
            xValidaz<-c(xValidaz,xx[i,j])
            yValidaz<-c(yValidaz,yy[i,j])
        }
    }
}

RealValues<-fs.test(xValidaz,yValidaz)

tGen<-seq(min(TimePoints),max(TimePoints),length.out = 4)

nValid=length(xValidaz)*length(tGen)


fsb=list(list(xbound,ybound))
names(fsb[[1]]) <- c("xData","yData")
knots <- data.frame(xData=xknot,yData=yknot)
fsbReal <- list(fs.boundary())

for(ntimes in 1:NTIMES)
{
    TempST<-NULL
    TempGAM<-NULL
    for(t in tGen)
    {
        #Ricavo i dati
        DataMatrix<-NULL
        for(i in TimePoints)
        {
            if(noise)
            {
                DataMatrix<-cbind(DataMatrix,Data*fun(i)+rnorm(length(xknot),0,0.05))
            } else
            {
                DataMatrix<-cbind(DataMatrix,Data*fun(i))
            }
        }
        DataVec<-NULL
        for(i in 1:(dim(DataMatrix)[1]))
        {
            DataVec<-c(DataVec,DataMatrix[i,])
        }
        
        ## MODELLO AUGUSTIN ##
        mod <-gam(DataVec~
                      te(xData,yData,TimeData,d=c(2,1),bs=c("tp","cr"))+
                      s(xData,yData,bs="so",xt=list(bnd=fsb)),
                  method="GCV.Cp",knots=knots)
        
        summary(mod)
        
        ## MODELLO SPAZIOTEMPO ##
        # Risolvo
        SolutionObj<-ST.Smooth(DataVec,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)
        
        
        ## PREVISIONE ##
        
        tvec=rep(t,length(xValidaz))
        
        # GAM
        NewPoints<-data.frame(xData=xValidaz,yData=yValidaz,TimeData=tvec)
        PredictionGAM<-as.vector(predict(mod,NewPoints))
        # ST
        PredictionST<-ST.Eval(xValidaz,yValidaz,tvec,SolutionObj)

        TempST<-c(TempST,sum((PredictionST-RealValues*cos(t))^2))
        TempGAM<-c(TempGAM,sum((PredictionGAM-RealValues*cos(t))^2))
        
    }

    # Studio degli errori
    ErrorGAM<-c(ErrorGAM,sum(TempGAM)/nValid)
    ErrorST<-c(ErrorST,sum(TempST)/nValid)
    
}

png(filename="Confronto tra i metodi.png")
plot(ErrorGAM,ylim=c(min(ErrorGAM,ErrorST),max(ErrorGAM,ErrorST)),main="Confronto tra metodi in griglia spazio-tempo",type='l')
points(ErrorST,type='l',col="red")
legend("bottomright",col=c("black","red"),legend=c("GAM","ST"),lty=1)
dev.off()

save(file="Risultati.RData",ErrorGAM,ErrorST)

