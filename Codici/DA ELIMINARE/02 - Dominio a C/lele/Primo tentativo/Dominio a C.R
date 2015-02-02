library(rgl)
library(mgcv)
library(fda)
source("SpazioTempo.R")

# Devo aggiungere del rumore al dato?
noise<-TRUE

# Per quanti istanti di tempo eseguire la stima?
TimePoints<-0:5

# Lambda
LambdaS=10^-6
LambdaT=10^-6



##### PERTURBAZIONE TEMPORALE DEI DATI #####


fun=function(x)
{
    y=(-1/6)*x*x+5/6*x
    return(y)
}


load("FrontieraC.RData")
DataMatrix<-NULL
# Perturbo secondo la funzione cos(x)
for(i in TimePoints)
{
    if(noise)
    {
        DataMatrix<-cbind(DataMatrix,(Data*fun(i)+rnorm(length(Data),0,0.05)))
    }
    else
    {
        DataMatrix<-cbind(DataMatrix,Data*fun(i))
    }
}
# Per i grafici che devo fare dopo, mi salvo il massimo e il minimo di questa matrice
MaxD=max(DataMatrix)
MinD=min(DataMatrix)
#Cerco il massimo
nx<-length(xvec)
ny<-length(yvec)
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))
true <- matrix(fs.test(xx,yy),nx,ny)
for(j in TimePoints)
{
    if(MaxD<max(true*fun(j),na.rm=TRUE))
    {
        MaxD=max(true*fun(j),na.rm=TRUE)
    }
    if(MinD>min(true*fun(j),na.rm=TRUE))
    {
        MinD=min(true*fun(j),na.rm=TRUE)
    }
}
##### RISOLVO IL SISTEMA #####

#Creo le basi in spazio e tempo
TimeBasisObj<-Create.Bspline.Time.Basis(TimePoints,3,T)
SpaceBasisObj<-Create.FEM.Space.Basis(cbind(x,y),Triang,TypePoint,1)

for(logS in -1:1)
{
    for(logT in -8:-2)
    {
        LambdaS=10^logS
        LambdaT=10^logT
        #Introduco anche la misurazione del tempo per la funzione chiave
        t1<-Sys.time()
        C<-smooth.ST.fd(DataMatrix,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)
        t2<-Sys.time()
        DeltaT<-t2-t1
        
        # Ora salvo i risultati
        #write.table(C,file="C.txt",row.names=FALSE,col.names=FALSE)
        #save(file="Risultati.RData",x,y,TypePoint,TimePoints,SpaceBasisObj,TimeBasisObj)
        
        
        ##### GRAFICI #####
        
        load("Validazione.RData")
        
        PlotMatrix=matrix(ncol=length(yvec),nrow=length(xvec))
        
        
        Max=MaxD
        Min=MinD
        # Voglio ricreare perfettamente le posizioni della image matrix tru
        for(j in TimePoints)
        {
            Time<-rep(j,length(xValid))
            Result<-eval.ST.fd(xValid,yValid,Time,C,SpaceBasisObj,TimeBasisObj)
            if(Max<max(Result,na.rm=TRUE))
            {
                Max=max(Result,na.rm=TRUE)
            }
            if(Min>min(Result,na.rm=TRUE))
            {
                Min=min(Result,na.rm=TRUE)
            }
        }
        
        zlim<-c(Min,Max)
        
        # Ora salvo tutti i grafici
        # Voglio ricreare perfettamente le posizioni della image matrix tru
        
        # Prima ricavo il grafico della funzione reale
        nx<-length(xvec)
        ny<-length(yvec)
        xx <- rep(xvec,ny)
        yy<-rep(yvec,rep(nx,ny))
        true <- matrix(fs.test(xx,yy),nx,ny)
        for(j in TimePoints)
        {
            png(filename=paste(logS, logT,j,"reale.png",sep=" "))
            image(xvec,yvec,true*fun(j),zlim=zlim,main=paste("Funzione reale",j,logS,logT,sep=" "))
            lines(x[TypePoint==FALSE],y[TypePoint==FALSE],lwd=3)
            contour(xvec,yvec,true*fun(j),levels=seq(-5,5,by=.25),add=TRUE)
            dev.off()
        }
        
        # Ora le stime
        for(j in TimePoints)
        {
            Time<-rep(j,length(xValid))
            Result<-eval.ST.fd(xValid,yValid,Time,C,SpaceBasisObj,TimeBasisObj)
            
            # Costruisco ora la matrice con l'image plot
            for(i in 1:(dim(PosMatrix)[1]))
            {
                PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=Result[i]
            }
            
            png(filename=paste(logS, logT,j,"stimata.png",sep=" "))
            image(xvec,yvec,PlotMatrix,zlim=zlim,main=paste("Funzione stimata",j,logS,logT,sep=" "))
            lines(x[TypePoint==FALSE],y[TypePoint==FALSE],lwd=3)
            contour(xvec,yvec,PlotMatrix,levels=seq(-5,5,by=.25),add=TRUE)
            dev.off()
        }
        
        print(paste("Execution of smooth.ST.fd: ",DeltaT," sec",sep=""))
    }
}

