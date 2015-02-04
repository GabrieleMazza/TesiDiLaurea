library(rgl)
library(mgcv)
library(fda)
source("SpazioTempo.R")

# Devo aggiungere del rumore al dato?
noise<-TRUE

# Per quanti istanti di tempo eseguire la stima?
TimePoints<-0:5

# Lambda
LambdaS=10^-3
LambdaT=10^-5
logS=-3
logT=-5

# Function for perturbation
fun=function(x)
{
    y=cos(x)
    return(y)
}


##### PERTURBAZIONE TEMPORALE DEI DATI #####

load("FrontieraC.RData")
DataMatrix<-NULL
# Perturbo secondo la funzione cos(x)
for(i in TimePoints)
{
    if(noise)
    {
        DataMatrix<-cbind(DataMatrix,(Data*fun(i)+rnorm(length(Data),0,0.005)))
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
load("Validazione.RData")

# RICERCA DEL MASSIMO
#PER AVERE UNA SCALA DI COLORI PER IL GRAFICO UNIFORME IN TUTTI GLI ISTANTI DI TEMPO
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
TimeBasisObj<-Create.Bspline.Time.Basis(TimePoints,4,PlotIt=T)
SpaceBasisObj<-Create.FEM.Space.Basis(cbind(x,y),Triang,TypePoint,1)

# for(logS in -2:-1)
# {
#     for(logT in -5:1)
#     {
#         LambdaS=10^logS
#         LambdaT=10^logT
        #Introduco anche la misurazione del tempo per la funzione chiave
        SolutionObj<-smooth.ST.fd(DataMatrix,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT,2)
        
        # Ora salvo i risultati
        write.table(SolutionObj$C,file=paste(logS, logT,"MatriceC.txt",sep=" "),row.names=FALSE,col.names=FALSE)
                
        ##### GRAFICI #####
        
        PlotMatrix=matrix(ncol=length(yvec),nrow=length(xvec))
           
        # Voglio ricreare perfettamente le posizioni della image matrix tru

        # RICERCA DEL MASSIMO
        #PER AVERE UNA SCALA DI COLORI PER IL GRAFICO UNIFORME IN TUTTI GLI ISTANTI DI TEMPO
    
        Max=MaxD
        Min=MinD
        
        for(j in TimePoints)
        {
            Time<-rep(j,length(xValid))
            Result<-eval.ST.fd(xValid,yValid,Time,SolutionObj)
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
        

        # GRAFICI

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
            Result<-eval.ST.fd(xValid,yValid,Time,SolutionObj)
            
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
        
        # Plot in un punto fisso
        png(filename="Plot per un punto fissato.png")
        FixedPointPlot(3,-0.5,SolutionObj)
        dev.off()
#     }
# }

