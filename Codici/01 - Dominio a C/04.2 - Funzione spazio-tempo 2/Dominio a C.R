library(rgl)
library(mgcv)
library(SDMTools)
library(fda)
source("SpazioTempo.R")
load("FrontieraC.RData")
load("Validazione.RData")

# Devo aggiungere del rumore al dato?
noise<-TRUE

# Per quanti istanti di tempo eseguire la stima?
TimePoints<-0:5

# Lambda
# LambdaS=10^-3
# LambdaT=10^-5
# logS=-3
# logT=-5

# Function for perturbation
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
            Result[i]=(x[i]+y[i]+t[i])^2 + log(t[i]+1)
        }
    }
    
    return(Result)
}

xbound<-x[!TypePoint]
ybound<-y[!TypePoint]

# Per il plot ad un punto fisso reale
xp<-3
yp<-(-0.5)
tseq<-seq(min(TimePoints),max(TimePoints),length.out = 100)
ftseq<-fun(rep(xp,100),rep(yp,100),tseq,xbound,ybound)

##### PERTURBAZIONE TEMPORALE DEI DATI #####

DataMatrix<-NULL

for(i in TimePoints)
{
    Time<-rep(i,length(x[TypePoint]))
    if(noise)
    {
        DataMatrix<-cbind(DataMatrix,(fun(x[TypePoint],y[TypePoint],Time,xbound,ybound)+rnorm(length(Data),0,0.005)))
    }
    else
    {
        DataMatrix<-cbind(DataMatrix,fun(x[TypePoint],y[TypePoint],Time,xbound,ybound))
    }
}


# CERCO I VALORI DELLA FUNZIONE REALE NEI TEMPI SCELTI

nx<-length(xvec)
ny<-length(yvec)
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))

ResultMatrix<-NULL
for(j in TimePoints)
{
    Time<-rep(j,length(xx))
    Result<-fun(xx,yy,Time,xbound,ybound)
    ResultMatrix<-cbind(ResultMatrix,Result)
}

# Aggiorno il massimo per il grafico
MaxD=max(ResultMatrix,na.rm=TRUE)
MinD=min(ResultMatrix,na.rm=TRUE)



##### RISOLVO IL SISTEMA #####

#Creo le basi in spazio e tempo
TimeBasisObj<-Create.Bspline.Time.Basis(TimePoints,4,PlotIt=F)
SpaceBasisObj<-Create.FEM.Space.Basis(cbind(x,y),Triang,TypePoint,1)

for(logS in -6:+1)
{
    for(logT in -6:1)
    {
        LambdaS=10^logS
        LambdaT=10^logT
        #Introduco anche la misurazione del tempo per la funzione chiave
        SolutionObj<-smooth.ST.fd(DataMatrix,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT,2)
        
        # Ora salvo i risultati
        write.table(SolutionObj$C,file=paste(logS, logT,"MatriceC.txt",sep=" "),row.names=FALSE,col.names=FALSE)
                
        ##### GRAFICI #####
        
        PlotMatrix=matrix(ncol=length(yvec),nrow=length(xvec))
           
        # Voglio ricreare perfettamente le posizioni della image matrix tru

        # RICERCA DEL MASSIMO
        #PER AVERE UNA SCALA DI COLORI PER IL GRAFICO UNIFORME IN TUTTI GLI ISTANTI DI TEMPO
    
        ResultFitted<-NULL
        
        for(j in TimePoints)
        {
            Time<-rep(j,length(xValid))
            Result<-eval.ST.fd(xValid,yValid,Time,SolutionObj)
            
            ResultFitted<-cbind(ResultFitted,Result)
        }
        
        MaxF=max(ResultFitted,na.rm=TRUE)
        MinF=min(ResultFitted,na.rm=TRUE)
        
        zlim<-c(min(MinD,MinF),max(MaxD,MaxF))
        

        # GRAFICI

        # Ora salvo tutti i grafici
        # Voglio ricreare perfettamente le posizioni della image matrix tru
        
        # Prima ricavo il grafico della funzione reale
        for(j in 1:length(TimePoints))
        {
            # Matrice con la reale
            Mat <- matrix(ResultMatrix[,j],nrow=nx,ncol=ny,byrow=F)
            # Plot
            png(filename=paste(logS, logT,TimePoints[j],"reale.png",sep=" "))
            image(xvec,yvec,Mat,zlim=zlim,main=paste("Funzione reale tempo ",TimePoints[j],", logS=",logS,", logT=",logT,sep=""))
            lines(x[TypePoint==FALSE],y[TypePoint==FALSE],lwd=3)
            contour(xvec,yvec,Mat,nlevels=10,add=TRUE)
            dev.off()
            
            
            # Matrice con la stimata
            for(i in 1:(dim(PosMatrix)[1]))
            {
                PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=ResultFitted[i,j]
            }
            # Plot
            png(filename=paste(logS, logT,TimePoints[j],"stimata.png",sep=" "))
            image(xvec,yvec,PlotMatrix,zlim=zlim,main=paste("Funzione stimata tempo ",TimePoints[j],", logS=",logS,", logT=",logT,sep=""))
            lines(x[TypePoint==FALSE],y[TypePoint==FALSE],lwd=3)
            contour(xvec,yvec,PlotMatrix,nlevels=10,add=TRUE)
            dev.off()
        }
        
        # Plot in un punto fisso
        png(filename=paste(logS, logT,"Plot per un punto fissato.png",sep=" "))
        FixedPointPlot(xp,yp,SolutionObj)
        points(tseq,ftseq,type='l',col="red")
        legend("bottomright", col = c("red","black"),legend = c("reale","stimata"),lty=1)
        dev.off()
    }
}

