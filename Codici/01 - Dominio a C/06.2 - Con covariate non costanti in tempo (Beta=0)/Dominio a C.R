library(rgl)
library(mgcv)
library(SDMTools)
source("SpazioTempo.R")
load("FrontieraC.RData")
load("Validazione.RData")

# Devo aggiungere del rumore al dato?
noise<-TRUE

# Per quanti istanti di tempo eseguire la stima?
TimePoints<-0:5

# Per la covariata
Mu<-5
Sigma<-1

# Plot Per un punto fissato
xp<-3
yp<--0.5
# Function for perturbation
fun=function(t)
{
    y=cos(t)
    return(t)
}



##### PERTURBAZIONE TEMPORALE DEI DATI #####

DataMatrix<-NULL

for(i in TimePoints)
{
    Time<-rep(i,length(x[TypePoint]))
    if(noise)
    {
        DataMatrix<-cbind(DataMatrix,Data*fun(i)+rnorm(length(Data),0,0.05))
    } else
    {
        DataMatrix<-cbind(DataMatrix,Data*fun(i))
    }
}

DataVec=NULL
for(i in 1:(dim(DataMatrix)[1]))
{
    DataVec<-c(DataVec,DataMatrix[i,])
}



##### MATRICE DISEGNO #####

# Voglio la covariata non fissa nel tempo
# Genero la covariata da una normale
DesMat<-rnorm(length(DataVec),Mu,Sigma)





##### GCV #####

#Creo le basi in spazio e tempo
TimeBasisObj<-Create.Bspline.Time.Basis(TimePoints,TimeOrder=4,DerivativeOrder=2,PlotIt=F)
SpaceBasisObj<-Create.FEM.Space.Basis(cbind(x,y),Triang,TypePoint,1)

# LogS<--8:+1
# LogT<--8:+1
# GCVResult<-ST.GCV.Covar(DataVec,DesMat,SpaceBasisObj,TimeBasisObj,LogS,LogT)
# 
# png(filename="GCV Matrix.png")
# image(LogS,LogT,GCVResult$GCVMatrix,col=heat.colors(100),main="GCV Matrix",xlab="logLambdaS",ylab="logLambdaT")
# dev.off()
# 
# LambdaS=10^GCVResult$Best[1]
# LambdaT=10^GCVResult$Best[2]
# 
# save(file="GCVResult.RData",GCVResult,LogS,LogT)

# Con Normale(0,5)
LambdaS=10^-3
LambdaT=10^1
# Con Normale(5,1)
LambdaS=10^-4
LambdaT=10^1

##### RISOLUZIONE DEL SISTEMA #####

SolutionObj<-ST.Smooth(DataVec,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)

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
        
MaxF=max(ResultFitted,na.rm=TRUE)
MinF=min(ResultFitted,na.rm=TRUE)
        
zlim<-c(MinF,MaxF)
        

# GRAFICI

# Ora salvo tutti i grafici
# Voglio ricreare perfettamente le posizioni della image matrix tru

nx<-length(xvec)
ny<-length(yvec)
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))

for(j in 1:length(TimePoints))
{
    # Matrice con la stimata
    for(i in 1:(dim(PosMatrix)[1]))
    {
        PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=ResultFitted[i,j]
    }
    # Plot
    png(filename=paste("Tempo",TimePoints[j],"stimata.png",sep=" "))
    image(xvec,yvec,PlotMatrix,zlim=zlim,main=paste("Funzione stimata tempo ",TimePoints[j],sep=""))
    lines(x[TypePoint==FALSE],y[TypePoint==FALSE],lwd=3)
    contour(xvec,yvec,PlotMatrix,nlevels=10,add=TRUE)
    dev.off()
}
        
# Plot in un punto fisso
png(filename=paste("Plot per un punto fissato.png",sep=" "))
FixedPointPlot(xp,yp,SolutionObj)
dev.off()



##### INTERVALLO DI CONFIDENZA #####
ICResult<-ST.IC(DataVec,DesMat,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)
save(file="ICResult.RData",ICResult)