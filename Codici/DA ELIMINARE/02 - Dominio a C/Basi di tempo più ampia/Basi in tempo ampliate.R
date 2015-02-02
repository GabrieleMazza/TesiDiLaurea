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
LambdaT=10^-3



##### PERTURBAZIONE TEMPORALE DEI DATI #####

load("FrontieraC.RData")
DataMatrix<-NULL
# Perturbo secondo la funzione cos(x)
for(i in TimePoints)
{
    if(noise)
    {
        DataMatrix<-cbind(DataMatrix,(Data*cos(i)+rnorm(length(Data),0,0.5)))
    }
    else
    {
        DataMatrix<-cbind(DataMatrix,Data*cos(i))
    }
}


##### RISOLVO IL SISTEMA #####

# Creo le basi in spazio e tempo
TimeBasisObj<-Create.Bspline.Time.Basis(TimePoints,4,10,T)
SpaceBasisObj<-Create.FEM.Space.Basis(cbind(x,y),Triang,TypePoint,1)

# Introduco anche la misurazione del tempo per la funzione chiave
C<-smooth.ST.fd(DataMatrix,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)

# Non faccio nient'altro, sono interessato alle misure di tempo