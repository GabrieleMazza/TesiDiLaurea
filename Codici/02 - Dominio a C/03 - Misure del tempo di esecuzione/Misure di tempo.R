# MISURAZIONI DEL TEMPO DI ESECUZIONE DEL CODICE

library(rgl)
library(mgcv)
library(fda)
library(Matrix)
source("SpazioTempo.R")

# Devo aggiungere del rumore al dato?
noise<-TRUE

title="sparseMatrix"

# Per quanti istanti di tempo eseguire la stima?
TimePoints<-0:5

# Lambda
LambdaS=10^-2
LambdaT=10^-3

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


##### RISOLVO IL SISTEMA #####

#Creo le basi in spazio e tempo
TimeBasisObj<-Create.Bspline.Time.Basis(TimePoints,4,F)
SpaceBasisObj<-Create.FEM.Space.Basis(cbind(x,y),Triang,TypePoint,1)

#Introduco anche la misurazione del tempo per la funzione chiave
SolutionObj<-smooth.ST.fd(DataMatrix,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT,2)

t1<-proc.time()
png(filename=paste("Tempo Fissato (",title,").png",sep=""))
FixedTimePlot(1,SolutionObj)
dev.off()
t2<-proc.time()
DeltaT<-t2-t1
print("Fixed Time Plot: ")
print(DeltaT)

t1<-proc.time()
png(filename=paste("Punto Fissato (",title,").png",sep=""))
FixedPointPlot(3,-0.5,SolutionObj)
dev.off()
t2<-proc.time()
DeltaT<-t2-t1
print("Fixed Point Plot: ")
print(DeltaT)


print(paste("Number of internal points: ",sum(TypePoint),"/",length(TypePoint),sep=""))
print(paste("Times: ",length(TimePoints),sep=""))
print(paste("Dimension of the matrix: ",length(TimePoints)*length(TypePoint),sep=""))
