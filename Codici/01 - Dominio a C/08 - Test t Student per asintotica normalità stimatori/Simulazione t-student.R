library(rgl)
library(mgcv)
library(SDMTools)
source("SpazioTempo.R")
load("FrontieraC.RData")
load("Validazione.RData")

set.seed(10)
# Devo aggiungere del rumore al dato?
noise<-TRUE

DOF<-5

NITER<-100

# Per quanti istanti di tempo eseguire la stima?
TimePoints<-seq(0,2*pi,length.out=9)
# Per la covariata
Mu<-0
Sigma<-1
Beta<-1

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




##### OGGETTI UTILI #####

TimeBasisObj<-Create.Bspline.Time.Basis(TimePoints,TimeOrder=4,DerivativeOrder=2,PlotIt=F)
SpaceBasisObj<-Create.FEM.Space.Basis(cbind(x,y),Triang,TypePoint,1)




##### GCV #####

# DataMatrix<-NULL
# for(t in TimePoints)
# {
#     if(noise)
#     {
#         DataMatrix<-cbind(DataMatrix,fun(xknot,yknot,rep(t,length(xknot)))+rt(length(xknot),DOF))
#     } else
#     {
#         DataMatrix<-cbind(DataMatrix,fun(xknot,yknot,rep(t,length(xknot))))
#     }
# }
# 
# Data<-NULL
# for(i in 1:(dim(DataMatrix)[1]))
# {
#     Data<-c(Data,DataMatrix[i,])
# }
# 
# LogS = seq(-2,1,by=0.5)
# LogT = seq(-5,-2,by=0.5)
# GCVResult<-ST.GCV(Data,SpaceBasisObj,TimeBasisObj,LogS,LogT)
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




##### SOLUZIONE #####

##### RISOLUZIONE DEL SISTEMA #####

ResultT<-NULL
xx<-rep(xknot,length(TimePoints))
yy<-rep(yknot,length(TimePoints))
tt<-rep(TimePoints,each=length(xknot))
ResultBeta<-NULL

for(iter in 1:NITER)
{
    print(paste("Iterazione ",iter,"/",NITER,sep=""))
    DataMatrix<-NULL
    for(t in TimePoints)
    {
        if(noise)
        {
            DataMatrix<-cbind(DataMatrix,fun(xknot,yknot,rep(t,length(xknot)))+rt(length(xknot),DOF))
        } else
        {
            DataMatrix<-cbind(DataMatrix,fun(xknot,yknot,rep(t,length(xknot))))
        }
    }
    Data<-NULL
    for(i in 1:(dim(DataMatrix)[1]))
    {
        Data<-c(Data,DataMatrix[i,])
    }
    DesMat=rnorm(length(Data),Mu,Sigma)
    Data<-Data+Beta*DesMat
    
    SolutionObj<-ST.Smooth.Covar(Data,DesMat,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)
    tmp<-ST.Eval(xx,yy,tt,SolutionObj)
    ResultT<-cbind(ResultT,tmp)
    ResultBeta<-c(ResultBeta,SolutionObj$BetaHat)
}


p<-NULL
for(i in 1:(dim(ResultT)[1]))
{
    p<-c(p,shapiro.test(ResultT[i,])$p.value)
}
p<-c(p,shapiro.test(ResultBeta)$p.value)

padj<-p.adjust(p,method = "fdr")
pf<-padj[1:(length(padj)-1)]
png("p-values.png")
plot(pf,xlab="Punti",ylab="p-values",main="p-values per le stime puntuali di f")
dev.off()
length(pf[pf>0.05])/length(pf))
padj[length(padj)]
save(file="Simulazione.RData",ResultT,ResultBeta)
