library(rgl)
library(mgcv)
library(SDMTools)
library(ggplot2)
source("ST-PDE.R")

# Read objects
# Boundary points
Bound=read.table(file="DomC_Boundary.txt",header=T)
# Space points (internal points of the domain, where I will have data)
SpacePoints=read.table(file="DomC_InternalPoints.txt",header=T)
# Triangulation
Triang=read.table(file="DomC_Triangulation.txt",header=F)
TriangPoints=read.table(file="DomC_TriangulationPoints.txt",header=T)

# Time points
TimePoints<-seq(0,2*pi,length.out=9)

# Noise sd
set.seed(1)
SigmaNoise=0.5
# Gaussian Covariate Parameters
Mu<-0
Sigma<-1
Beta<-1

# Real Function 
fun=function(x,y,t)
{
    y=fs.test(x,y)*cos(t)
    return(y)
}



##### BASIS CREATION #####
# Time basis 
TimeOrder=4                 # Cubic B-Splines (use order +1)
NumBasis=length(TimePoints)  # How many time basis? I use the same number as TimePoints
TimeBasisObj = create.bspline.basis(c(min(TimePoints),max(TimePoints)), NumBasis, TimeOrder)
# Use fda function!
plot(TimeBasisObj,main=paste(NumBasis,"Time B-Spline Basis"),xlab="Time",ylab="Basis")

# Space basis
FEOrder=1   #Linear Finite Elements
SpaceBasisObj<-create.FEM.basis(TriangPoints,e=NULL,Triang,FEOrder)


##### DATA AND DESIGN MATRIX CREATION ####
# Data -> Real function + gaussian Noise
# Data Vector
# If i is Space index (in SpacePoints), and j Time Index (in TimePoints)
# DATA MUST BE ORDERED FIRST VARYING TIME INDEX, THAN SPACE INDEX
# z11, z12, z13, ..., z1m, z21, z22, ..., z2m, ..., znm
DataMatrix<-NULL
for(t in TimePoints)
{
    DataMatrix<-cbind(DataMatrix,fun(SpacePoints$xInt,SpacePoints$yInt,rep(t,length(SpacePoints$xInt)))+rnorm(length(SpacePoints$xInt),0,SigmaNoise))
}
Data<-NULL
for(i in 1:(dim(DataMatrix)[1]))
{
    Data<-c(Data,DataMatrix[i,])
}


# Design Matrix
# Use gaussian distribution
DesMat<-rnorm(length(Data),Mu,Sigma)
# In covariate case, I will add it to Data using Beta



##### MODEL WITHOUT COVARIATES #####

LogLambdaS=-0.5
LogLambdaT=-3.25
SolutionObj<-ST.Smooth(Data=Data,
                       #DESMAT MUST BE A NA
                       DesMat=NA,
                       SpacePoints=SpacePoints,
                       SpaceBasisObj=SpaceBasisObj,
                       TimePoints=TimePoints,
                       TimeBasisObj=TimeBasisObj,
                       LogLambdaS=LogLambdaS,
                       LogLambdaT=LogLambdaT)

# Save results
write.table(SolutionObj$cHat,file="cHat (without covariate).txt",row.names=FALSE,col.names=FALSE)

## FUNCTION PLOT AT GIVEN TIMES

nx<-200
ny<-200
xvec<-seq(min(Bound$xBound),max(Bound$xBound),length.out=nx)
yvec<-seq(min(Bound$yBound),max(Bound$yBound),length.out=ny)
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))

ResultFitted<-NULL

ntot=length(TimePoints)+1
k=1
pb <- txtProgressBar(min = 1, max = ntot, style = 3)

for(j in TimePoints)
{
    if(k==1)
    {
        setTxtProgressBar(pb, k)
    }
    Time<-rep(j,length(xx))
    Result<-ST.Eval(xx,yy,Time,SpaceBasisObj,TimeBasisObj,SolutionObj)
    
    ResultFitted<-cbind(ResultFitted,Result)
    k=k+1
    setTxtProgressBar(pb, k)
}

save(file="ResultFitted.RData",ResultFitted)
load(file="ResultFitted.RData")

zlim=c(min(ResultFitted,na.rm=T),max(ResultFitted,na.rm=T))

for(j in 1:length(TimePoints))
{
    # Matrix
    Mat <- matrix(ResultFitted[,j],nrow=nx,ncol=ny,byrow=F)
    Mat<-cbind(Mat,numeric(dim(Mat)[1]))
    Mat[1,ny+1]=zlim[1]
    Mat[2,ny+1]=zlim[2]
    
    # Plot
    png(filename=paste("Time",round(TimePoints[j],2),".png",sep=""))
    image(xvec,c(yvec,50),Mat,col=heat.colors(100),main=paste("Function in time ",round(TimePoints[j],2),sep=""),ylim=c(yvec[1],yvec[ny]),xlab="",ylab="")
    lines(Bound$xBound,Bound$yBound,lwd=1)
    lines(c(Bound$xBound[1],Bound$xBound[length(Bound$xBound)]),c(Bound$yBound[1],Bound$yBound[length(Bound$yBound)]),lwd=1)
    contour(xvec,c(yvec,50),Mat,nlevels=10,add=TRUE)
    dev.off()
}

# Try also fixed time plot
FixedTimePlot(pi,SpaceBasisObj,TimeBasisObj,SolutionObj)

## FUNCTION IN GIVEN POINTS
Selected<-c(10,28,67)
for(i in Selected)
{
    xP<-SpacePoints[i,1]
    yP<-SpacePoints[i,2]
    png(filename=paste("Plot in fixed point ",i,".png",sep=" "))
    FixedPointPlot(xP,yP,SpaceBasisObj,TimeBasisObj,SolutionObj,lwd=2,NameLocation = paste("(",round(xP,2),",",round(yP,2),")",sep=""),ylim=c(-4.4,4.4))
    points(TimePoints,DataMatrix[i,],pch=16,col="red")
    points(seq(min(TimePoints),max(TimePoints),length.out=100),fs.test(xP,yP)*cos(seq(min(TimePoints),max(TimePoints),length.out=100)),type='l',col="blue",lwd=2)
    legend("bottomleft",c("real", "fitted"), lty = c(1,1),col=c("blue","black"),lwd=2)
    dev.off()
}







##### MODEL WITH COVARIATES #####

# Add DesMat contribute
Data=Data+Beta*DesMat
LogLambdaS=-0.5
LogLambdaT=-3.25
SolutionObjCovar<-ST.Smooth(Data=Data,
                       DesMat=DesMat,
                       SpacePoints=SpacePoints,
                       SpaceBasisObj=SpaceBasisObj,
                       TimePoints=TimePoints,
                       TimeBasisObj=TimeBasisObj,
                       LogLambdaS=LogLambdaS,
                       LogLambdaT=LogLambdaT)

# Beta?
SolutionObjCovar$BetaHat

# Save results
write.table(SolutionObjCovar$cHat,file="cHat (with covariate).txt",row.names=FALSE,col.names=FALSE)
write.table(SolutionObjCovar$BetaHat,file="BetaHat (with covariate).txt",row.names=FALSE,col.names=FALSE)


## Confidence interval
IC=ST.IC(Data,DesMat,SpaceBasisObj,TimeBasisObj,LogLambdaS,LogLambdaT,Alpha=0.05)


## FUNCTION PLOT AT GIVEN TIMES

nx<-200
ny<-200
xvec<-seq(min(Bound$xBound),max(Bound$xBound),length.out=nx)
yvec<-seq(min(Bound$yBound),max(Bound$yBound),length.out=ny)
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))

ResultFitted<-NULL

ntot=length(TimePoints)+1
k=1
pb <- txtProgressBar(min = 1, max = ntot, style = 3)

for(j in TimePoints)
{
    if(k==1)
    {
        setTxtProgressBar(pb, k)
    }
    Time<-rep(j,length(xx))
    Result<-ST.Eval(xx,yy,Time,SpaceBasisObj,TimeBasisObj,SolutionObjCovar)
    
    ResultFitted<-cbind(ResultFitted,Result)
    k=k+1
    setTxtProgressBar(pb, k)
}

save(file="ResultFitted.RData",ResultFitted)
load(file="ResultFitted.RData")

zlim=c(min(ResultFitted,na.rm=T),max(ResultFitted,na.rm=T))

for(j in 1:length(TimePoints))
{
    # Matrix
    Mat <- matrix(ResultFitted[,j],nrow=nx,ncol=ny,byrow=F)
    Mat<-cbind(Mat,numeric(dim(Mat)[1]))
    Mat[1,ny+1]=zlim[1]
    Mat[2,ny+1]=zlim[2]
    
    # Plot
    png(filename=paste("Time",round(TimePoints[j],2),".png",sep=""))
    image(xvec,c(yvec,50),Mat,col=heat.colors(100),main=paste("Function in time ",round(TimePoints[j],2),sep=""),ylim=c(yvec[1],yvec[ny]),xlab="",ylab="")
    lines(Bound$xBound,Bound$yBound,lwd=1)
    lines(c(Bound$xBound[1],Bound$xBound[length(Bound$xBound)]),c(Bound$yBound[1],Bound$yBound[length(Bound$yBound)]),lwd=1)
    contour(xvec,c(yvec,50),Mat,nlevels=10,add=TRUE)
    dev.off()
}

# Try also fixed time plot
FixedTimePlot(pi,SpaceBasisObj,TimeBasisObj,SolutionObjCovar)

## FUNCTION IN GIVEN POINTS
Selected<-c(10,28,67)
for(i in Selected)
{
    xP<-SpacePoints[i,1]
    yP<-SpacePoints[i,2]
    png(filename=paste("Plot in fixed point ",i,".png",sep=" "))
    FixedPointPlot(xP,yP,SpaceBasisObj,TimeBasisObj,SolutionObjCovar,lwd=2,NameLocation = paste("(",round(xP,2),",",round(yP,2),")",sep=""),ylim=c(-4.4,4.4))
    points(TimePoints,DataMatrix[i,],pch=16,col="red")
    points(seq(min(TimePoints),max(TimePoints),length.out=100),fs.test(xP,yP)*cos(seq(min(TimePoints),max(TimePoints),length.out=100)),type='l',col="blue",lwd=2)
    legend("bottomleft",c("real", "fitted"), lty = c(1,1),col=c("blue","black"),lwd=2)
    dev.off()
}
