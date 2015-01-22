# Simulate spatial data on a C-shaped domain
# Last modified 24 July 2010 by laura


library(rgl)
library(mgcv)


nrep=50

######################################################################
######################################################################
######################   DATA GENERATION   ###########################
######################################################################
######################################################################




# Plot the function, and its boundary
fsb <- list(fs.boundary())
nx<-250
ny<-100 
xvec <- seq(-1,4,length=nx)
yvec<-seq(-1,1,length=ny)
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))
iY=1:ny
tru <- matrix(fs.test(xx,yy),nx,ny) ## truth
image(xvec,yvec,tru,col=heat.colors(100),xlab="x",ylab="y",asp=1)
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=3)
contour(xvec,yvec,tru,levels=seq(-5,5,by=.25),add=TRUE)

persp3d(xvec,yvec,tru, col='red', alpha=1)


 
 
 
for (rep in 1:nrep) 
{ 
## Simulate some fitting data, inside boundary...
np = 200
n <-500
x <- runif(n)*5-1
y <-runif(n)*2-1
muvec <- fs.test(x,y,b=1)
ind <- inSide(fsb,x=x,y=y) ## remove outsiders
muvec <- muvec[ind]
x <- x[ind]
y <- y[ind] 

muvec = muvec[1:np]
x = x[1:np]
y = y[1:np] 

sigma=0.5
W <- rnorm(np,muvec,sigma) ## add noise


##
## Generate design matrix
##


desmat=matrix(0,nrow=np,ncol=2)

desmat[,1]=rnorm(np,3,1.5)
desmat[,2]=rnorm(np,7,5)

beta1=-1/2
beta2=1/5

par(mfrow=c(2,2))
hist(W)
hist(beta1*desmat[,1])
hist(beta2*desmat[,2])
hist(W+beta1*desmat[,1]+beta2*desmat[,2])



# the data are obtained as:  W + beta1 * desmat[,1] + beta2 * desmat[,2]



data=data.frame(x=x,y=y,dataNOCovar=W,dataCovar=(W+beta1*desmat[,1]+beta2*desmat[,2]),Covariate1=desmat[,1],Covariate2=desmat[,2])

###########write.table(data,paste("rep",rep,"_data_TimWood.txt",sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE)


}



# Plot sampled values, mean and residuals  
# high values = white
# low  values = black
#
#
#
#
#windows()
#par(mfrow=c(2,2))
#
#plot(fsb[[1]]$x,fsb[[1]]$y,type="l")
#title("Sampled values of surface")
#for (k in 1:np)
#{points(s[k,1],s[k,2],pch=19,col=gray((dataNOCovar[k]-min(dataNOCovar))/(max(dataNOCovar)-min(dataNOCovar))))  
#}
#
#
#plot(fsb[[1]]$x,fsb[[1]]$y,type="l")
#title("Sampled data (including covariates)")
#for (k in 1:np)
#{points(s[k,1],s[k,2],pch=19,col=gray((dataCovar[k,2]-min(dataCovar[,2]))/(max(dataCovar[,2])-min(dataCovar[,2]))))  
#}
#
#
#plot(fsb[[1]]$x,fsb[[1]]$y,type="l")
#title("true value surface")
#for (k in 1:np)
#{points(s[k,1],s[k,2],pch=19,col=gray((muvec[k]-min(muvec))/(max(muvec)-min(muvec))))
#}
#
#
#residuals=(dataNOCovar-muvec)
#
#plot(fsb[[1]]$x,fsb[[1]]$y,type="l")
#title("Sampled values of surface - true values of surface")
#for (k in 1:np)
#{points(s[k,1],s[k,2],pch=19,col=gray((residuals[k]-min(residuals))/(max(residuals)-min(residuals))))
#}
#





######################################################################
######################################################################
##################   END DATA GENERATION   ###########################
######################################################################
######################################################################

setwd("E:/Dati/canada/fdaR")

library(soap)
library(fields)
library(soap)

nrep=50


# Plot the function, and its boundary
fsb <- fs.boundary()
nx<-250
ny<-100 
xvec <- seq(-1,4,length=nx)
yvec<-seq(-1,1,length=ny)
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))
iY=1:ny
tru <- matrix(fs.test(xx,yy),nx,ny) ## truth
image(xvec,yvec,tru,col=heat.colors(100),xlab="x",ylab="y",asp=1)
lines(fsb$x,fsb$y,lwd=3)
contour(xvec,yvec,tru,levels=seq(-5,5,by=.25),add=TRUE)

persp3d(xvec,yvec,tru, col='red', alpha=1)



## create a boundary...
fsb <- list(fs.boundary())
 
 
# create some internal knots for soap film
knots <- data.frame(x=rep(seq(-.5,3,by=.5),4),
                    y=rep(c(-.6,-.3,.3,.6),rep(8,4)))
 
 
 



betahat.Tps   = matrix(0,nrow=nrep,ncol=2)
sigmahat.Tps  = numeric(nrep)
MSEf.Tps      = numeric(nrep)



betahat.Krig   = matrix(0,nrow=nrep,ncol=2)
sigmahat.Krig  = numeric(nrep)
MSEf.Krig      = numeric(nrep)



betahat.soap   = matrix(0,nrow=nrep,ncol=2)
sigmahat.soap  = numeric(nrep)
MSEf.soap      = numeric(nrep)


 
for (rep in 1:10) 
{ 
print(rep)

data = as.matrix(read.table(paste("rep",rep,"_data_TimWood.txt",sep=""), header=T))


np           = length(data[,1])
p            = data[,1:2]
dataNOCovar  =  data[,3]
dataCovar    =  cbind(1:np,data[,4])
desmat       =  data[,5:6]




### Tps

fit.Tps.Covar=Tps(x=p,Y=dataCovar[,2],Z=desmat)

betahat.Tps[rep,] = fit.Tps.Covar$d[4:5]
sigmahat.Tps[rep] = fit.Tps.Covar$shat.GCV


XYmat.Tps.Covar.error = matrix(NA,nrow=length(xvec),ncol=length(yvec))
for (indexX in 1:length(xvec))
{
for (indexY in iY[inSide(fsb,x=rep(xvec[indexX],ny),y=yvec)])
{
XYmat.Tps.Covar.error[indexX,indexY]          = predict( fit.Tps.Covar, t(c(xvec[indexX],yvec[indexY])),drop.Z=T) - tru[indexX,indexY]
}
}
MSEf.Tps[rep] = mean(na.omit(as.vector(XYmat.Tps.Covar.error^2)))




### Krig


fit.Krig.Covar=Krig(x=p,Y=dataCovar[,2],Z=desmat)

betahat.Krig[rep,] = fit.Krig.Covar$d[4:5]
sigmahat.Krig[rep] = fit.Krig.Covar$shat.GCV

XYmat.Krig.Covar.error = matrix(NA,nrow=length(xvec),ncol=length(yvec))
for (indexX in 1:length(xvec))
{
for (indexY in iY[inSide(fsb,x=rep(xvec[indexX],ny),y=yvec)])
{
XYmat.Krig.Covar.error[indexX,indexY] = predict( fit.Krig.Covar, t(c(xvec[indexX],yvec[indexY])),drop.Z=T) - tru[indexX,indexY]
}
}
MSEf.Krig[rep] = mean(na.omit(as.vector(XYmat.Krig.Covar.error^2)))




### Soap film smoother


x=p[,1]
y=p[,2]


fit.soap.Covar <- gam(dataCovar[,2]~desmat+s(x,y,k=40,bs="so",xt=list(bnd=fsb)),knots=knots)

betahat.soap[rep,] = fit.soap.Covar$coefficients[c(2,3)]
sigmahat.soap[rep] = sqrt(fit.soap.Covar$sig2)


sob <- smooth.construct2(s(x,y,bs="so",k=40,xt=list(bnd=fsb)),data=data.frame(x=x,y=y),knots=knots)

matpred=Predict.matrix2(sob,data=list(x=xx,y=yy))

XYmat.soap.Covar=matrix((matpred%*%fit.soap.Covar$coefficients[-c(2,3)]),nx,ny)
XYmat.soap.Covar.error= XYmat.soap.Covar- tru

MSEf.soap[rep] = mean(na.omit(as.vector(XYmat.soap.Covar.error^2)))

}


results.Tps  = data.frame(betahat=betahat.Tps,sigmahat=sigmahat.Tps,MSEf=MSEf.Tps)
results.Krig = data.frame(betahat=betahat.Krig,sigmahat=sigmahat.Krig,MSEf=MSEf.Krig)
results.soap = data.frame(betahat=betahat.soap,sigmahat=sigmahat.soap,MSEf=MSEf.soap)



#######write.table(results.Tps,"results.Tps.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)
#######write.table(results.Krig,"results.Krig.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)
#######write.table(results.soap,"results.soap.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)





#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################


### Finite elements splines

source("AllFunctionsFEM.R")
library(fda)


betahat.FES   = matrix(0,nrow=nrep,ncol=2)
sigmahat.FES  = numeric(nrep)
MSEf.FES      = numeric(nrep)


Cconstredges=read.table("boundaryedges_TimWood.txt", header=F)
colnames(Cconstredges)=c("x","y")
 
 
 
#loglambdaseq = seq(-1,0.6,by=0.2)

#loglambdaseq = seq(-1.6,0.6,by=0.2) # for rep=c(13,30,35,37,40)

loglambdaseq = seq(-1,1.4,by=0.2)   # for rep =23
 
for (rep in 23) 
{ 
print(rep)

data = as.matrix(read.table(paste("rep",rep,"_data_TimWood.txt",sep=""), header=T))


np           = length(data[,1])
p            = data[,1:2]
dataNOCovar  =  data[,3]
dataCovar    =  cbind(1:np,data[,4])
desmat       =  data[,5:6]

p2col        = rbind(p,Cconstredges)

Tri = as.matrix(read.table(paste("rep",rep,"_t_TimWood.txt",sep=""), header=F))



#  The triangle index matrix


nt = dim(Tri)[[1]]

windows() 
plot(fsb[[1]]$x,fsb[[1]]$y,type="l",asp=1,main=rep)
title("Triangulation")
for (ne in 1:nt)
{polygon(c(p2col[Tri[ne,1],1],p2col[Tri[ne,2],1],p2col[Tri[ne,3],1]),c(p2col[Tri[ne,1],2],p2col[Tri[ne,2],2],p2col[Tri[ne,3],2]))
}

#  we don't need an edge matrix for these analyses

e = NULL

#  set up the FEM basis object and plot it

basisobj = create.FEM.basis(p2col, e, Tri)


#  set up the nodes, consisting of vertices plus edge midpoints
mknodes= makenodes(p2col,Tri)
nodes=mknodes$nodes
nodemesh=mknodes$nodemesh
nNodes = dim(nodes)[[1]]




#  set up a dummy FEM functional data object

simfd = fd(numeric(nNodes),basisobj)




#  smooth the data using covariates

GCVres = smooth.FEM.fd.Covar.GCV(dataCovar,desmat,simfd,loglambdaseq,0.95)
windows() 
plot(loglambdaseq, GCVres$GCVseq,main=rep)



loglambdaCovaroptim = loglambdaseq[which.min(GCVres$GCVseq)]

betahat.FES[rep,] = as.vector(GCVres$approxCIbeta[,"betahat"])
sigmahat.FES[rep] = GCVres$sigmahatseq[which.min(GCVres$GCVseq)]



lambdaCovar = 10^loglambdaCovaroptim
fit.FES.Covar = smooth.FEM.fd.Covar(dataCovar,desmat,simfd,lambdaCovar)


XYmat.FES.Covar       = matrix(NA,nrow=length(xvec),ncol=length(yvec))
XYmat.FES.Covar.error = XYmat.FES.Covar
for (indexX in 1:length(xvec))
{
for (indexY in iY[inSide(fsb,x=rep(xvec[indexX],ny),y=yvec)])
{
XYmat.FES.Covar[indexX,indexY]      = eval.FEM.fd(xvec[indexX], yvec[indexY], fit.FES.Covar$newfdobj) 
XYmat.FES.Covar.error[indexX,indexY] = XYmat.FES.Covar[indexX,indexY] - tru[indexX,indexY]
}
}

MSEf.FES[rep] = mean(na.omit(as.vector(XYmat.FES.Covar.error^2)))
#graphics.off()

}


results.FES = data.frame(betahat=betahat.FES,sigmahat=sigmahat.FES,MSEf=MSEf.FES)
#######write.table(results.FES,"results.FES.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)



























#########################################################
#########################################################
#########################################################
#########################################################


results.Tps=read.table("results.Tps.txt", header=TRUE)
results.Krig=read.table("results.Krig.txt", header=TRUE)
results.soap=read.table("results.soap.txt", header=TRUE)
results.FES=read.table("results.FES.txt", header=TRUE)




boxplot(sqrt(results.Krig$MSEf),sqrt(results.Tps$MSEf),sqrt(results.soap$MSEf),sqrt(results.FES$MSEf),names=c("Krig", "Tps", "Soap", "FES"),main="Root Mean Square Error f")


#boxplot(log10(sqrt(results.Krig$MSEf)),log10(sqrt(results.Tps$MSEf)),log10(sqrt(results.soap$MSEf)),log10(sqrt(results.FES$MSEf)),names=c("Krig", "Tps", "Soap", "FES"),main="Root Mean Square Error f")


sqrt(mean((results.Krig$betahat.1+0.5)^2))
sqrt(mean((results.Tps$betahat.1+0.5)^2))
sqrt(mean((results.soap$betahat.1+0.5)^2))
sqrt(mean((results.FES$betahat.1+0.5)^2))


sqrt(mean((results.Krig$betahat.2-0.2)^2))
sqrt(mean((results.Tps$betahat.2-0.2)^2))
sqrt(mean((results.soap$betahat.2-0.2)^2))
sqrt(mean((results.FES$betahat.2-0.2)^2))


sqrt(mean((results.Krig$sigmahat-0.5)^2))
sqrt(mean((results.Tps$sigmahat-0.5)^2))
sqrt(mean((results.soap$sigmahat-0.5)^2))
sqrt(mean((results.FES$sigmahat-0.5)^2))






#########################################################
#########################################################
#### Figures draft


rep = 2

data = as.matrix(read.table(paste("rep",rep,"_data_TimWood.txt",sep=""), header=T))


np           = length(data[,1])
p            = data[,1:2]
dataNOCovar  =  data[,3]
dataCovar    =  cbind(1:np,data[,4])
desmat       =  data[,5:6]

p2col        = rbind(p,Cconstredges)

Tri = as.matrix(read.table(paste("rep",rep,"_t_TimWood.txt",sep=""), header=F))


plot(fsb[[1]]$x,fsb[[1]]$y,xlab="x",ylab="y",type="l",asp=1,lwd=3)
title("Sampled data")
for (k in 1:np)
{points(p[k,1],p[k,2],pch=19,cex=(dataCovar[k,2]-min(dataCovar[,2]))/(max(dataCovar[,2])-min(dataCovar[,2])))
}



windows()
par(mfrow=c(1,2))
plot(desmat[,1],dataCovar[,2],ylab="response variable", xlab="covariate 1")
plot(desmat[,2],dataCovar[,2],ylab="response variable", xlab="covariate 2")



##### Tps

library(fields)
fit.Tps.Covar=Tps(x=p,Y=dataCovar[,2],Z=desmat)


dvec.Tps.Covar = numeric(np*(np-1)/2)
rvec.Tps.Covar = dvec.Tps.Covar

mm = 0
for (i in 2:np)
 {   
   for (j in 1:(i-1))
     {  mm = mm + 1
        dvec.Tps.Covar[mm] = sum((p[i,]-p[j,])^2)
        rvec.Tps.Covar[mm] = fit.Tps.Covar$residuals[i] * fit.Tps.Covar$residuals[j]
     }
 }


windows()
plot(log10(dvec.Tps.Covar), rvec.Tps.Covar,xlab="log_{10} Squared distance")
abline(h=0,col=2,lwd=2,lty="dashed")

XYmat.Tps.Covar.trend.nonparam = matrix(NA,nrow=length(xvec),ncol=length(yvec))
XYmat.Tps.Covar.error          = XYmat.Tps.Covar.trend.nonparam 
for (indexX in 1:length(xvec))
{
for (indexY in iY[inSide(fsb,x=rep(xvec[indexX],ny),y=yvec)])
{
XYmat.Tps.Covar.trend.nonparam[indexX,indexY] = predict( fit.Tps.Covar, t(c(xvec[indexX],yvec[indexY])),drop.Z=T)                 # smoothing surface, trend + nonparametric
XYmat.Tps.Covar.error[indexX,indexY]          = XYmat.Tps.Covar.trend.nonparam[indexX,indexY] - tru[indexX,indexY]
}
}


##### Krig


fit.Krig.Covar=Krig(x=p,Y=dataCovar[,2],Z=desmat)

dvec.Krig.Covar = numeric(np*(np-1)/2)
rvec.Krig.Covar = dvec.Krig.Covar

mm = 0
for (i in 2:np)
 {   
   for (j in 1:(i-1))
     {  mm = mm + 1
        dvec.Krig.Covar[mm] = sum((p[i,]-p[j,])^2)
        rvec.Krig.Covar[mm] = fit.Krig.Covar$residuals[i] * fit.Krig.Covar$residuals[j]
     }
 }


windows()
plot(log10(dvec.Krig.Covar), rvec.Krig.Covar,xlab="log_{10} Squared distance")
abline(h=0,col=2,lwd=2,lty="dashed")


XYmat.Krig.Covar.trend.nonparam = matrix(NA,nrow=length(xvec),ncol=length(yvec))
XYmat.Krig.Covar.error          = XYmat.Krig.Covar.trend.nonparam 
for (indexX in 1:length(xvec))
{
for (indexY in iY[inSide(fsb,x=rep(xvec[indexX],ny),y=yvec)])
{
XYmat.Krig.Covar.trend.nonparam[indexX,indexY] = predict( fit.Krig.Covar, t(c(xvec[indexX],yvec[indexY])),drop.Z=T)                 # smoothing surface, trend + nonparametric
XYmat.Krig.Covar.error[indexX,indexY]          = XYmat.Krig.Covar.trend.nonparam[indexX,indexY] - tru[indexX,indexY]
}
}




### Soap film smoother


x=p[,1]
y=p[,2]

knots <- data.frame(x=rep(seq(-.5,3,by=.5),4),
                    y=rep(c(-.6,-.3,.3,.6),rep(8,4)))

fit.soap.Covar <- gam(dataCovar[,2]~desmat+s(x,y,k=40,bs="so",xt=list(bnd=fsb)),knots=knots)

sob <- smooth.construct2(s(x,y,bs="so",k=40,xt=list(bnd=fsb)),data=data.frame(x=x,y=y),knots=knots)

matpred=Predict.matrix2(sob,data=list(x=xx,y=yy))

XYmat.soap.Covar=matrix((matpred%*%fit.soap.Covar$coefficients[-c(2,3)]),nx,ny)
XYmat.soap.Covar.error= XYmat.soap.Covar- tru


dvec.soap.Covar = numeric(np*(np-1)/2)
rvec.soap.Covar = dvec.soap.Covar

mm = 0
for (i in 2:np)
 {   
   for (j in 1:(i-1))
     {  mm = mm + 1
        dvec.soap.Covar[mm] = sum((p[i,]-p[j,])^2)
        rvec.soap.Covar[mm] = fit.soap.Covar$residuals[i] * fit.soap.Covar$residuals[j]
     }
 }


windows()
plot(log10(dvec.soap.Covar), rvec.soap.Covar,xlab="log_{10} Squared distance")
abline(h=0,col=2,lwd=2,lty="dashed")





###### FES
library(fda)


Cconstredges=read.table("boundaryedges_TimWood.txt", header=F)
colnames(Cconstredges)=c("x","y")
 
 
Tri = as.matrix(read.table(paste("rep",rep,"_t_TimWood.txt",sep=""), header=F))

#  The triangle index matrix

nt = dim(Tri)[[1]]

windows() 
#plot(fsb[[1]]$x,fsb[[1]]$y,xlab="x",ylab="y",type="l",col="red",asp=1)
plot(fsb[[1]]$x,fsb[[1]]$y,xlab="x",ylab="y",type="n",col="red",asp=1)
title("Triangulation")
for (ne in 1:nt)
{polygon(c(p2col[Tri[ne,1],1],p2col[Tri[ne,2],1],p2col[Tri[ne,3],1]),c(p2col[Tri[ne,1],2],p2col[Tri[ne,2],2],p2col[Tri[ne,3],2]))
}

#  we don't need an edge matrix for these analyses

e = NULL

#  set up the FEM basis object and plot it

basisobj = create.FEM.basis(p2col, e, Tri)

#  set up the nodes, consisting of vertices plus edge midpoints
mknodes= makenodes(p2col,Tri)
nodes=mknodes$nodes
nodemesh=mknodes$nodemesh
nNodes = dim(nodes)[[1]]

#  set up a dummy FEM functional data object

simfd = fd(numeric(nNodes),basisobj)

loglambdaseq = seq(-1,0.6,by=0.2)

GCVres = smooth.FEM.fd.Covar.GCV(dataCovar,desmat,simfd,loglambdaseq,0.95)
windows() 
plot(loglambdaseq, GCVres$GCVseq,main=rep)



loglambdaCovaroptim = loglambdaseq[which.min(GCVres$GCVseq)]

shapiro.test(Res.fit.FES.Covar[,2])
GCVres$approxCIbeta

lambdaCovar = 10^loglambdaCovaroptim
fit.FES.Covar = smooth.FEM.fd.Covar(dataCovar,desmat,simfd,lambdaCovar)

fhat.Covar=fit.FES.Covar$newfdobj$coef
betahat.Covar= ( solve( t(desmat) %*% desmat ) ) %*% t(desmat) %*% (dataCovar[,2]-fhat.Covar[1:np,])

#  evaluate the solution on the points

eval.fit.FES.Covar = eval.FEM.fd(p[,1], p[,2], fit.FES.Covar$newfdobj) +  desmat %*% betahat.Covar

#  compute the residuals

Res.fit.FES.Covar = cbind(dataCovar[,1], dataCovar[,2] - eval.fit.FES.Covar)

#  compute the squared distances and the residual products

dvec.FES.Covar = numeric(np*(np-1)/2)
rvec.FES.Covar = dvec.FES.Covar

m = 0
for (i in 2:np)
 {   
   for (j in 1:(i-1))
     {  m = m + 1
        dvec.FES.Covar[m] = sum((p[i,]-p[j,])^2)
        rvec.FES.Covar[m] = Res.fit.FES.Covar[i,2] * Res.fit.FES.Covar[j,2]
     }
 }


windows()
plot(log10(dvec.FES.Covar), rvec.FES.Covar,xlab="log_{10} Squared distance")
abline(h=0,col=2,lwd=2,lty="dashed")



XYmat.FES.Covar       = matrix(NA,nrow=length(xvec),ncol=length(yvec))
XYmat.FES.Covar.error = XYmat.FES.Covar
for (indexX in 1:length(xvec))
{
for (indexY in iY[inSide(fsb,x=rep(xvec[indexX],ny),y=yvec)])
{
XYmat.FES.Covar[indexX,indexY]      = eval.FEM.fd(xvec[indexX], yvec[indexY], fit.FES.Covar$newfdobj) 
XYmat.FES.Covar.error[indexX,indexY] = XYmat.FES.Covar[indexX,indexY] - tru[indexX,indexY]
}
}






###############################################################
###############################################################


beta1=-1/2
beta2=1/5


 
windows() 
zlimits=c(min(c(tru[which.min(tru)],XYmat.soap.Covar[which.min(XYmat.soap.Covar)],XYmat.Krig.Covar.trend.nonparam[which.min(XYmat.Krig.Covar.trend.nonparam)],XYmat.Tps.Covar.trend.nonparam[which.min(XYmat.Tps.Covar.trend.nonparam)],XYmat.FES.Covar[which.min(XYmat.FES.Covar)])),max(c(tru[which.max(tru)],XYmat.soap.Covar[which.max(XYmat.soap.Covar)],XYmat.Krig.Covar.trend.nonparam[which.max(XYmat.Krig.Covar.trend.nonparam)],XYmat.Tps.Covar.trend.nonparam[which.max(XYmat.Tps.Covar.trend.nonparam)],XYmat.FES.Covar[which.max(XYmat.FES.Covar)])))

par(mfrow=c(2,2))

betahatKrig=round(1000*fit.Krig.Covar$d[(length(fit.Krig.Covar$d)-1):length(fit.Krig.Covar$d)])/1000
image(xvec,yvec,XYmat.Krig.Covar.trend.nonparam,main=paste(c("Filtered kriging, betahat = " , betahatKrig)),zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.Krig.Covar.trend.nonparam,zlim=zlimits,add=T)
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=3)

betahatTps=round(1000*fit.Tps.Covar$d[4:5])/1000
image(xvec,yvec,XYmat.Tps.Covar.trend.nonparam,main=paste(c("Thin plate splines, betahat = " , betahatTps)),zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.Tps.Covar.trend.nonparam,zlim=zlimits,add=T)
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=3)

betahatsoap=round(1000*as.vector(fit.soap.Covar$coefficients[c(2,3)]))/1000
image(xvec,yvec,XYmat.soap.Covar,main=paste(c("Soap film smoothing, betahat = " , betahatsoap )),zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.soap.Covar,zlim=zlimits,add=T)
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=3)

betahatFES=round(1000*betahat.Covar[,1])/1000
image(xvec,yvec,XYmat.FES.Covar,main=paste(c("Finite Elements splines, betahat = " , betahatFES)),zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.FES.Covar,zlim=zlimits,add=T)   
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=3)










XYmat.Tps.Covar.error.abs=abs(XYmat.Tps.Covar.error)
XYmat.Krig.Covar.error.abs=abs(XYmat.Krig.Covar.error)
XYmat.FES.Covar.error.abs=abs(XYmat.FES.Covar.error)
XYmat.soap.Covar.error.abs=abs(XYmat.soap.Covar.error)


windows()
par(mfrow=c(2,2))
zlimits=c(min(c(XYmat.Krig.Covar.error.abs[which.min(XYmat.Krig.Covar.error.abs)],XYmat.soap.Covar.error.abs[which.min(XYmat.soap.Covar.error.abs)],XYmat.Tps.Covar.error.abs[which.min(XYmat.Tps.Covar.error.abs)],XYmat.FES.Covar.error.abs[which.min(XYmat.FES.Covar.error.abs)])),max(c(XYmat.Krig.Covar.error.abs[which.max(XYmat.Krig.Covar.error.abs)],XYmat.soap.Covar.error.abs[which.max(XYmat.soap.Covar.error.abs)],XYmat.Tps.Covar.error.abs[which.max(XYmat.Tps.Covar.error.abs)],XYmat.FES.Covar.error.abs[which.max(XYmat.FES.Covar.error.abs)])))

image(xvec,yvec,XYmat.Tps.Covar.error.abs,main="Thin plate splines",zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.Tps.Covar.error.abs,zlim=zlimits,add=T)

image(xvec,yvec,XYmat.Krig.Covar.error.abs,main="Filtered kriging",zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.Krig.Covar.error.abs,zlim=zlimits,add=T)

image(xvec,yvec,XYmat.soap.Covar.error.abs,main="Soap film smoothing",zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.soap.Covar.error.abs,zlim=zlimits,add=T)   

image(xvec,yvec,XYmat.FES.Covar.error.abs,main="Finite Elements splines",zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.FES.Covar.error.abs,zlim=zlimits,add=T)   



windows()
par(mfrow=c(2,2))
ylimits=c(min(c(rvec.soap.Covar,rvec.Tps.Covar,rvec.FES.Covar,rvec.Krig.Covar)),max(c(rvec.soap.Covar,rvec.Tps.Covar,rvec.FES.Covar,rvec.Krig.Covar)))

plot(log10(dvec.Krig.Covar), rvec.Krig.Covar,xlab="log_{10} Squared distance",main="Filtered kriging",ylab="",ylim=ylimits)
abline(h=0,col=2,lwd=2,lty="dashed")

plot(log10(dvec.Tps.Covar), rvec.Tps.Covar,xlab="log_{10} Squared distance",main="Thin plate splines",ylab="",ylim=ylimits)
abline(h=0,col=2,lwd=2,lty="dashed")

plot(log10(dvec.soap.Covar), rvec.soap.Covar,xlab="log_{10} Squared distance",main="Soap film smoothing",ylab="",ylim=ylimits)
abline(h=0,col=2,lwd=2,lty="dashed")

plot(log10(dvec.FES.Covar), rvec.FES.Covar,xlab="log_{10} Squared distance",main="Finite Elements splines",ylab="",ylim=ylimits)
abline(h=0,col=2,lwd=2,lty="dashed")









citation(package = "fields")












############ 
#### Pictures for IWFOS short paper

graphics.off()

par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.7,pty="m", font.main=1.1,font.lab=1.1, font.axis=1.1,cex.lab=1.1,cex.axis=0.9)


zlimits2=c(min(c(tru[which.min(tru)],XYmat.Tps.Covar.trend.nonparam[which.min(XYmat.Tps.Covar.trend.nonparam)],XYmat.FES.Covar[which.min(XYmat.FES.Covar)])),max(c(tru[which.max(tru)],XYmat.Tps.Covar.trend.nonparam[which.max(XYmat.Tps.Covar.trend.nonparam)],XYmat.FES.Covar[which.max(XYmat.FES.Covar)])))

betahatTps=round(1000*fit.Tps.Covar$d[4:5])/1000
image(xvec,yvec,XYmat.Tps.Covar.trend.nonparam,main="Thin-plate spline, surface estimate",zlim=zlimits2,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.Tps.Covar.trend.nonparam,zlim=zlimits2,add=T)
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=2)

betahatFES=round(1000*betahat.Covar[,1])/1000
image(xvec,yvec,XYmat.FES.Covar,main="Finite element L-spline, surface estimate",zlim=zlimits2,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.FES.Covar,zlim=zlimits2,add=T)   
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=2)


zlimits2e=c(min(c(XYmat.Tps.Covar.error.abs[which.min(XYmat.Tps.Covar.error.abs)],XYmat.FES.Covar.error.abs[which.min(XYmat.FES.Covar.error.abs)])),max(c(XYmat.Tps.Covar.error.abs[which.max(XYmat.Tps.Covar.error.abs)],XYmat.FES.Covar.error.abs[which.max(XYmat.FES.Covar.error.abs)])))

image(xvec,yvec,XYmat.Tps.Covar.error.abs,main="Thin-plate spline, absolute residuals",zlim=zlimits2e,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.Tps.Covar.error.abs,zlim=zlimits2e,add=T)

image(xvec,yvec,XYmat.FES.Covar.error.abs,main="Finite element L-spline, absolute residuals",zlim=zlimits2e,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.FES.Covar.error.abs,zlim=zlimits2e,add=T)   



dev.copy(device=postscript,file="IWFOS_TimWood_results3.ps",horizontal=FALSE,width=7.1, height=4.3)
dev.off()

windows()
par (layout(rbind(c(1,2,2), c(3,4,5)), widths=c(2,1,1), heights=c(1,1)),mar=c(6,5,2,1),mex=0.9,pty="m", font.main=2,font.lab=1.3, font.axis=1.3,cex.lab=1.3,cex.axis=1.1)

image(xvec,yvec,tru,main="True surface",zlim=zlimits2,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,tru,zlim=zlimits2,add=T)   
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=2)



plot(fsb[[1]]$x,fsb[[1]]$y,xlab="",ylab="",type="l",asp=1,lwd=2,xlim=c(-0.85,3.85))
title("Sampled data")
for (k in 1:np)
{points(p[k,1],p[k,2],pch=19,cex=(dataCovar[k,2]-min(dataCovar[,2]))/(max(dataCovar[,2])-min(dataCovar[,2])))
}



#plot(fsb[[1]]$x,fsb[[1]]$y,xlab="x",ylab="y",type="l",col="red",asp=1,lwd=2)
plot(fsb[[1]]$x,fsb[[1]]$y,xlab="",ylab="",type="n",col="red",asp=1,xlim=c(-0.85,3.85))
title("Triangulation")
for (ne in 1:nt)
{polygon(c(p2col[Tri[ne,1],1],p2col[Tri[ne,2],1],p2col[Tri[ne,3],1]),c(p2col[Tri[ne,1],2],p2col[Tri[ne,2],2],p2col[Tri[ne,3],2]))
}


plot(desmat[,1],dataCovar[,2],ylab="response variable", xlab="covariate 1")
plot(desmat[,2],dataCovar[,2],ylab="response variable", xlab="covariate 2")

dev.copy(device=postscript,file="IWFOS_TimWood_test2.ps",horizontal=FALSE,width=7, height=4.5)
dev.off()





par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.9,pty="m", font.main=1.1,font.lab=1.5, font.axis=1.3,cex.lab=1.5,cex.axis=1)
boxplot(sqrt(results.Tps$MSEf),sqrt(results.FES$MSEf),names=c("Thin-plate spline", "FELspline"), main="Root Mean Square Error f")
