# Simulate spatial data on a C-shaped domain
# Last modified 11 May 2011 by laura


setwd("E:/Dati/canada/fdaR")
library(rgl)
library(soap)

### Finite elements splines

source("SFDA_AllFunctions.R")
library(fda)



nrep=50


fs.test=function (x, y, r0 = 0.1, r = 0.5, l = 3, b = 1, exclude = TRUE) 
{
    q <- pi * r/2
    a <- d <- x * 0
    ind <- x >= 0 & y > 0
    a[ind] <- q + x[ind]
    d[ind] <- y[ind] - r
    ind <- x >= 0 & y <= 0
    a[ind] <- -q - x[ind]
    d[ind] <- -r - y[ind]
    ind <- x < 0
    a[ind] <- -atan(y[ind]/x[ind]) * r
    d[ind] <- sqrt(x[ind]^2 + y[ind]^2) - r
    ind <- abs(d) > r - r0 | (x > l & (x - l)^2 + d^2 > (r - 
        r0)^2)
    f <- a * b + d^2
    if (exclude) 
        f[ind] <- NA
    attr(f, "exclude") <- ind
    f
}


fs.boundary=function (r0 = 0.1, r = 0.5, l = 3, n.theta = 20) 
{
    rr <- r + (r - r0)
    theta <- seq(pi, pi/2, length = n.theta)
    x <- rr * cos(theta)
    y <- rr * sin(theta)
    theta <- seq(pi/2, -pi/2, length = 2 * n.theta)
    x <- c(x, (r - r0) * cos(theta) + l)
    y <- c(y, (r - r0) * sin(theta) + r)
    theta <- seq(pi/2, pi, length = n.theta)
    x <- c(x, r0 * cos(theta))
    y <- c(y, r0 * sin(theta))
    n <- length(x)
    x <- c(x, x[n:1])
    y <- c(y, -y[n:1])
    return(list(x = x, y = y))
}

inSide=function (bnd, x, y) 
{
    if (is.null(bnd$x)) {
        bn <- list(x = bnd[[1]]$x, y = bnd[[1]]$y)
        if (length(bnd) > 1) 
            for (i in 2:length(bnd)) {
                bn$x <- c(bn$x, NA, bnd[[i]]$x)
                bn$y <- c(bn$y, NA, bnd[[i]]$y)
            }
        bnd <- bn
    }
    lowLim <- min(c(bnd$x, bnd$y), na.rm = TRUE) - 1
    ind <- is.na(bnd$x) | is.na(bnd$y)
    bnd$x[ind] <- bnd$y[ind] <- lowLim - 1
    n <- length(bnd$x)
    if (n != length(bnd$y)) 
        stop("x and y must be same length")
    um <- .C("in_out", bx = as.double(bnd$x), by = as.double(bnd$y), 
        break.code = as.double(lowLim), x = as.double(x), y = as.double(y), 
        inside = as.integer(y * 0), nb = as.integer(n), n = as.integer(length(x)), 
        PACKAGE = "soap")
    as.logical(um$inside)
}






# Plot the function, and its boundary
fsb <- list(fs.boundary())
nx<-500
ny<-200 
xvec <- seq(-1,4,length=nx)
yvec<-seq(-1,1,length=ny)
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))
iY=1:ny
tru <- matrix(fs.test(xx,yy),nx,ny) ## truth


image(xvec,yvec,tru,col=heat.colors(100),xlab="x",ylab="y",zlab="",asp=1)
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=3)
contour(xvec,yvec,tru,levels=seq(-5,5,by=.25),add=TRUE)

persp3d(xvec,yvec,tru, col='red', xlab="",ylab="",zlab="", xlim=c(-1.1,3.5),alpha=1)



Cconstredges=read.table("boundaryedges_TimWood.txt", header=F)
colnames(Cconstredges)=c("x","y")
 











##################################
##################################
########  NO COVARIATES  #########
##################################
##################################





order=2 

sigmahat.FES  = numeric(nrep)
MSEf.FES      = numeric(nrep)
GCV.FES          = numeric(nrep)
EDF.FES          = numeric(nrep)
loglambdaoptim.FES  = numeric(nrep)

 
loglambdaseq = seq(-1.8,0.8,by=0.2)

#loglambdaseq = seq(-1.6,0.6,by=0.2) # for rep=c(13,30,35,37,40)

#loglambdaseq = seq(-1,1.4,by=0.2)   # for rep =23
 
for (rep in 1:nrep) 
{ 
print(rep)

data = as.matrix(read.table(paste("rep",rep,"_data_TimWood.txt",sep=""), header=T))


np           = length(data[,1])
p            = data[,1:2]
data = cbind(1:np,data[,3])

p2col        = rbind(p,Cconstredges)

Tri = as.matrix(read.table(paste("rep",rep,"_t_TimWood.txt",sep=""), header=F))

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

basisobj = create.FEM.basis(p2col, e, Tri, order)


#  set up a dummy FEM functional data object

simfd = fd(numeric(basisobj$nbasis),basisobj)




#  smooth the data without using covariates

GCVres = smooth.FEM.fd.GCV(data,simfd,loglambdaseq)

loglambdaoptim.FES[rep] = loglambdaseq[which.min(GCVres$GCVseq)]
GCV.FES[rep]          = min(GCVres$GCVseq)
EDF.FES[rep]          = GCVres$EDFseq[which.min(GCVres$GCVseq)]


sigmahat.FES[rep] = GCVres$sigmahatseq[which.min(GCVres$GCVseq)]


windows() 
plot(loglambdaseq, GCVres$GCVseq,main=rep)
abline(v=loglambdaoptim.FES[rep])


lambdaCovar = 10^loglambdaoptim.FES[rep] 
fit.FES.Covar = smooth.FEM.fd(data,simfd,lambdaCovar)


coefres=data.frame(coeffelsplobj=fit.FES.Covar$felsplobj$coef,coeflaplacefd=fit.FES.Covar$laplacefd$coef)
colnames(coefres)=c("coef.felsplobj","coef.laplacefd")

write.table(coefres,paste("rep",rep,"_results.FES.NOcovar.order2.txt",sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE)



XYmat.FES.Covar       = matrix(NA,nrow=length(xvec),ncol=length(yvec))
XYmat.FES.Covar.error = XYmat.FES.Covar
for (indexX in 1:length(xvec))
{
for (indexY in iY[inSide(fsb,x=rep(xvec[indexX],ny),y=yvec)])
{
XYmat.FES.Covar[indexX,indexY]      = eval.FEM.fd(xvec[indexX], yvec[indexY], fit.FES.Covar$felsplobj) 
XYmat.FES.Covar.error[indexX,indexY] = XYmat.FES.Covar[indexX,indexY] - tru[indexX,indexY]
}
}

MSEf.FES[rep] = mean(na.omit(as.vector(XYmat.FES.Covar.error^2)))
graphics.off()

}


results.FES = data.frame(loglambdaoptim=loglambdaoptim.FES, GCV=GCV.FES, EDF=EDF.FES, sigmahat=sigmahat.FES, MSEf=MSEf.FES)
write.table(results.FES,"results.FES.NOcovar.order2.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)







####################################
####################################
####################################




order=1 # diff

sigmahat.FES  = numeric(nrep)
MSEf.FES      = numeric(nrep)
GCV.FES          = numeric(nrep)
EDF.FES          = numeric(nrep)
loglambdaoptim.FES  = numeric(nrep)

 
loglambdaseq = seq(-1.8,0.8,by=0.2)

#loglambdaseq = seq(-1.6,0.6,by=0.2) # for rep=c(13,30,35,37,40)

#loglambdaseq = seq(-1,1.4,by=0.2)   # for rep =23
 
for (rep in 1:nrep) 
{ 
print(rep)

data = as.matrix(read.table(paste("rep",rep,"_data_TimWood.txt",sep=""), header=T))


np           = length(data[,1])
p            = data[,1:2]
data = cbind(1:np,data[,3])

p2col        = rbind(p,Cconstredges)

Tri = as.matrix(read.table(paste("rep",rep,"_t_TimWood.txt",sep=""), header=F))

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

basisobj = create.FEM.basis(p2col, e, Tri, order)


#  set up a dummy FEM functional data object

simfd = fd(numeric(basisobj$nbasis),basisobj)




#  smooth the data without using covariates

GCVres = smooth.FEM.fd.GCV(data,simfd,loglambdaseq)

loglambdaoptim.FES[rep] = loglambdaseq[which.min(GCVres$GCVseq)]
GCV.FES[rep]          = min(GCVres$GCVseq)
EDF.FES[rep]          = GCVres$EDFseq[which.min(GCVres$GCVseq)]


sigmahat.FES[rep] = GCVres$sigmahatseq[which.min(GCVres$GCVseq)]


windows() 
plot(loglambdaseq, GCVres$GCVseq,main=rep)
abline(v=loglambdaoptim.FES[rep])


lambdaCovar = 10^loglambdaoptim.FES[rep] 
fit.FES.Covar = smooth.FEM.fd(data,simfd,lambdaCovar)


coefres=data.frame(coeffelsplobj=fit.FES.Covar$felsplobj$coef)
colnames(coefres)=c("coef.felsplobj")

write.table(coefres,paste("rep",rep,"_results.FES.NOcovar.order1.txt",sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE)



XYmat.FES.Covar       = matrix(NA,nrow=length(xvec),ncol=length(yvec))
XYmat.FES.Covar.error = XYmat.FES.Covar
for (indexX in 1:length(xvec))
{
for (indexY in iY[inSide(fsb,x=rep(xvec[indexX],ny),y=yvec)])
{
XYmat.FES.Covar[indexX,indexY]      = eval.FEM.fd(xvec[indexX], yvec[indexY], fit.FES.Covar$felsplobj) 
XYmat.FES.Covar.error[indexX,indexY] = XYmat.FES.Covar[indexX,indexY] - tru[indexX,indexY]
}
}

MSEf.FES[rep] = mean(na.omit(as.vector(XYmat.FES.Covar.error^2)))
graphics.off()

}


results.FES = data.frame(loglambdaoptim=loglambdaoptim.FES, GCV=GCV.FES, EDF=EDF.FES, sigmahat=sigmahat.FES, MSEf=MSEf.FES)
write.table(results.FES,"results.FES.NOcovar.order1.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)





####################################
####################################
####################################





order=1 # NOdiff!!!!!!!!

sigmahat.FES  = numeric(nrep)
MSEf.FES      = numeric(nrep)
GCV.FES          = numeric(nrep)
EDF.FES          = numeric(nrep)
loglambdaoptim.FES  = numeric(nrep)

 
loglambdaseq = seq(-1.8,0.8,by=0.2)

#loglambdaseq = seq(-1.6,0.6,by=0.2) # for rep=c(13,30,35,37,40)

#loglambdaseq = seq(-1,1.4,by=0.2)   # for rep =23
 
for (rep in 1:nrep) 
{ 
print(rep)

data = as.matrix(read.table(paste("rep",rep,"_data_TimWood.txt",sep=""), header=T))


np           = length(data[,1])
p            = data[,1:2]
data = cbind(1:np,data[,3])

p2col        = rbind(p,Cconstredges)

Tri = as.matrix(read.table(paste("rep",rep,"_t_TimWood.txt",sep=""), header=F))

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

basisobj = create.FEM.basis(p2col, e, Tri, order)


#  set up a dummy FEM functional data object

simfd = fd(numeric(basisobj$nbasis),basisobj)




#  smooth the data without using covariates

GCVres = smooth.FEM.fd.GCV.NOdiff(data,simfd,loglambdaseq)

loglambdaoptim.FES[rep] = loglambdaseq[which.min(GCVres$GCVseq)]
GCV.FES[rep]          = min(GCVres$GCVseq)
EDF.FES[rep]          = GCVres$EDFseq[which.min(GCVres$GCVseq)]


sigmahat.FES[rep] = GCVres$sigmahatseq[which.min(GCVres$GCVseq)]


windows() 
plot(loglambdaseq, GCVres$GCVseq,main=rep)
abline(v=loglambdaoptim.FES[rep])


lambdaCovar = 10^loglambdaoptim.FES[rep] 
fit.FES.Covar = smooth.FEM.fd.NOdiff(data,simfd,lambdaCovar)


coefres=data.frame(coeffelsplobj=fit.FES.Covar$felsplobj$coef,coeflaplacefd=fit.FES.Covar$laplacefd$coef)
colnames(coefres)=c("coef.felsplobj","coef.laplacefd")

write.table(coefres,paste("rep",rep,"_results.FES.NOcovar.order1.NOdiff.txt",sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE)



XYmat.FES.Covar       = matrix(NA,nrow=length(xvec),ncol=length(yvec))
XYmat.FES.Covar.error = XYmat.FES.Covar
for (indexX in 1:length(xvec))
{
for (indexY in iY[inSide(fsb,x=rep(xvec[indexX],ny),y=yvec)])
{
XYmat.FES.Covar[indexX,indexY]      = eval.FEM.fd(xvec[indexX], yvec[indexY], fit.FES.Covar$felsplobj) 
XYmat.FES.Covar.error[indexX,indexY] = XYmat.FES.Covar[indexX,indexY] - tru[indexX,indexY]
}
}

MSEf.FES[rep] = mean(na.omit(as.vector(XYmat.FES.Covar.error^2)))
graphics.off()

}


results.FES = data.frame(loglambdaoptim=loglambdaoptim.FES, GCV=GCV.FES, EDF=EDF.FES, sigmahat=sigmahat.FES, MSEf=MSEf.FES)
write.table(results.FES,"results.FES.NOcovar.order1.NOdiff.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)






#########################################################
#########################################################
#########################################################
#########################################################


results.order2=read.table("results.FES.NOcovar.order2.txt", header=TRUE)
results.order1=read.table("results.FES.NOcovar.order1.txt", header=TRUE)
results.order1.NOdiff=read.table("results.FES.NOcovar.order1.NOdiff.txt", header=TRUE)




boxplot(sqrt(results.order2$MSEf),sqrt(results.order1$MSEf),sqrt(results.order1.NOdiff$MSEf),names=c("order2", "order1", "order1.NOdiff"),main="Root Mean Square Error f")




#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################
#########################################################



##################################
##################################
############ COVARIATES ##########
##################################
##################################


order = 2 
betahat.FES   = matrix(0,nrow=nrep,ncol=2)
sigmahat.FES  = numeric(nrep)
MSEf.FES      = numeric(nrep)
GCV.FES          = numeric(nrep)
EDF.FES          = numeric(nrep)
loglambdaoptim.FES  = numeric(nrep)


 
loglambdaseq =  seq(-1.8,0.8,by=0.2)

#loglambdaseq = seq(-1.6,0.6,by=0.2) # for rep=c(13,30,35,37,40)

#loglambdaseq = seq(-1,1.4,by=0.2)   # for rep =23
 
for (rep in 1:nrep) 
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

######windows() 
######plot(fsb[[1]]$x,fsb[[1]]$y,type="l",asp=1,main=rep)
######title("Triangulation")
######for (ne in 1:nt)
######{polygon(c(p2col[Tri[ne,1],1],p2col[Tri[ne,2],1],p2col[Tri[ne,3],1]),c(p2col[Tri[ne,1],2],p2col[Tri[ne,2],2],p2col[Tri[ne,3],2]))
######}

#  we don't need an edge matrix for these analyses

e = NULL

#  set up the FEM basis object and plot it

basisobj = create.FEM.basis(p2col, e, Tri, order)


#  set up a dummy FEM functional data object

simfd = fd(numeric(basisobj$nbasis),basisobj)




#  smooth the data using covariates

GCVres = smooth.FEM.fd.Covar.GCV(dataCovar,desmat,simfd,loglambdaseq,0.95)

loglambdaoptim.FES[rep] = loglambdaseq[which.min(GCVres$GCVseq)]
GCV.FES[rep]          = min(GCVres$GCVseq)
EDF.FES[rep]          = GCVres$EDFseq[which.min(GCVres$GCVseq)]

betahat.FES[rep,] = as.vector(GCVres$approxCIbeta[,"betahat"])

sigmahat.FES[rep] = GCVres$sigmahatseq[which.min(GCVres$GCVseq)]


windows() 
plot(loglambdaseq, GCVres$GCVseq,main=rep)
abline(v=loglambdaoptim.FES[rep])





lambdaoptim = 10^loglambdaoptim.FES[rep]
fit.FES.Covar = smooth.FEM.fd.Covar(dataCovar,desmat,simfd,lambdaoptim)


coefres=data.frame(coeffelsplobj=fit.FES.Covar$felsplobj$coef,coeflaplacefd=fit.FES.Covar$laplacefd$coef)
colnames(coefres)=c("coef.felsplobj","coef.laplacefd")

write.table(coefres,paste("rep",rep,"_results.FES.covar.order2.txt",sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE)


#plot.FEM(fit.FES.Covar$felsplobj)
#windows()
#plot.FEM(fd(coefres$coef.felsplobj,basisobj))
#windows()
#plot.FEM(fit.FES.Covar$laplacefd)
#windows()
#plot.FEM(fd(coefres$coef.laplacefd,basisobj))



XYmat.FES.Covar       = matrix(NA,nrow=length(xvec),ncol=length(yvec))
XYmat.FES.Covar.error = XYmat.FES.Covar
for (indexX in 1:length(xvec))
{
for (indexY in iY[inSide(fsb,x=rep(xvec[indexX],ny),y=yvec)])
{
XYmat.FES.Covar[indexX,indexY]      = eval.FEM.fd(xvec[indexX], yvec[indexY], fit.FES.Covar$felsplobj) 
XYmat.FES.Covar.error[indexX,indexY] = XYmat.FES.Covar[indexX,indexY] - tru[indexX,indexY]
}
}

MSEf.FES[rep] = mean(na.omit(as.vector(XYmat.FES.Covar.error^2)))
graphics.off()

}


results.FES = data.frame(loglambdaoptim=loglambdaoptim.FES, GCV=GCV.FES, EDF=EDF.FES, betahat=betahat.FES, sigmahat=sigmahat.FES, MSEf=MSEf.FES)
write.table(results.FES,"results.FES.covar.order2.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)






###########################################
###########################################
###########################################


order = 1
betahat.FES      = matrix(0,nrow=nrep,ncol=2)
sigmahat.FES     = numeric(nrep)
MSEf.FES         = numeric(nrep)
GCV.FES          = numeric(nrep)
EDF.FES          = numeric(nrep)
loglambdaoptim.FES  = numeric(nrep)



 
 
loglambdaseq = seq(-1,1,by=0.2)

#loglambdaseq = seq(-1.6,0.6,by=0.2) # for rep=c(13,30,35,37,40)

#loglambdaseq = seq(-1,1.4,by=0.2)   # for rep =23
 
for (rep in 1:nrep) 
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

####windows() 
####plot(fsb[[1]]$x,fsb[[1]]$y,type="l",asp=1,main=rep)
####title("Triangulation")
####for (ne in 1:nt)
####{polygon(c(p2col[Tri[ne,1],1],p2col[Tri[ne,2],1],p2col[Tri[ne,3],1]),c(p2col[Tri[ne,1],2],p2col[Tri[ne,2],2],p2col[Tri[ne,3],2]))
####}

#  we don't need an edge matrix for these analyses

e = NULL

#  set up the FEM basis object and plot it

basisobj = create.FEM.basis(p2col, e, Tri, order)


#  set up a dummy FEM functional data object

simfd = fd(numeric(basisobj$nbasis),basisobj)




#  smooth the data using covariates

GCVres = smooth.FEM.fd.Covar.GCV(dataCovar,desmat,simfd,loglambdaseq,0.95)

loglambdaoptim.FES[rep] = loglambdaseq[which.min(GCVres$GCVseq)]
GCV.FES[rep]          = min(GCVres$GCVseq)
EDF.FES[rep]          = GCVres$EDFseq[which.min(GCVres$GCVseq)]


betahat.FES[rep,] = as.vector(GCVres$approxCIbeta[,"betahat"])

sigmahat.FES[rep] = GCVres$sigmahatseq[which.min(GCVres$GCVseq)]
 



windows() 
plot(loglambdaseq, GCVres$GCVseq,main=rep)
abline(v=loglambdaoptim.FES[rep])





lambdaoptim = 10^loglambdaoptim.FES[rep]
fit.FES.Covar = smooth.FEM.fd.Covar(dataCovar,desmat,simfd,lambdaoptim)

coefres=data.frame(coeffelsplobj=fit.FES.Covar$felsplobj$coef)
colnames(coefres)=c("coef.felsplobj")

write.table(coefres,paste("rep",rep,"_results.FES.covar.order1.txt",sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE)


XYmat.FES.Covar       = matrix(NA,nrow=length(xvec),ncol=length(yvec))
XYmat.FES.Covar.error = XYmat.FES.Covar
for (indexX in 1:length(xvec))
{
for (indexY in iY[inSide(fsb,x=rep(xvec[indexX],ny),y=yvec)])
{
XYmat.FES.Covar[indexX,indexY]      = eval.FEM.fd(xvec[indexX], yvec[indexY], fit.FES.Covar$felsplobj) 
XYmat.FES.Covar.error[indexX,indexY] = XYmat.FES.Covar[indexX,indexY] - tru[indexX,indexY]
}
}

MSEf.FES[rep] = mean(na.omit(as.vector(XYmat.FES.Covar.error^2)))
graphics.off()

}


results.FES = data.frame(loglambdaoptim=loglambdaoptim.FES, GCV=GCV.FES, EDF=EDF.FES, betahat=betahat.FES,sigmahat=sigmahat.FES,MSEf=MSEf.FES)
write.table(results.FES,"results.FES.covar.order1.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)







###########################################
###########################################
###########################################


order = 1 # NO diff
betahat.FES   = matrix(0,nrow=nrep,ncol=2)
sigmahat.FES  = numeric(nrep)
MSEf.FES      = numeric(nrep)
GCV.FES          = numeric(nrep)
EDF.FES          = numeric(nrep)
loglambdaoptim.FES  = numeric(nrep)


 
loglambdaseq =  seq(-1.8,0.8,by=0.2)

#loglambdaseq = seq(-1.6,0.6,by=0.2) # for rep=c(13,30,35,37,40)

#loglambdaseq = seq(-1,1.4,by=0.2)   # for rep =23
 
for (rep in 1:nrep) 
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

####windows() 
####plot(fsb[[1]]$x,fsb[[1]]$y,type="l",asp=1,main=rep)
####title("Triangulation")
####for (ne in 1:nt)
####{polygon(c(p2col[Tri[ne,1],1],p2col[Tri[ne,2],1],p2col[Tri[ne,3],1]),c(p2col[Tri[ne,1],2],p2col[Tri[ne,2],2],p2col[Tri[ne,3],2]))
####}

#  we don't need an edge matrix for these analyses

e = NULL

#  set up the FEM basis object and plot it

basisobj = create.FEM.basis(p2col, e, Tri, order)


#  set up a dummy FEM functional data object

simfd = fd(numeric(basisobj$nbasis),basisobj)




#  smooth the data using covariates

GCVres = smooth.FEM.fd.Covar.GCV.NOdiff(dataCovar,desmat,simfd,loglambdaseq,0.95)

loglambdaoptim.FES[rep] = loglambdaseq[which.min(GCVres$GCVseq)]
GCV.FES[rep]          = min(GCVres$GCVseq)
EDF.FES[rep]          = GCVres$EDFseq[which.min(GCVres$GCVseq)]

betahat.FES[rep,] = as.vector(GCVres$approxCIbeta[,"betahat"])

sigmahat.FES[rep] = GCVres$sigmahatseq[which.min(GCVres$GCVseq)]



windows() 
plot(loglambdaseq, GCVres$GCVseq,main=rep)
abline(v=loglambdaoptim.FES[rep])




lambdaoptim = 10^loglambdaoptim.FES[rep]
fit.FES.Covar = smooth.FEM.fd.Covar.NOdiff(dataCovar,desmat,simfd,lambdaoptim)


coefres=data.frame(coeffelsplobj=fit.FES.Covar$felsplobj$coef,coeflaplacefd=fit.FES.Covar$laplacefd$coef)
colnames(coefres)=c("coef.felsplobj","coef.laplacefd")

write.table(coefres,paste("rep",rep,"_results.FES.covar.order1.NOdiff.txt",sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE)


XYmat.FES.Covar       = matrix(NA,nrow=length(xvec),ncol=length(yvec))
XYmat.FES.Covar.error = XYmat.FES.Covar
for (indexX in 1:length(xvec))
{
for (indexY in iY[inSide(fsb,x=rep(xvec[indexX],ny),y=yvec)])
{
XYmat.FES.Covar[indexX,indexY]      = eval.FEM.fd(xvec[indexX], yvec[indexY], fit.FES.Covar$felsplobj) 
XYmat.FES.Covar.error[indexX,indexY] = XYmat.FES.Covar[indexX,indexY] - tru[indexX,indexY]
}
}

MSEf.FES[rep] = mean(na.omit(as.vector(XYmat.FES.Covar.error^2)))
graphics.off()

}


results.FES = data.frame(loglambdaoptim=loglambdaoptim.FES, GCV=GCV.FES, EDF=EDF.FES, betahat=betahat.FES, sigmahat=sigmahat.FES, MSEf=MSEf.FES)
write.table(results.FES,"results.FES.covar.order1.NOdiff.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)


#######################
#######################
#######################
#######################









results.old=read.table("results.FES.txt", header=TRUE)

results.order2=read.table("results.FES.covar.order2.txt", header=TRUE)
results.order1=read.table("results.FES.covar.order1.txt", header=TRUE)
results.order1.NOdiff=read.table("results.FES.covar.order1.NOdiff.txt", header=TRUE)



boxplot(sqrt(results.old$MSEf),sqrt(results.order2$MSEf),sqrt(results.order1$MSEf),sqrt(results.order1.NOdiff$MSEf),names=c("order2.old","order2", "order1", "order1.NOdiff"),main="Root Mean Square Error f")

boxplot(sqrt(results.order2$MSEf),sqrt(results.order1$MSEf),sqrt(results.order1.NOdiff$MSEf),names=c("order2", "order1", "order1.NOdiff"),main="Root Mean Square Error f")




sqrt(mean((results.order2$betahat.1+0.5)^2))
sqrt(mean((results.order1$betahat.1+0.5)^2))
sqrt(mean((results.order1.NOdiff$betahat.1+0.5)^2))


sqrt(mean((results.order2$betahat.2-0.2)^2))
sqrt(mean((results.order1$betahat.2-0.2)^2))
sqrt(mean((results.order1.NOdiff$betahat.2-0.2)^2))


sqrt(mean((results.order2$sigmahat-0.5)^2))
sqrt(mean((results.order1$sigmahat-0.5)^2))
sqrt(mean((results.order1.NOdiff$sigmahat-0.5)^2))

######
######
######rep=2
######
######
######graphics.off()
######rep.results.order2         = read.table(paste("rep",rep,"_results.FES.covar.order2.txt",sep=""), header=T)
######rep.results.order1         = read.table(paste("rep",rep,"_results.FES.covar.order1.txt",sep=""), header=T)
######rep.results.order1.NOdiff  = read.table(paste("rep",rep,"_results.FES.covar.order1.NOdiff.txt",sep=""), header=T)
######
######
######
######data = as.matrix(read.table(paste("rep",rep,"_data_TimWood.txt",sep=""), header=T))
######
######
######np           = length(data[,1])
######p            = data[,1:2]
######p2col        = rbind(p,Cconstredges)
######Tri = as.matrix(read.table(paste("rep",rep,"_t_TimWood.txt",sep=""), header=F))
######nt = dim(Tri)[[1]]
######
######e = NULL
######
######
######basisobj.order2 = create.FEM.basis(p2col, e, Tri, 2)
######felsplobj.order2=fd(rep.results.order2$coef.felsplobj,basisobj.order2)
######laplacefd.order2=fd(rep.results.order2$coef.laplacefd,basisobj.order2)
######
######
######basisobj.order1 = create.FEM.basis(p2col, e, Tri, 1)
######
######felsplobj.order1=fd(rep.results.order1$coef.felsplobj,basisobj.order1)
######
######felsplobj.order1.NOdiff=fd(rep.results.order1.NOdiff$coef.felsplobj,basisobj.order1)
######laplacefd.order1.NOdiff=fd(rep.results.order1.NOdiff$coef.laplacefd,basisobj.order1)
######
######
######
######
######
######plot.FEM(felsplobj.order2)
######title("felsplobj.order2")
######
######
######windows()
######plot.FEM(felsplobj.order1)
######title("felsplobj.order1")
######
######
######windows()
######plot.FEM(felsplobj.order1.NOdiff)
######title("felsplobj.order1.NOdiff")
######
######
######
######
######windows()
######plot.FEM(laplacefd.order2)
######title("laplacefd.order2")
######
######
######windows()
######plot.FEM(laplacefd.order1.NOdiff)
######title("laplacefd.order1.NOdiff")
######
######
######
######
######
######
######
######
######
######
######
######XYmat.FES.order2       = matrix(NA,nrow=length(xvec),ncol=length(yvec))
######XYmat.FES.order2.error = XYmat.FES.order2
######
######XYmat.FES.order1       = matrix(NA,nrow=length(xvec),ncol=length(yvec))
######XYmat.FES.order1.error = XYmat.FES.order1
######
######XYmat.FES.order1.NOdiff       = matrix(NA,nrow=length(xvec),ncol=length(yvec))
######XYmat.FES.order1.NOdiff.error = XYmat.FES.order1.NOdiff
######
######
######for (indexX in 1:length(xvec))
######{
######for (indexY in iY[inSide(fsb,x=rep(xvec[indexX],ny),y=yvec)])
######{
######XYmat.FES.order2[indexX,indexY]      = eval.FEM.fd(xvec[indexX], yvec[indexY],felsplobj.order2) 
######XYmat.FES.order2.error[indexX,indexY] = XYmat.FES.order2[indexX,indexY] - tru[indexX,indexY]
######
######XYmat.FES.order1[indexX,indexY]      = eval.FEM.fd(xvec[indexX], yvec[indexY],felsplobj.order1) 
######XYmat.FES.order1.error[indexX,indexY] = XYmat.FES.order1[indexX,indexY] - tru[indexX,indexY]
######
######XYmat.FES.order1.NOdiff[indexX,indexY]      = eval.FEM.fd(xvec[indexX], yvec[indexY],felsplobj.order1.NOdiff) 
######XYmat.FES.order1.NOdiff.error[indexX,indexY] = XYmat.FES.order1.NOdiff[indexX,indexY] - tru[indexX,indexY]
######}
######}
######
######
######
######
######
######windows()
######par(mfrow=c(2,2))
######zlimits=c(min(c(XYmat.FES.order1[which.min(XYmat.FES.order1)],XYmat.FES.order2[which.min(XYmat.FES.order2)],XYmat.FES.order1.NOdiff[which.min(XYmat.FES.order1.NOdiff)])),max(c(XYmat.FES.order1[which.max(XYmat.FES.order1)],XYmat.FES.order2[which.max(XYmat.FES.order2)],XYmat.FES.order1.NOdiff[which.max(XYmat.FES.order1.NOdiff)])))
######
######image(xvec,yvec,XYmat.FES.order2,main="order2",zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
######contour(xvec,yvec,XYmat.FES.order2,zlim=zlimits,add=T)
######
######image(xvec,yvec,XYmat.FES.order1,main="order1",zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
######contour(xvec,yvec,XYmat.FES.order1,zlim=zlimits,add=T)
######
######image(xvec,yvec,XYmat.FES.order1.NOdiff,main="order1.NOdiff",zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
######contour(xvec,yvec,XYmat.FES.order1.NOdiff,zlim=zlimits,add=T)   
######
######
######
######
######XYmat.FES.order2.error.abs=abs(XYmat.FES.order2.error)
######XYmat.FES.order1.error.abs=abs(XYmat.FES.order1.error)
######XYmat.FES.order1.NOdiff.error.abs=abs(XYmat.FES.order1.NOdiff.error)
######
######windows()
######par(mfrow=c(2,2))
######zlimits=c(min(c(XYmat.FES.order1.error.abs[which.min(XYmat.FES.order1.error.abs)],XYmat.FES.order2.error.abs[which.min(XYmat.FES.order2.error.abs)],XYmat.FES.order1.NOdiff.error.abs[which.min(XYmat.FES.order1.NOdiff.error.abs)])),max(c(XYmat.FES.order1.error.abs[which.max(XYmat.FES.order1.error.abs)],XYmat.FES.order2.error.abs[which.max(XYmat.FES.order2.error.abs)],XYmat.FES.order1.NOdiff.error.abs[which.max(XYmat.FES.order1.NOdiff.error.abs)])))
######
######image(xvec,yvec,XYmat.FES.order2.error.abs,main="order2.error.abs",zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
######contour(xvec,yvec,XYmat.FES.order2.error.abs,zlim=zlimits,add=T)
######
######image(xvec,yvec,XYmat.FES.order1.error.abs,main="order1.error.abs",zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
######contour(xvec,yvec,XYmat.FES.order1.error.abs,zlim=zlimits,add=T)
######
######image(xvec,yvec,XYmat.FES.order1.NOdiff.error.abs,main="order1.NOdiff.error.abs",zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
######contour(xvec,yvec,XYmat.FES.order1.NOdiff.error.abs,zlim=zlimits,add=T)   
######
######
######
######
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################

#### Pictures for draft 2011_semiparam...


############### Plot the function, and its boundary
##############fsb <- list(fs.boundary())
##############nx<-750
##############ny<-300 
##############xvec <- seq(-1,4,length=nx)
##############yvec<-seq(-1,1,length=ny)
##############xx <- rep(xvec,ny)
##############yy<-rep(yvec,rep(nx,ny))
##############iY=1:ny
##############tru <- matrix(fs.test(xx,yy),nx,ny) ## truth
##############
##############persp3d(xvec,yvec,tru, col='red', xlab="",ylab="",zlab="", xlim=c(-1.1,3.5),alpha=1)
##############


rep=2


results.FES.order2=read.table("results.FES.covar.order2.txt", header=TRUE)
results.FES.order1=read.table("results.FES.covar.order1.NOdiff.txt", header=TRUE)




graphics.off()
rep.results.order2         = read.table(paste("rep",rep,"_results.FES.covar.order2.txt",sep=""), header=T)
rep.results.order1  = read.table(paste("rep",rep,"_results.FES.covar.order1.NOdiff.txt",sep=""), header=T)



data = as.matrix(read.table(paste("rep",rep,"_data_TimWood.txt",sep=""), header=T))


np           = length(data[,1])
p            = data[,1:2]
p2col        = rbind(p,Cconstredges)
Tri = as.matrix(read.table(paste("rep",rep,"_t_TimWood.txt",sep=""), header=F))
nt = dim(Tri)[[1]]

e = NULL


basisobj.order2 = create.FEM.basis(p2col, e, Tri, 2)
felsplobj.order2=fd(rep.results.order2$coef.felsplobj,basisobj.order2)
laplacefd.order2=fd(rep.results.order2$coef.laplacefd,basisobj.order2)


basisobj.order1 = create.FEM.basis(p2col, e, Tri, 1)
felsplobj.order1=fd(rep.results.order1$coef.felsplobj,basisobj.order1)
laplacefd.order1=fd(rep.results.order1$coef.laplacefd,basisobj.order1)





plot.FEM(felsplobj.order2)
title("felsplobj.order2")


windows()
plot.FEM(felsplobj.order1)
title("felsplobj.order1")



windows()
plot.FEM(laplacefd.order2)
title("laplacefd.order2")


windows()
plot.FEM(laplacefd.order1)
title("laplacefd.order1")











XYmat.FES.order2       = matrix(NA,nrow=length(xvec),ncol=length(yvec))
XYmat.FES.order2.error = XYmat.FES.order2

XYmat.FES.order1       = matrix(NA,nrow=length(xvec),ncol=length(yvec))
XYmat.FES.order1.error = XYmat.FES.order1


for (indexX in 1:length(xvec))
{
for (indexY in iY[inSide(fsb,x=rep(xvec[indexX],ny),y=yvec)])
{
XYmat.FES.order2[indexX,indexY]      = eval.FEM.fd(xvec[indexX], yvec[indexY],felsplobj.order2) 
XYmat.FES.order2.error[indexX,indexY] = XYmat.FES.order2[indexX,indexY] - tru[indexX,indexY]

XYmat.FES.order1[indexX,indexY]      = eval.FEM.fd(xvec[indexX], yvec[indexY],felsplobj.order1) 
XYmat.FES.order1.error[indexX,indexY] = XYmat.FES.order1[indexX,indexY] - tru[indexX,indexY]

}
}




library(fields)

dataCovar    =  cbind(1:np,data[,4])
desmat       =  data[,5:6]

fit.Tps.Covar=Tps(x=p,Y=dataCovar[,2],Z=desmat)

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



persp3d(xvec,yvec,XYmat.Tps.Covar.trend.nonparam, col='red', xlab="",ylab="",zlab="",alpha=1)




##### Krig


fit.Krig.Covar=Krig(x=p,Y=dataCovar[,2],Z=desmat)


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



**********

results.Tps=read.table("results.Tps.txt", header=TRUE)
results.Krig=read.table("results.Krig.txt", header=TRUE)
results.soap=read.table("results.soap.bis.txt", header=TRUE)

results.FES.order2=read.table("results.FES.covar.order2.txt", header=TRUE)
results.FES.order1=read.table("results.FES.covar.order1.NOdiff.txt", header=TRUE)



zlimits=c(min(c(tru[which.min(tru)],XYmat.Tps.Covar.trend.nonparam[which.min(XYmat.Tps.Covar.trend.nonparam)],XYmat.Krig.Covar.trend.nonparam[which.min(XYmat.Krig.Covar.trend.nonparam)],XYmat.soap.Covar[which.min(XYmat.soap.Covar)],XYmat.FES.order1[which.min(XYmat.FES.order1)],XYmat.FES.order2[which.min(XYmat.FES.order2)])),max(c(tru[which.max(tru)],XYmat.Tps.Covar.trend.nonparam[which.max(XYmat.Tps.Covar.trend.nonparam)],XYmat.Krig.Covar.trend.nonparam[which.max(XYmat.Krig.Covar.trend.nonparam)],XYmat.soap.Covar[which.max(XYmat.soap.Covar)],XYmat.FES.order1[which.max(XYmat.FES.order1)],XYmat.FES.order2[which.max(XYmat.FES.order2)])))

graphics.off()

par (layout(rbind(c(1,2,2), c(3,4,5)), widths=c(2,1,1), heights=c(1,1)),mex=1,pty="m", font.main=3,font.lab=1.5, font.axis=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)

image(xvec,yvec,tru,main="True surface",zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,tru,zlim=zlimits,add=T)   
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






par (mfrow=c(3,2),mar=c(6,5,2,1),mex=0.7,pty="m", font.main=1.1,font.lab=1.1, font.axis=1.1,cex.lab=1.1,cex.axis=0.9)
zlimits=c(min(c(tru[which.min(tru)],XYmat.Tps.Covar.trend.nonparam[which.min(XYmat.Tps.Covar.trend.nonparam)],XYmat.Krig.Covar.trend.nonparam[which.min(XYmat.Krig.Covar.trend.nonparam)],XYmat.soap.Covar[which.min(XYmat.soap.Covar)],XYmat.FES.order1[which.min(XYmat.FES.order1)],XYmat.FES.order2[which.min(XYmat.FES.order2)])),max(c(tru[which.max(tru)],XYmat.Tps.Covar.trend.nonparam[which.max(XYmat.Tps.Covar.trend.nonparam)],XYmat.Krig.Covar.trend.nonparam[which.max(XYmat.Krig.Covar.trend.nonparam)],XYmat.soap.Covar[which.max(XYmat.soap.Covar)],XYmat.FES.order1[which.max(XYmat.FES.order1)],XYmat.FES.order2[which.max(XYmat.FES.order2)])))

########image(xvec,yvec,tru,main="True surface",zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
########contour(xvec,yvec,tru,zlim=zlimits,add=T)
########lines(fsb[[1]]$x,fsb[[1]]$y,lwd=2)


image(xvec,yvec,XYmat.Krig.Covar.trend.nonparam,main="KRIG, surface estimate",zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.Krig.Covar.trend.nonparam,zlim=zlimits,add=T)
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=2)


image(xvec,yvec,XYmat.Tps.Covar.trend.nonparam,main="TPS, surface estimate",zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.Tps.Covar.trend.nonparam,zlim=zlimits,add=T)
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=2)

image(xvec,yvec,XYmat.soap.Covar,main="SOAP, surface estimate",zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.soap.Covar,zlim=zlimits,add=T)
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=2)

image(xvec,yvec,XYmat.FES.order1,main="SSR1, surface estimate",zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.FES.order1,zlim=zlimits,add=T)   
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=2)

##########
##########image(xvec,yvec,XYmat.FES.order2,main="FELspline quadratic, surface estimate",zlim=zlimits,xlab="",ylab="",col=heat.colors(100),asp=1)
##########contour(xvec,yvec,XYmat.FES.order2,zlim=zlimits,add=T)   
##########lines(fsb[[1]]$x,fsb[[1]]$y,lwd=2)

#dev.copy(device=postscript,file="fig_TimWood_4surface-estimates.ps",horizontal=FALSE,width=7, height=7)
#dev.off()






graphics.off()

par (mfrow=c(3,2),mar=c(6,5,2,1),mex=1,pty="m", font.main=1.5,font.lab=1.5, font.axis=1.1,cex.lab=1.5,cex.axis=1)

XYmat.Tps.Covar.error.abs=abs(XYmat.Tps.Covar.error)
XYmat.Krig.Covar.error.abs=abs(XYmat.Krig.Covar.error)
XYmat.soap.Covar.error.abs=abs(XYmat.soap.Covar.error)
XYmat.FES.order1.error.abs=abs(XYmat.FES.order1.error)
XYmat.FES.order2.error.abs=abs(XYmat.FES.order2.error)


zlimitserr=c(min(c(XYmat.Tps.Covar.error.abs[which.min(XYmat.Tps.Covar.error.abs)],XYmat.Krig.Covar.error.abs[which.min(XYmat.Krig.Covar.error.abs)],XYmat.soap.Covar.error.abs[which.min(XYmat.soap.Covar.error.abs)],XYmat.FES.order1.error.abs[which.min(XYmat.FES.order1.error.abs)],XYmat.FES.order2.error.abs[which.min(XYmat.FES.order2.error.abs)])),max(c(XYmat.Tps.Covar.error.abs[which.max(XYmat.Tps.Covar.error.abs)],XYmat.Krig.Covar.error.abs[which.max(XYmat.Krig.Covar.error.abs)],XYmat.soap.Covar.error.abs[which.max(XYmat.soap.Covar.error.abs)],XYmat.FES.order1.error.abs[which.max(XYmat.FES.order1.error.abs)],XYmat.FES.order2.error.abs[which.max(XYmat.FES.order2.error.abs)])))

######plot("n",frame.plot=F,axes = F,pty="n")

image(xvec,yvec,XYmat.Krig.Covar.error.abs,main="KRIG, absolute residuals",zlim=zlimitserr,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.Krig.Covar.error.abs,zlim=zlimitserr,add=T)

image(xvec,yvec,XYmat.Tps.Covar.error.abs,main="TPS, absolute residuals",zlim=zlimitserr,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.Tps.Covar.error.abs,zlim=zlimitserr,add=T)

image(xvec,yvec,XYmat.soap.Covar.error.abs,main="SOAP, absolute residuals",zlim=zlimitserr,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.soap.Covar.error.abs,zlim=zlimitserr,add=T)

image(xvec,yvec,XYmat.FES.order1.error.abs,main="SSR1, absolute residuals",zlim=zlimitserr,xlab="",ylab="",col=heat.colors(100),asp=1)
contour(xvec,yvec,XYmat.FES.order1.error.abs,zlim=zlimitserr,add=T)   


######image(xvec,yvec,XYmat.FES.order2.error.abs,main="FELspline quadratic, absolute residuals",zlim=zlimitserr,xlab="",ylab="",col=heat.colors(100),asp=1)
######contour(xvec,yvec,XYmat.FES.order2.error.abs,zlim=zlimitserr,add=T)   


dev.copy(device=postscript,file="fig_TimWood_4surface-abs-error.ps",horizontal=FALSE,width=7, height=7)
dev.off()


####dev.copy(device=postscript,file="fig_TimWood_surface-abs-error.ps",horizontal=FALSE,width=7, height=7)
####dev.off()



par (mfrow=c(2,2),mar=c(6,5,2,1),mex=1,pty="m", font.main=1.5,font.lab=1.2, font.axis=1.1,cex.lab=1.2,cex.axis=1.1)

boxplot(sqrt(results.Krig$MSEf),sqrt(results.Tps$MSEf),sqrt(results.soap$MSEf),sqrt(results.FES.order1$MSEf),sqrt(results.FES.order2$MSEf),names=c("KRIG","TPS","SOAP", "SSR1", "SSR2"), ylab="RMSE surface estimate",las=3)
##
##dev.copy(device=postscript,file="fig_TimWood_surface-RMSE.ps",horizontal=FALSE,width=7, height=7)
##dev.off()



###############################
########## SLIDES #############

par (mfrow=c(2,2),mar=c(6,5,2,1),mex=1,pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.5,cex.axis=1.5,cex.main=2)
boxplot(sqrt(results.Tps$MSEf),sqrt(results.FES.order2$MSEf),names=c("TPS", "SSR"), main="RMSE f")
boxplot(sqrt(results.Krig$MSEf),sqrt(results.Tps$MSEf),sqrt(results.soap$MSEf),sqrt(results.FES.order1$MSEf),sqrt(results.FES.order2$MSEf),names=c("Krig","TPS","Soap", "SSR1", "SSR2"),las=3, main="RMSE f")

dev.copy(device=postscript,file="fig_TimWood_surface-RMSE_slides.ps",horizontal=FALSE,width=7, height=7)
dev.off()


###############################





graphics.off()

#par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.8,pty="m", font.main=1.5,font.lab=1.5, font.axis=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)

par (layout(rbind(c(1,2), c(3,4)), widths=c(1,1), heights=c(1,1)),mar=c(6,5,2,1),mex=0.9,pty="m", font.main=3,font.lab=1.3, font.axis=1.3,cex.lab=1.3,cex.axis=1.3,cex.main=1.4)


zlimits2=c(min(c(tru[which.min(tru)],XYmat.Tps.Covar.trend.nonparam[which.min(XYmat.Tps.Covar.trend.nonparam)],XYmat.FES.order2[which.min(XYmat.FES.order2)])),max(c(tru[which.max(tru)],XYmat.Tps.Covar.trend.nonparam[which.max(XYmat.Tps.Covar.trend.nonparam)],XYmat.FES.order2[which.max(XYmat.FES.order2)])))

image(xvec,yvec,XYmat.Tps.Covar.trend.nonparam,main="TPS surface estimate",zlim=zlimits2,xlab="",ylab="",col=heat.colors(100),xlim=c(-1.2,3.7),asp=1)
contour(xvec,yvec,XYmat.Tps.Covar.trend.nonparam,zlim=zlimits2,add=T)
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=2)

image(xvec,yvec,XYmat.FES.order2,main="SSR surface estimate",zlim=zlimits2,xlab="",ylab="",col=heat.colors(100),xlim=c(-1.2,3.7),asp=1)
contour(xvec,yvec,XYmat.FES.order2,zlim=zlimits2,add=T)   
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=2)


zlimits2e=c(min(c(XYmat.Tps.Covar.error.abs[which.min(XYmat.Tps.Covar.error.abs)],XYmat.FES.order2.error.abs[which.min(XYmat.FES.order2.error.abs)])),max(c(XYmat.Tps.Covar.error.abs[which.max(XYmat.Tps.Covar.error.abs)],XYmat.FES.order2.error.abs[which.max(XYmat.FES.order2.error.abs)])))

image(xvec,yvec,XYmat.Tps.Covar.error.abs,main="TPS absolute residuals",zlim=zlimits2e,xlab="",ylab="",col=heat.colors(100),xlim=c(-1.2,3.7),asp=1)
contour(xvec,yvec,XYmat.Tps.Covar.error.abs,zlim=zlimits2e,add=T)

image(xvec,yvec,XYmat.FES.order2.error.abs,main="SSR absolute residuals",zlim=zlimits2e,xlab="",ylab="",col=heat.colors(100),xlim=c(-1.2,3.7),asp=1)
contour(xvec,yvec,XYmat.FES.order2.error.abs,zlim=zlimits2e,add=T)   


####dev.copy(device=postscript,file="fig_poster_TimWood_results.ps",horizontal=FALSE,width=7, height=5)
####dev.off()





windows()
par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.9,pty="m", font.main=1.5,font.lab=1.5, font.axis=1.3,cex.lab=1.5,cex.axis=1.5)
boxplot(sqrt(results.Tps$MSEf),sqrt(results.FES.order2$MSEf),names=c("TPS", "SSR"), ylab="RMSE f")


####dev.copy(device=postscript,file="fig_poster_TimWood_boxplots.ps",horizontal=FALSE,width=7, height=6)
####dev.off()




#par (layout(rbind(c(1,2,2), c(3,4,5)), widths=c(2,1,1), heights=c(1,1)),mar=c(6,5,2,1),mex=0.9,pty="m", font.main=3,font.lab=1.5, font.axis=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)

par (layout(rbind(c(1,3,5), c(2,4,6)), widths=c(1,1/2,1/2), heights=c(1,1)),mar=c(6,5,2,1),mex=0.9,pty="m", font.main=3,font.lab=1.5, font.axis=1.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)



image(xvec,yvec,tru,main="True surface",zlim=zlimits2,xlab="",ylab="",col=heat.colors(100),xlim=c(-1.2,3.7),asp=1)
contour(xvec,yvec,tru,zlim=zlimits2,add=T)   
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=2)



plot(fsb[[1]]$x,fsb[[1]]$y,xlab="",ylab="",type="l",asp=1,lwd=2,xlim=c(-1.1,3.6))
title("Sampled data")
for (k in 1:np)
{points(p[k,1],p[k,2],pch=19,cex=(dataCovar[k,2]-min(dataCovar[,2]))/(max(dataCovar[,2])-min(dataCovar[,2])))
}




###########plot(fsb[[1]]$x,fsb[[1]]$y,xlab="x",ylab="y",type="l",col="red",asp=1,lwd=2)
##########plot(fsb[[1]]$x,fsb[[1]]$y,xlab="",ylab="",type="n",col="red",asp=1,xlim=c(-0.85,3.85))
##########title("Triangulation")
##########for (ne in 1:nt)
##########{polygon(c(p2col[Tri[ne,1],1],p2col[Tri[ne,2],1],p2col[Tri[ne,3],1]),c(p2col[Tri[ne,1],2],p2col[Tri[ne,2],2],p2col[Tri[ne,3],2]))
##########}


plot(desmat[,1],dataCovar[,2],ylab="response variable", xlab="covariate 1")
plot(desmat[,2],dataCovar[,2],ylab="response variable", xlab="covariate 2")

boxplot(sqrt(results.Tps$MSEf),sqrt(results.FES.order2$MSEf),names=c("TPS", "SSR"), main="RMSE f")


####dev.copy(device=postscript,file="fig_poster_TimWood_data.ps",horizontal=FALSE,width=7, height=5)
####dev.off()



dev.copy(device=postscript,file="fig_SCo2011_TimWood_data.ps",horizontal=FALSE,width=7, height=5)
dev.off()






sqrt(mean((results.Krig$betahat.1+0.5)^2))
sqrt(mean((results.Tps$betahat.1+0.5)^2))
sqrt(mean((results.soap$betahat.1+0.5)^2))
sqrt(mean((results.FES.order1$betahat.1+0.5)^2))
sqrt(mean((results.FES.order2$betahat.1+0.5)^2))


sqrt(mean((results.Krig$betahat.2-0.2)^2))
sqrt(mean((results.Tps$betahat.2-0.2)^2))
sqrt(mean((results.soap$betahat.2-0.2)^2))
sqrt(mean((results.FES.order2$betahat.2-0.2)^2))
sqrt(mean((results.FES.order1$betahat.2-0.2)^2))

sqrt(mean((results.Krig$sigmahat-0.5)^2))
sqrt(mean((results.Tps$sigmahat-0.5)^2))
sqrt(mean((results.soap$sigmahat-0.5)^2))
sqrt(mean((results.FES.order2$sigmahat-0.5)^2))
sqrt(mean((results.FES.order1$sigmahat-0.5)^2))







#######################
#######################
#######################
#######################

# soap-film smoothing

library(soap)

 
# create some internal knots for soap film
knots <- data.frame(x=rep(seq(-.5,3,by=.5),4),
                    y=rep(c(-.6,-.3,.3,.6),rep(8,4)))
 
 
 

betahat.soap   = matrix(0,nrow=nrep,ncol=2)
sigmahat.soap  = numeric(nrep)
MSEf.soap      = numeric(nrep)


 
for (rep in 1:nrep) 
{ 
print(rep)

data = as.matrix(read.table(paste("rep",rep,"_data_TimWood.txt",sep=""), header=T))


np           = length(data[,1])
p            = data[,1:2]
dataNOCovar  =  data[,3]
dataCovar    =  cbind(1:np,data[,4])
desmat       =  data[,5:6]




### Soap film smoother


x=p[,1]
y=p[,2]


###########fit.soap.Covar <- gam(dataCovar[,2]~desmat+s(x,y,k=40,bs="so",xt=list(bnd=fsb)),knots=knots)
###########
###########sob <- smooth.construct2(s(x,y,bs="so",k=40,xt=list(bnd=fsb)),data=data.frame(x=x,y=y),knots=knots)
###########
###########matpred=Predict.matrix2(sob,data=list(x=xx,y=yy))
###########
###########XYmat.soap.Covar=matrix((matpred%*%fit.soap.Covar$coefficients[-c(2,3)]),nx,ny)

Covariate1=desmat[,1]

Covariate2=desmat[,2]
fit.soap.Covar <- gam(dataCovar[,2]~Covariate1+Covariate2+s(x,y,k=40,bs="so",xt=list(bnd=fsb)),knots=knots)

betahat.soap[rep,] = fit.soap.Covar$coefficients[c(2,3)]
sigmahat.soap[rep] = sqrt(fit.soap.Covar$sig2)

fit.soap.Covar.predict <- predict(fit.soap.Covar,newdata=data.frame(Covariate1=numeric(length(xx)),Covariate2=numeric(length(xx)),x=xx,y=yy),block.size=-1)

XYmat.soap.Covar.error= matrix(fit.soap.Covar.predict,nx,ny)- tru

MSEf.soap[rep] = mean(na.omit(as.vector(XYmat.soap.Covar.error^2)))

}


results.soap = data.frame(betahat=betahat.soap,sigmahat=sigmahat.soap,MSEf=MSEf.soap)


write.table(results.soap,"results.soap.bis.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)












results.FES.order2=read.table("results.FES.covar.order2.txt", header=TRUE)
results.FES.order1=read.table("results.FES.covar.order1.NOdiff.txt", header=TRUE)
results.FES.NOCovar.order2=read.table("results.FES.NOcovar.order2.txt", header=TRUE)
results.FES.NOCovar.order1=read.table("results.FES.NOcovar.order1.NOdiff.txt", header=TRUE)



boxplot(sqrt(results.FES.NOCovar.order1$MSEf),sqrt(results.FES.NOCovar.order2$MSEf),sqrt(results.FES.order1$MSEf),sqrt(results.FES.order2$MSEf),names=c("order1 NOcovar","order2 NOcovar", "FELspline lin.", "FELspline quad."), main="RMSE surface estimate",las=3)
