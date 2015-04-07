########################################
#  Data on Isle of Montreal
########################################

#  Last modified 24 February by Laura Sangalli


setwd("E:/Dati/canada/fdaR")  
source("SFDA_AllFunctions.R")
library(fda)
library(rgl)

Montreal_income=as.matrix(read.table('Montreal_income.txt', header=FALSE))
Montreal_popden=as.matrix(read.table('Montreal_popden.txt', header=FALSE))
Montreal_incden=as.matrix(read.table('Montreal_incden.txt', header=FALSE))


t = as.matrix(read.table('Montreal_Tri.txt', header=F))

#  define the vertices as points where there are data

p = as.matrix(read.table('Montreal_p.txt', header=F))

np = dim(Montreal_income)[[1]]
 
 
#  The triangle index matrix

nt = dim(t)[[1]]




plot(sqrt(Montreal_income),sqrt(Montreal_popden))


data = matrix(0,nrow=np,ncol=2)
data[,1] = 1:np 
data[,2] = sqrt(Montreal_popden)


desmat = matrix(1,nrow=np,ncol=1)
desmat[,1] =  sqrt(Montreal_income)



e = NULL

order=1


#  set up the FEM basis object and plot it

basisobj = create.FEM.basis(p, e, t, order)



#  set up a dummy FEM functional data object

Montrealfd = fd(numeric(basisobj$nbasis),basisobj)




#  smooth the data using covariates

lambda = 0.1
Montrealfelsplobj = smooth.FEM.fd.Covar(data,desmat,Montrealfd,lambda)

#  plot smoothing surface computed using covariates
windows()
plot.FEM(Montrealfelsplobj$felsplobj)



fhat=Montrealfelsplobj$felsplobj$coef
betahat= ( solve( t(desmat) %*% desmat ) ) %*% t(desmat) %*% (data[,2]-fhat[1:np,])
betahat
