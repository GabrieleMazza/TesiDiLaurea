#  Smoothing of the Meuse data in Bivand, et al, 
#  
#  Last modified on 4 February 2010 by laura

#setwd("E:/Dati/canada/Spatial Data Analysis")  


setwd("D:/Dati/canada/fdaR")  
source("2013_SSR_AllFunctions.R")
library(fda)
library(rgl)

#  The Meuse data are displayed in R, where they are included 
#  with the geoR package, loaded by command data(meuse)

#  The variable labels in the R file are:

#  row  x  y cadmium copper lead zinc   elev   dist   om ffreq soil

#  Here we shall directly read the .txt also used in the analyses carried out in Matlab 

MeuseData=read.table('MeuseData.txt', header=FALSE)

#  define the vertices as points where there are data

p = MeuseData[,2:3]
np = dim(p)[[1]]


#  plot these points

plot(p[,1],p[,2])

# #  use Delaunay triangulation to define computational mesh
# #  this mesh will be convex, which not appropriate for these data
# 
# t = delaunay(p(:,1),p(:,2));
# 
# disp(['Number of unedited triangles = ',num2str(size(t,1))])
# 
# #  plot of unedited mesh
# 
# figure(2)
# triplot(t, p(:,1),p(:,2));
# hold on
# for i=1:np
#     text(p(i,1),p(i,2),num2str(i))
# end
# hold off
# xlabel('\fontsize{13} X')
# ylabel('\fontsize{13} Y')

#  mesh is editted to remove triangles inside convex hull that
#  don't include data

MeuseTri = as.matrix(read.table('MeuseTri.txt', header=F))



#  The triangle index matrix

t = MeuseTri
nt = dim(t)[[1]]

plot(p[,1],p[,2],type="n")
polygon(c(p[MeuseTri[1,1],1],p[MeuseTri[1,2],1],p[MeuseTri[1,3],1]),c(p[MeuseTri[1,1],2],p[MeuseTri[1,2],2],p[MeuseTri[1,3],2]))
for (ne in 1:nt)
{polygon(c(p[MeuseTri[ne,1],1],p[MeuseTri[ne,2],1],p[MeuseTri[ne,3],1]),c(p[MeuseTri[ne,1],2],p[MeuseTri[ne,2],2],p[MeuseTri[ne,3],2]))
}

print(paste("Number of edited triangles = ",nt))

#  we don't need an edge matrix for these analyses

e = NULL

#  set up the FEM basis object and plot it

order=2
#order=1

basisobj = create.FEM.basis(p, e, t, order)



#### plot(basisobj)  # ATTENZIONE!! yet to be modified!

#########
######### PENSO CHE SIA INUTILE:
##########  set up the nodes, consisting of vertices plus edge midpoints
#########mknodes= makenodes(p,t,order)
#########nodes=mknodes$nodes
#########nodeindex=mknodes$nodeindex
#########nNodes = dim(nodes)[[1]]




#  set up a dummy FEM functional data object

Meusefd = fd(numeric(basisobj$nbasis),basisobj)

#  set up the data array using zinc concentration

data = matrix(0,nrow=np,ncol=2)
data[,1] = 1:np
######## CAMBIATO SOLO PER OGGI (4 Febbr 2011):
########data[,2] = log(MeuseData[,7])
data[,2] = MeuseData[,7]


#  use as covariate the distance from river; for the moment we use as regressors
#  sqrt(distance) and 
#  elev

desmat = matrix(1,nrow=np,ncol=2)
desmat[,1] = sqrt(MeuseData[,9])
desmat[,2] = (MeuseData[,8])

# Projection matrix
# ATTENZIONE: is it possible to get the projection matrix without having to compute it as below here?
#
#H= desmat %*% ( solve( t(desmat) %*% desmat ) ) %*% t(desmat)
#
## check
#library(geoR)
#data(meuse)
#dist2=(meuse$dist)^2
#zn.lm <- lm(log(zinc) ~ 0 + sqrt(dist)+dist2, meuse)
#summary(zn.lm)
#znbis.lm <- lm(data[,2] ~ 0 + desmat)
#max(abs(zn.lm$fitted- H %*% data[,2]))
#max(abs(zn.lm$coeff - ((( solve( t(desmat) %*% desmat ) ) %*% t(desmat)) %*% data[,2])))



####### CHECKS
#######
#######ZincMeusefd = smooth.FEM.fd(data,Meusefd,lambda)
#######u=ZincMeusefd$felsplobj$coef
#######
#######
#######bigH=matrix(0,nNodes,nNodes)
#######bigH[1:np,1:np]=diag(1,np)-H
#######
########bigH is simmetric and idempotent
#######
#######b = matrix(numeric(nNodes),ncol=1)
#######b[data[,1],] = data[,2]
#######
#######
#######t(u) %*% bigH %*% b 
#######
#######i=1
#######eq=((u[i]- (H[i,]%*%u[1:np,]))*(data[i,2]- (H[i,]%*%data[1:np,2])))
#######for (i in 2:np)
#######{
#######eq=eq+((u[i]- (H[i,]%*%u[1:np,]))*(data[i,2]- (H[i,]%*%data[1:np,2])))
#######}
#######eq
#######
#######f=ZincMeusefd$laplacefd$coef
#######
#######t(u) %*% bigH  %*% f 
#######
#######i=1
#######eq2=((u[i]- (H[i,]%*%u[1:np,]))*(f[i]- (H[i,]%*%f[1:np,])))
#######for (i in 2:np)
#######{
#######eq2=eq2+((u[i]- (H[i,]%*%u[1:np,]))*(f[i]- (H[i,]%*%f[1:np,])))
#######}
#######eq2




#  smooth the data using covariates

lambdaCovar = 10^(3.5)
ZincMeusefdCovar = smooth.FEM.fd.Covar(data,desmat,Meusefd,lambdaCovar)



fhat=ZincMeusefdCovar$felsplobj$coef
betahat= ( solve( t(desmat) %*% desmat ) ) %*% t(desmat) %*% (data[,2]-fhat[1:np,])


#  plot smoothing surface computed using covariates

plot.FEM(ZincMeusefdCovar$felsplobj)





#  smooth and plot data without using covariates (NO Covar)

######## CAMBIATO SOLO PER OGGI (4 Febbr 2011):
######## lambda = 10^(3)

lambda = 1
ZincMeusefd = smooth.FEM.fd(data,Meusefd,lambda)

plot.FEM(ZincMeusefd$felsplobj)


#  evaluate the solution on the points

MeuseDataFitCovar = eval.FEM.fd(p[,1], p[,2], ZincMeusefdCovar$felsplobj) +  desmat %*% betahat



#  evaluate the solution on the points NO Covar

MeuseDataFit = eval.FEM.fd(p[,1], p[,2], ZincMeusefd$felsplobj)



#  compute the residuals

MeuseResCovar = cbind(data[,1], data[,2] - MeuseDataFitCovar)
min(MeuseResCovar[,2])
max(MeuseResCovar[,2])


#  compute the residuals NO Covar

MeuseRes = cbind(data[,1], data[,2] - MeuseDataFit)
min(MeuseRes[,2])
max(MeuseRes[,2])




#  smooth the residuals

# ATTENZIONE: with what value of lambda?

ZincResMeusefdCovar = smooth.FEM.fd(MeuseResCovar,Meusefd,lambda)

windows()
plot.FEM(ZincResMeusefdCovar$felsplobj)
#title('\fontsize{16} Residuals for log10 of zinc concentration')



#  smooth the residuals NO Covar

ZincResMeusefd = smooth.FEM.fd(MeuseRes,Meusefd,lambda)

windows()
plot.FEM(ZincResMeusefd$felsplobj)
#title('\fontsize{16} Residuals for log10 of zinc concentration')






#  analysis of residuals

#  compute the squared distances and the residual products

dvec.FES = numeric(np*(np-1)/2)
rvec.FES = dvec.FES

m = 0
for (i in 2:np)
 {   
   for (j in 1:(i-1))
     {  m = m + 1
        dvec.FES[m] = sum((p[i,]-p[j,])^2)
        rvec.FES[m] = MeuseResCovar[i,2] * MeuseResCovar[j,2]
     }
 }


windows()
plot(log10(dvec.FES), rvec.FES,xlab="log_{10} Squared distance")
abline(h=0,col=2,lwd=2,lty="dashed")



#  analysis of residuals NO Covar

#  compute the squared distances and the residual products

dvec = numeric(np*(np-1)/2)
rvec = dvec

m = 0
for (i in 2:np)
 {   
   for (j in 1:(i-1))
     {  m = m + 1
        dvec[m] = sum((p[i,]-p[j,])^2)
        rvec[m] = MeuseRes[i,2] * MeuseRes[j,2]
     }
 }


windows()
plot(log10(dvec), rvec,xlab="log_{10} Squared distance")
abline(h=0,col=2,lwd=2,lty="dashed")













#  Evaluation:  lambda = 1e5 and 1e4 are clearly over-smoothing ... a
#               lot of the trend in the data is left in the residuals.
#               lambda = 1e3 has residuals with little trend, and some
#               isolated peaks and valleys. Range +-0.25
#               lambda = 1e1 has highly localized peaks, and seems
#               undersmoothed,  Range +-0.03

#  compute cross-validated error sum of squares

loglam = seq(1,5,by=0.5)
nloglam = length(loglam)
SEsaveCovar = numeric(nloglam)
for (ilam in 1:nloglam)
 {
    print(c(ilam,loglam[ilam]))
    lami = 10^loglam[ilam]
    ressave = numeric(np)
    for (idrop in 1:np)
     {
        datai = data[-idrop,]
        desmati = desmat[-idrop,]
        Meusefdobji   = smooth.FEM.fd.Covar(datai,desmati,Meusefd,lami)
        fhati = Meusefdobji$felsplobj$coef
        betahati= ( solve( t(desmati) %*% desmati ) ) %*% t(desmati) %*% (datai[,2]-fhati[1:(np-1),])
        MeuseDataFiti = eval.FEM.fd(p[idrop,1], p[idrop,2], Meusefdobji$felsplobj) + desmat[idrop,] %*% betahati
        ressave[idrop] = data[idrop,2] - MeuseDataFiti
     }
    SEsaveCovar[ilam] = sqrt(mean(ressave^2))
}

print(rbind(loglam, SEsaveCovar))

plot(loglam,SEsaveCovar)
points(loglam,SEsaveCovar,type="l")




#  compute cross-validated error sum of squares NO Covar

loglam = seq(1,5,by=0.5)
nloglam = length(loglam)
SEsave = numeric(nloglam)
for (ilam in 1:nloglam)
 {
    print(c(ilam,loglam[ilam]))
    lami = 10^loglam[ilam]
    fullind = 1:np
    ressave = numeric(np)
    for (idrop in 1:np)
     {
        datai = data[-idrop,]
        Meusefdobji   = smooth.FEM.fd(datai,Meusefd,lami)
        MeuseDataFiti = eval.FEM.fd(p[idrop,1], p[idrop,2], Meusefdobji$felsplobj)
        ressave[idrop] = data[idrop,2] - MeuseDataFiti
     }
    SEsave[ilam] = sqrt(mean(ressave^2))
}

print(rbind(loglam, SEsave))

plot(loglam,SEsave)
points(loglam,SEsave,type="l")














########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################
########################################


########################################
#  !!!!!!! HERE we use the full Meuse bank area
########################################


setwd("E:/Dati/canada/fdaR")  
source("SFDA_AllFunctions.R")
library(fda)
library(rgl)

MeuseData=read.table('MeuseData.txt', header=FALSE)

MeuseTri = as.matrix(read.table('MeuseTri_constr.txt', header=F))

#  define the vertices as points where there are data

p = as.matrix(read.table('pMeuse_constr.txt', header=F))

np = dim(MeuseData)[[1]]
 
p_3col = p
p_2col = p[,2:3]

#  The triangle index matrix

t = MeuseTri
nt = dim(t)[[1]]

data = matrix(0,nrow=np,ncol=2)
data[,1] = 1:np
data[,2] = log(MeuseData[,7])


#  use as covariate the distance from river; for the moment we use as regressors
#  sqrt(distance) and 
#  elev

desmat = matrix(1,nrow=np,ncol=2)
desmat[,1] = sqrt(MeuseData[,9])
desmat[,2] = (MeuseData[,8])


library(sp)
data(meuse.grid)

plot(p[,1],p[,2],type="n",ylim=c(min(meuse.grid[,2]),max(meuse.grid[,2])),xlim=c(min(meuse.grid[,1]),max(meuse.grid[,1])))
polygon(c(p_2col[MeuseTri[1,1],1],p_2col[MeuseTri[1,2],1],p_2col[MeuseTri[1,3],1]),c(p_2col[MeuseTri[1,1],2],p_2col[MeuseTri[1,2],2],p_2col[MeuseTri[1,3],2]))
for (ne in 1:nt)
{polygon(c(p_2col[MeuseTri[ne,1],1],p_2col[MeuseTri[ne,2],1],p_2col[MeuseTri[ne,3],1]),c(p_2col[MeuseTri[ne,1],2],p_2col[MeuseTri[ne,2],2],p_2col[MeuseTri[ne,3],2]))
}

print(paste("Number of edited triangles = ",nt))

#  we don't need an edge matrix for these analyses

e = NULL

order=1


#  set up the FEM basis object and plot it

basisobj = create.FEM.basis(p_2col, e, t, order)



#  set up a dummy FEM functional data object

Meusefd = fd(numeric(basisobj$nbasis),basisobj)




#  smooth the data using covariates

lambdaCovar = 10^(3)
ZincMeusefdCovar = smooth.FEM.fd.Covar(data,desmat,Meusefd,lambdaCovar)



fhat=ZincMeusefdCovar$felsplobj$coef
betahat= ( solve( t(desmat) %*% desmat ) ) %*% t(desmat) %*% (data[,2]-fhat[1:np,])


#  plot smoothing surface computed using covariates
windows()
plot.FEM(ZincMeusefdCovar$felsplobj)




#  evaluate the solution on the points

MeuseDataFitCovar = eval.FEM.fd(p[,1], p[,2], ZincMeusefdCovar$felsplobj) +  desmat %*% betahat



#  compute the residuals

MeuseResCovar = cbind(data[,1], data[,2] - MeuseDataFitCovar)
min(MeuseResCovar[,2])
max(MeuseResCovar[,2])


#  smooth the residuals

# ATTENZIONE: with what value of lambda?
lambda = 10^(3)
ZincResMeusefdCovar = smooth.FEM.fd(MeuseResCovar,Meusefd,lambda)

windows()
plot.FEM(ZincResMeusefdCovar$felsplobj)
title('Residuals for log10 of zinc concentration')




#  analysis of residuals

#  compute the squared distances and the residual products

dvec.FES = numeric(np*(np-1)/2)
rvec.FES = dvec.FES

m = 0
for (i in 2:np)
 {   
   for (j in 1:(i-1))
     {  m = m + 1
        dvec.FES[m] = sum((p[i,]-p[j,])^2)
        rvec.FES[m] = MeuseResCovar[i,2] * MeuseResCovar[j,2]
     }
 }


windows()
plot(log10(dvec.FES), rvec.FES,xlab="log_{10} Squared distance")
abline(h=0,col=2,lwd=2,lty="dashed")



#  Evaluation:  lambda = 1e5 and 1e4 are clearly over-smoothing ... a
#               lot of the trend in the data is left in the residuals.
#               lambda = 1e3 has residuals with little trend, and some
#               isolated peaks and valleys. Range +-0.25
#               lambda = 1e1 has highly localized peaks, and seems
#               undersmoothed,  Range +-0.03

#  compute cross-validated error sum of squares

loglam = seq(2.6,3.4,by=0.1)
nloglam = length(loglam)
SEsaveCovar = numeric(nloglam)
for (ilam in 1:nloglam)
 {
    print(c(ilam,loglam[ilam]))
    lami = 10^loglam[ilam]
    ressave = numeric(np)
    for (idrop in 1:np)
     {
        datai = data[-idrop,]
        desmati = desmat[-idrop,]
        Meusefdobji   = smooth.FEM.fd.Covar(datai,desmati,Meusefd,lami)
        fhati = Meusefdobji$felsplobj$coef
        betahati= ( solve( t(desmati) %*% desmati ) ) %*% t(desmati) %*% (datai[,2]-fhati[1:(np-1),])
        MeuseDataFiti = eval.FEM.fd(p[idrop,1], p[idrop,2], Meusefdobji$felsplobj) + desmat[idrop,] %*% betahati
        ressave[idrop] = data[idrop,2] - MeuseDataFiti
     }
    SEsaveCovar[ilam] = sqrt(mean(ressave^2))
}

print(rbind(loglam, SEsaveCovar))

plot(loglam,SEsaveCovar)
points(loglam,SEsaveCovar,type="l")










#################################
#################################
#################################
#################################
#################################
#################################
#################################
####### Smoothing kriging #######

library(fields)
library(sp)

trendorder=3

fitmeuse=Krig(x=p,Y=data[,2],Z=desmat,m=trendorder)

surface( fitmeuse, type="C")



#############

errori=data[,2]-fitmeuse$fitted.values 
errori-fitmeuse$residuals


predict(fitmeuse)-fitmeuse$fitted.values

(cbind(fields.mkpoly(p, m = trendorder),desmat)%*%t(t(fitmeuse$d)))-predict(fitmeuse,  just.fixed=T)

fields.mkpoly(p, m = trendorder)%*%t(t(fitmeuse$d[1:6]))-predict(fitmeuse,  drop.Z=T, just.fixed=TRUE)
 


predict(fitmeuse,  drop.Z=T)+fitmeuse$d[(length(fitmeuse$d)-1)]*desmat[,1]+ fitmeuse$d[length(fitmeuse$d)]*desmat[,2] -fitmeuse$fitted.values
 


#  compute the squared distances and the residual products

dvec.Krig = numeric(np*(np-1)/2)
rvec.Krig = dvec.Krig


mm = 0
for (i in 2:np)
 {   
   for (j in 1:(i-1))
     {  mm = mm + 1
        dvec.Krig[mm] = sum((p[i,]-p[j,])^2)
        rvec.Krig[mm] = fitmeuse$residuals[i] * fitmeuse$residuals[j]
     }
 }


windows()
plot(log10(dvec.Krig), rvec.Krig,xlab="log_{10} Squared distance")
abline(h=0,col=2,lwd=2,lty="dashed")





np=155
ressave=numeric(np)
    for (idrop in 1:np)
     {  
        print(idrop)
        modelp=p[-idrop,]
        modeldata = data[-idrop,]
        modeldesmat = desmat[-idrop,]
        modelfitmeuse <- Krig(x=modelp,Y=modeldata[,2],Z=modeldesmat,m=trendorder)
        validp=p[idrop,]
        validdesmat=desmat[idrop,]
        ressave[idrop] <- data[idrop,2] - predict(modelfitmeuse, validp, Z=t(validdesmat))   
     }
    SE.Krig = sqrt(mean(ressave^2))
    SE.Krig
    

data(meuse.grid)

Zpredict.Krig<- predict( fitmeuse, meuse.grid[,1:2],drop.Z=T)  # smoothing surface, trend + nonparametric 


#plot(meuse.grid[,1:2],type="n")
#for (k in 1:length(meuse.grid[,1]))
#{points(meuse.grid[k,1],meuse.grid[k,2],pch=19,col=gray((Zpredict.Krig[k]-min(Zpredict.Krig))/(max(Zpredict.Krig)-min(Zpredict.Krig))))
#}




xvec= seq(min(meuse.grid[,1]),max(meuse.grid[,1]),by=40)
yvec=seq(min(meuse.grid[,2]),max(meuse.grid[,2]),by=40)
XYmat.Krig.trend.nonparam = matrix(NA,nrow=length(xvec),ncol=length(yvec))
XYmat.Krig.trend          = XYmat.Krig.trend.nonparam 
XYmat.Krig.nonparam       = XYmat.Krig.trend.nonparam 
XYmatbis=matrix(NA,nrow=length(xvec),ncol=length(yvec))
for (indexX in 1:length(xvec))
{
yfill=meuse.grid[which(meuse.grid[,1]==xvec[indexX]),2]
for (indexyfill in 1:length(yfill))
{indexY=which(yvec==yfill[indexyfill])
XYmat.Krig.trend.nonparam[indexX,indexY] = predict( fitmeuse, t(c(xvec[indexX],yvec[indexY])),drop.Z=T)                 # smoothing surface, trend + nonparametric
XYmat.Krig.trend[indexX,indexY]          = predict( fitmeuse, t(c(xvec[indexX],yvec[indexY])),drop.Z=T, just.fixed=T)   # trend 
XYmat.Krig.nonparam[indexX,indexY]       = XYmat.Krig.trend.nonparam[indexX,indexY]-XYmat.Krig.trend[indexX,indexY]     # nonparametric 
XYmatbis[indexX,indexY]=sum(fitmeuse$d[1:6]*c(1,xvec[indexX],yvec[indexY],xvec[indexX]^2,xvec[indexX]*yvec[indexY],yvec[indexY]^2))
}
}

windows()
image(xvec,yvec,XYmat.Krig.trend.nonparam)
contour(xvec,yvec,XYmat.Krig.trend.nonparam,add=T)

persp3d(xvec,yvec,XYmat.Krig.trend.nonparam, col='red', alpha=1)

windows()
image(xvec,yvec,XYmat.Krig.trend)
contour(xvec,yvec,XYmat.Krig.trend,add=T)

persp3d(xvec,yvec,XYmat.Krig.trend, col='red', alpha=1)


windows()
image(xvec,yvec,XYmat.Krig.nonparam)
contour(xvec,yvec,XYmat.Krig.nonparam,add=T)

persp3d(xvec,yvec,XYmat.Krig.nonparam, col='red', alpha=1)





#####################
#####################
#####################
#####################
#####################
#####################

#
#plot(Zpredict.Krig)
#plot(Zpredict.Krig[2101:2800])
#abline(h=8.42)
#which(Zpredict.Krig[2101:2800]<8.42)
#points(which(Zpredict.Krig[2101:2800]<8.42),Zpredict.Krig[2101:2800][which(Zpredict.Krig[2101:2800]<8.42)],col=2)
#
#
#plot(Zpredict.Krig)
#points((2100+which(Zpredict.Krig[2101:2800]<8.42)),Zpredict.Krig[(2100+which(Zpredict.Krig[2101:2800]<8.42))],col=2)
#
#meuse.grid[(2100+which(Zpredict.Krig[2101:2800]<8.42)),1:2]
#predict( fitmeuse, meuse.grid[(2100+which(Zpredict.Krig[2101:2800]<8.42)),1:2],drop.Z=T)
#abline(h=predict( fitmeuse, meuse.grid[(2100+which(Zpredict.Krig[2101:2800]<8.42)),1:2],drop.Z=T))
#
#
#aggiunta=meuse.grid[(2100+which(Zpredict.Krig[2101:2800]<8.42)),1:2]
#rownames(aggiunta)=c("pa","pb","pc")
#meuse.grid.aug=rbind(meuse.grid[,1:2],aggiunta)
#
#dist.meuse.grid.aug=as.matrix(dist(meuse.grid.aug))
#
#which(dist.meuse.grid.aug[,length(meuse.grid[,1])+1]==0)
#
#pgrida=meuse.grid[which(dist.meuse.grid.aug[,length(meuse.grid[,1])+1]==0)[1],1:2]
#
#
#
#which(dist.meuse.grid.aug[,length(meuse.grid[,1])+2]==0)
#
#pgridb=meuse.grid[which(dist.meuse.grid.aug[,length(meuse.grid[,1])+2]==0)[1],1:2]
#
#
#
#which(dist.meuse.grid.aug[,length(meuse.grid[,1])+3]==0)
#
#pgridc=meuse.grid[which(dist.meuse.grid.aug[,length(meuse.grid[,1])+3]==0)[1],1:2]
#
#
#names(p)=c("x","y")
#
#Ngrid=length(meuse.grid[,1])
#
#meuse.grid.aug.p=rbind(meuse.grid[,1:2],p)
#dist.meuse.grid.aug.p=as.matrix(dist(meuse.grid.aug.p))
#
#
#dist.meuse.grid.aug.p.diag=dist.meuse.grid.aug.p+diag(100000000,nrow=length(dist.meuse.grid.aug.p[,1]))
#
#minpointdist=numeric(np)
#for (prova in 1:np)
#{
#minpointdist[prova]=min(dist.meuse.grid.aug.p.diag[,Ngrid+prova])
#}
#
#
#
#
#p[95,] 
#pgridc
#
#p[102,] 
#pgrida
#
#p[101,] 
#pgridb
#
#
#
#XYmatmod=XYmat.Krig
#XYmatmod[which(xvec==pgrida[,1]),which(yvec==pgrida[,2])]=NA
#XYmatmod[which(xvec==pgridb[,1]),which(yvec==pgridb[,2])]=NA
#XYmatmod[which(xvec==pgridc[,1]),which(yvec==pgridc[,2])]=NA
#
#
#image(xvec,yvec,XYmat.Krig)
#contour(xvec,yvec,XYmat.Krig,add=T)
#
#
#
#persp3d(xvec,yvec,XYmatmod, col='red', alpha=1)

#####################
#####################
#####################
#####################
#####################
#####################

### Thin Plate Splines

lambdafix=0.0001

fitmeuseTps=Tps(x=p,Y=data[,2],Z=desmat,lambda=lambdafix)
surface( fitmeuseTps)




erroriTps=data[,2]-fitmeuseTps$fitted.values 
erroriTps-fitmeuseTps$residuals

predict(fitmeuseTps)-fitmeuseTps$fitted.values

predict(fitmeuseTps,  drop.Z=T)+fitmeuseTps$d[(length(fitmeuseTps$d)-1)]*desmat[,1]+ fitmeuseTps$d[length(fitmeuseTps$d)]*desmat[,2] -fitmeuseTps$fitted.values
 



#  compute the squared distances and the residual products

dvec.Tps = numeric(np*(np-1)/2)
rvec.Tps = dvec.Tps

mm = 0
for (i in 2:np)
 {   
   for (j in 1:(i-1))
     {  mm = mm + 1
        dvec.Tps[mm] = sum((p[i,]-p[j,])^2)
        rvec.Tps[mm] = fitmeuseTps$residuals[i] * fitmeuseTps$residuals[j]
     }
 }


windows()
plot(log10(dvec.Tps), rvec.Tps,xlab="log_{10} Squared distance")
abline(h=0,col=2,lwd=2,lty="dashed")



np=155
ressave=numeric(np)
    for (idrop in 1:np)
     {  
        print(idrop)
        modelp=p[-idrop,]
        modeldata = data[-idrop,]
        modeldesmat = desmat[-idrop,]
        modelfitmeuse <- Tps(x=modelp,Y=modeldata[,2],Z=modeldesmat,lambda=lambdafix)
        validp=p[idrop,]
        validdesmat=desmat[idrop,]
        ressave[idrop] <- data[idrop,2] - predict(modelfitmeuse, validp, Z=t(validdesmat))  
     }
    SE.Tps = sqrt(mean(ressave^2))
    SE.Tps
    

Zpredict.Tps<- predict( fitmeuseTps, meuse.grid[,1:2],drop.Z=T)  # smoothing surface, trend + nonparametric





xvec= seq(min(meuse.grid[,1]),max(meuse.grid[,1]),by=40)
yvec=seq(min(meuse.grid[,2]),max(meuse.grid[,2]),by=40)
XYmat.Tps.trend.nonparam = matrix(NA,nrow=length(xvec),ncol=length(yvec))
XYmat.Tps.trend          = XYmat.Tps.trend.nonparam 
XYmat.Tps.nonparam       = XYmat.Tps.trend.nonparam 
for (indexX in 1:length(xvec))
{
yfill=meuse.grid[which(meuse.grid[,1]==xvec[indexX]),2]
for (indexyfill in 1:length(yfill))
{indexY=which(yvec==yfill[indexyfill])
XYmat.Tps.trend.nonparam[indexX,indexY] = predict( fitmeuseTps, t(c(xvec[indexX],yvec[indexY])),drop.Z=T)                 # smoothing surface, trend + nonparametric
XYmat.Tps.trend[indexX,indexY]          = predict( fitmeuseTps, t(c(xvec[indexX],yvec[indexY])),drop.Z=T, just.fixed=T)   # trend 
XYmat.Tps.nonparam[indexX,indexY]       = XYmat.Tps.trend.nonparam[indexX,indexY]-XYmat.Tps.trend[indexX,indexY]         # nonparametric 
}
}





windows()
image(xvec,yvec,XYmat.Tps.trend.nonparam)
contour(xvec,yvec,XYmat.Tps.trend.nonparam,add=T)

persp3d(xvec,yvec,XYmat.Tps.trend.nonparam, col='red', alpha=1)

windows()
image(xvec,yvec,XYmat.Tps.trend)
contour(xvec,yvec,XYmat.Tps.trend,add=T)

persp3d(xvec,yvec,XYmat.Tps.trend, col='red', alpha=1)


windows()
image(xvec,yvec,XYmat.Tps.nonparam)
contour(xvec,yvec,XYmat.Tps.nonparam,add=T)

persp3d(xvec,yvec,XYmat.Tps.nonparam, col='red', alpha=1)











##
indtri   = matrix(1:nt,ncol=1)

xmin = min(p_2col[,1])
xmax = max(p_2col[,1])
nx   = 101
X    = matrix(seq(xmin, xmax, len=nx),ncol=1)
 
ymin = min(p_2col[,2])
ymax = max(p_2col[,2])
ny   = 101
Y    = matrix(seq(ymin, ymax, len=ny),ncol=1)    
  
Xmat = X %*% matrix(1,nrow=1,ncol=ny)
Ymat = matrix(1,nrow=nx,ncol=1) %*% t(Y)
Xvec = NULL
for (numc in 1:nx)
 {Xvec=c(Xvec,Xmat[,numc])}
Yvec = NULL
for (numc in 1:ny)
 {Yvec=c(Yvec,Ymat[,numc])}
###

evalmat = eval.FEM.fd(Xvec, Yvec, ZincMeusefdCovar$felsplobj)
evalmati = matrix(evalmat[,1],nrow=nx, ncol=ny, byrow=F)
 
par(mfrow=c(2,2))
zlimits=c(min(c(XYmat.Tps.trend.nonparam[which.min(XYmat.Tps.trend.nonparam)],evalmati[which.min(evalmati)])),max(c(XYmat.Tps.trend.nonparam[which.max(XYmat.Tps.trend.nonparam)],evalmati[which.max(evalmati)])))

image(xvec,yvec,XYmat.Tps.trend.nonparam,main=paste(c("Thin plate splines, betahat= " ,fitmeuseTps$d[4:5])),zlim=zlimits,xlab="",ylab="")
contour(xvec,yvec,XYmat.Tps.trend.nonparam,add=T)

image(X,Y,evalmati,main=paste(c("Finite Elements splines,  betahat= " , betahat[,1])),ylim=c(min(meuse.grid[,2]),max(meuse.grid[,2])),xlim=c(min(meuse.grid[,1]),max(meuse.grid[,1])),zlim=zlimits,xlab="",ylab="")
contour(X,Y,evalmati,add=T)   

image(xvec,yvec,XYmat.Krig.trend.nonparam,main=paste(c("Filtered kriging,  betahat= " , fitmeuse$d[7:8])),zlim=zlimits,xlab="",ylab="")
contour(xvec,yvec,XYmat.Krig.trend.nonparam,add=T)







par(mfrow=c(2,2))
ylimits=c(min(c(rvec.Tps,rvec.FES,rvec.Krig)),max(c(rvec.Tps,rvec.FES,rvec.Krig)))


plot(log10(dvec.Tps), rvec.Tps,xlab="log_{10} Squared distance",main="Thin plate splines",ylab="",ylim=ylimits)
abline(h=0,col=2,lwd=2,lty="dashed")



plot(log10(dvec.FES), rvec.FES,xlab="log_{10} Squared distance",main="Finite Elements splines",ylab="",ylim=ylimits)
abline(h=0,col=2,lwd=2,lty="dashed")

   

plot(log10(dvec.Krig), rvec.Krig,xlab="log_{10} Squared distance",main="Filtered kriging",ylab="",ylim=ylimits)
abline(h=0,col=2,lwd=2,lty="dashed")
















data(meuse.grid)

minmat=NULL
maxmat=NULL
plot(meuse.grid$x,meuse.grid$y)
xvec= seq(min(meuse.grid[,1]),max(meuse.grid[,1]),by=40)
yvec=seq(min(meuse.grid[,2]),max(meuse.grid[,2]),by=40)
for (indexX in 1:length(xvec))
{
yfill=meuse.grid[which(meuse.grid[,1]==xvec[indexX]),2]
minmat=rbind(minmat,c(xvec[indexX],min(yfill)))
maxmat=rbind(maxmat,c(xvec[indexX],max(yfill)))
points(xvec[indexX],min(yfill),col=2)
points(xvec[indexX],max(yfill),col=2)
}


plot(meuse.grid$x,meuse.grid$y)
points(rbind(minmat,maxmat),col=2)
abline(h=minmat[58,2])
abline(h=minmat[57,2])
abline(v=minmat[58,1])
abline(v=180300)


diffmat=NULL

for (indeyY in (which(yvec==minmat[57,2])+1):(which(yvec==minmat[58,2])-1))
{
xfill=meuse.grid[which(meuse.grid[,2]==yvec[indeyY]),1]
diffmat=rbind(diffmat,c(max(xfill),yvec[indeyY]))
points(max(xfill),yvec[indeyY],col=3)
}


fullmat=rbind(minmat[1:57,],diffmat,minmat[58:length(minmat[,1]),],maxmat[length(maxmat[,1]):1,])


plot(meuse.grid$x,meuse.grid$y)
points(fullmat,col=2)

plot(fullmat[,1])


plot(fullmat[,2])


ng=length(fullmat[,1])

ac=numeric(ng)
for (index in 2:ng)
{ac[index]=ac[(index-1)]+sqrt((fullmat[index,1]-fullmat[(index-1),1])^2+(fullmat[index,2]-fullmat[(index-1),2])^2)
}


plot(ac)



plot(ac,fullmat[,1])
plot(ac,fullmat[,2])

fullmat=cbind(fullmat,numeric(ng))

nkiniziale=round(ng/10)
passo=round(ng/(nkiniziale+1))
kstart=numeric(nkiniziale)
for (nodiniz in 1:nkiniziale)
{kstart[nodiniz]=nodiniz*passo
}



source("p3Dfks_functions.R")

C=500
m=2

risultati=SARSmain3D(t(fullmat),ac,kstart,m,1,C)


knots.data=risultati[[1]]

nk=length(knots.data)


splinehatM=matrix(0,nrow=2,ncol=ng)
beta=matrix(0,nrow=2,ncol=(nk+m))
beta1=matrix(0,nrow=2,ncol=(nk+m-1))
beta2=matrix(0,nrow=2,ncol=(nk+m-2))
beta3=matrix(0,nrow=2,ncol=(nk+m-3))

for (k in 1:2)
 {
 splinehat=SARSregressionspline(fullmat[,k],ac,knots.data,m)
 splinehatM[k,]=splinehat[[2]]
 beta[k,]=splinehat[[1]]
 }


par(mfrow=c(1,2))
plot(ac,fullmat[,1])
points(ac,splinehatM[1,],col=2,type="l")
abline(v=ac[knots.data])

plot(ac,fullmat[,2])
points(ac,splinehatM[2,],col=2,type="l")
abline(v=ac[knots.data])



par(mfrow=c(1,1))
plot(fullmat,col="gray")
points(splinehatM[1,],splinehatM[2,],col=2,type="l")
points(splinehatM[1,knots.data],splinehatM[2,knots.data])


knots.data.mod=knots.data
knots.data.mod[17]=107
knots.data.mod=c(1,knots.data.mod)


knots.data.mod=c(knots.data.mod,119)
knots.data.mod=sort(knots.data.mod)


knots.data.mod=c(knots.data.mod,c(126))
knots.data.mod=sort(knots.data.mod)

knots.data.mod[23]=129


knots.data.mod=c(knots.data.mod,c(131))
knots.data.mod=sort(knots.data.mod)


knots.data.mod=c(knots.data.mod,c(134,136,138,139))
knots.data.mod=sort(knots.data.mod)



knots.data.mod=c(knots.data.mod,c(143))
knots.data.mod=sort(knots.data.mod)


knots.data.mod=c(knots.data.mod,c(146))
knots.data.mod=sort(knots.data.mod)


knots.data.mod[33]=149


knots.data.mod=c(knots.data.mod,c(152))
knots.data.mod=sort(knots.data.mod)


knots.data.mod=c(knots.data.mod,c(155,157))
knots.data.mod=sort(knots.data.mod)


knots.data.mod=c(knots.data.mod,c(159))
knots.data.mod=sort(knots.data.mod)


knots.data.mod=c(knots.data.mod,c(162))
knots.data.mod=sort(knots.data.mod)


knots.data.mod=c(knots.data.mod,c(165,166,167,168))
knots.data.mod=sort(knots.data.mod)


knots.data.mod=c(knots.data.mod,c(171,173))
knots.data.mod=sort(knots.data.mod)

knots.data.mod=c(knots.data.mod,c(176,178))
knots.data.mod=sort(knots.data.mod)


knots.data.mod=c(knots.data.mod,c(180))
knots.data.mod=sort(knots.data.mod)


knots.data.mod=c(knots.data.mod,c(185,188))
knots.data.mod=sort(knots.data.mod)


knots.data.mod=c(knots.data.mod,c(7))
knots.data.mod=sort(knots.data.mod)

knots.data.mod=c(knots.data.mod,c(42))
knots.data.mod=sort(knots.data.mod)

knots.data.mod=c(knots.data.mod,c(32))
knots.data.mod=sort(knots.data.mod)


par(mfrow=c(1,1))
plot(fullmat,col="gray")
points(splinehatM[1,],splinehatM[2,],col=2,type="l")
points(splinehatM[1,knots.data.mod],splinehatM[2,knots.data.mod])
points(p[,1],p[,2],col="blue")

points(splinehatM[1,44],splinehatM[2,44],col="green")




nk=length(knots.data.mod)


splinehatM=matrix(0,nrow=2,ncol=ng)
beta=matrix(0,nrow=2,ncol=(nk+m))
beta1=matrix(0,nrow=2,ncol=(nk+m-1))
beta2=matrix(0,nrow=2,ncol=(nk+m-2))
beta3=matrix(0,nrow=2,ncol=(nk+m-3))

for (k in 1:2)
 {
 splinehat=SARSregressionspline(fullmat[,k],ac,knots.data.mod,m)
 splinehatM[k,]=splinehat[[2]]
 beta[k,]=splinehat[[1]]
 }




par(mfrow=c(1,2))
plot(ac,fullmat[,1])
points(ac,splinehatM[1,],col=2,type="l")
abline(v=ac[knots.data.mod])

plot(ac,fullmat[,2])
points(ac,splinehatM[2,],col=2,type="l")
abline(v=ac[knots.data.mod])



par(mfrow=c(1,1))
plot(fullmat,col="gray")
points(splinehatM[1,],splinehatM[2,],col=2,type="l")
points(splinehatM[1,knots.data.mod],splinehatM[2,knots.data.mod])




fullmat[knots.data.mod,1:2]
write.table(fullmat[knots.data.mod,1:2],"polygon_meuse.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)


write.table(knots.data,"knots.data_meuse.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(knots.data,"knots.data.mod_meuse.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)



polygon_meuse=read.table('polygon_meuse.txt', header=FALSE)
knots.data_meuse=read.table('knots.data_meuse.txt', header=FALSE)
knots.data.mod_meuse=read.table('knots.data.mod_meuse.txt', header=FALSE)
