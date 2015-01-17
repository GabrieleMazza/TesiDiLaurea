#  Smoothing of the Meuse data in Bivand et al. 
#  
#  Last modified on 26 March 2014 by Laura Sangalli
# !!! Attenzione: c'è una versione molto più estesa di questo esempio, 
#     con anche analisi dei residui, varie mesh, kriging e thin-plate splines



setwd("C:/Users/Laura/Desktop/SR-PDE_short-tutorials_R")  
source("2013_SSR_AllFunctions.R")
library(fda)
library(rgl)

#  The Meuse data are displayed in R, where they are included 
#  with the geoR package, loaded by command data(meuse)

library(geoR)
data(meuse)
help(meuse)

#  The variable labels in the R file are:

#  row  x  y cadmium copper lead zinc   elev   dist   om ffreq soil

#  Here we shall directly read the .txt also used in the analyses carried out in Matlab 

MeuseData=read.table('MeuseData.txt', header=FALSE)

#  define the vertices as points where there are data

p = MeuseData[,2:3]
np = dim(p)[[1]]


#  plot these points

plot(p[,1],p[,2])

# upload the mesh (already computed with Matlab)

MeuseTri = as.matrix(read.table('MeuseTri.txt', header=F))

#  The triangle index matrix

t = MeuseTri
nt = dim(t)[[1]]

# plot the mesh

plot(p[,1],p[,2],type="n")
polygon(c(p[MeuseTri[1,1],1],p[MeuseTri[1,2],1],p[MeuseTri[1,3],1]),c(p[MeuseTri[1,1],2],p[MeuseTri[1,2],2],p[MeuseTri[1,3],2]))
for (ne in 1:nt)
{polygon(c(p[MeuseTri[ne,1],1],p[MeuseTri[ne,2],1],p[MeuseTri[ne,3],1]),c(p[MeuseTri[ne,1],2],p[MeuseTri[ne,2],2],p[MeuseTri[ne,3],2]))
}

print(paste("Number of edited triangles = ",nt))

#  we don't need an edge matrix for these analyses

e = NULL

#  set up the FEM basis object and plot it

order=1

basisobj = create.FEM.basis(p, e, t, order)


#  set up a dummy FEM functional data object

Meusefd = fd(numeric(basisobj$nbasis),basisobj)

#  set up the data array using zinc concentration

data = matrix(0,nrow=np,ncol=2)
data[,1] = 1:np
data[,2] = log(MeuseData[,7])



#  use as covariate the distance from river; for the moment we use as regressors
#  sqrt(distance) and 
#  elev

desmat = matrix(1,nrow=np,ncol=2)
desmat[,1] = sqrt(MeuseData[,9])
desmat[,2] = (MeuseData[,8])


#  smooth the data using covariates

lambdaCovar = 10^(3.5)
ZincMeusefdCovar = smooth.FEM.fd.Covar(data,desmat,Meusefd,lambdaCovar)



fhat=ZincMeusefdCovar$felsplobj$coef
fhat

betahat= ( solve( t(desmat) %*% desmat ) ) %*% t(desmat) %*% (data[,2]-fhat[1:np,])
betahat

#  plot smoothing surface computed using covariates

plot.FEM(ZincMeusefdCovar$felsplobj)




########################################
########################################
########################################
########################################
########################################
########################################


########################################
#  !!!!!!! HERE we use the full Meuse bank area
########################################



MeuseTri = as.matrix(read.table('MeuseTri_constr.txt', header=F))
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

order=2


#  set up the FEM basis object and plot it

basisobj = create.FEM.basis(p_2col, e, t, order)



#  set up a dummy FEM functional data object

Meusefd = fd(numeric(basisobj$nbasis),basisobj)




#  smooth the data using covariates

lambdaCovar = 10^(3.5)
ZincMeusefdCovar = smooth.FEM.fd.Covar(data,desmat,Meusefd,lambdaCovar)



fhat=ZincMeusefdCovar$felsplobj$coef
fhat
betahat= ( solve( t(desmat) %*% desmat ) ) %*% t(desmat) %*% (data[,2]-fhat[1:np,])
betahat

#  plot smoothing surface computed using covariates
windows()
plot.FEM(ZincMeusefdCovar$felsplobj)
