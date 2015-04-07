#######################################################
#  Data on Isle of Montreal
#  Upload GOOGLE maps and plot on them!
#######################################################

#  Last modified 04 December 2011 by Laura Sangalli



setwd("E:/Dati/canada/fdaR") 
source("2013_SSR_AllFunctions.R")


library(RgoogleMaps)
library(rgl)
library(fda)


# Upload Montreal_tractcenter and go back to real longitude and latitude:


Montreal_tractcenter     = 
    as.matrix(read.table('Montreal_tractcenter.txt',     header=FALSE))
Montreal_Tri_abconstrtot = 
    as.matrix(read.table('Montreal_Tri_abconstrtot.txt', header=FALSE))
Montreal_bound_diric = 
    as.matrix(read.table('laura_Montreal_Boundary.txt', header=FALSE))


# in Montreal_Tri_abconstrtot:
# the first 100 edges are for external boundary
# from 101 to 118 are for Dorval
# from 119 to 133 are for Eastend


tractcenter = Montreal_tractcenter + cbind(rep(-73.6939,493),rep(45.5100,493))

#  do the same for boundaries

p = as.matrix(read.table('Montreal_p.txt', header=F))
p = p + cbind(rep(-73.6939,(dim(p)[[1]])),rep(45.5100,(dim(p)[[1]])))


# Upload Triangulation

t = as.matrix(read.table('Montreal_Tri.txt', header=F))

nt = dim(t)[[1]]


# Upload Isle of Montreal income and density data

Montreal_income=as.matrix(read.table('Montreal_income.txt', header=FALSE))
Montreal_popden=as.matrix(read.table('Montreal_popden.txt', header=FALSE))
Montreal_incden=as.matrix(read.table('Montreal_incden.txt', header=FALSE))

np = dim(Montreal_income)[[1]]




#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################
#####################################################################################################





# Create a png image that we can use as background image where to overly specific plots and images

# Extrema of the map

centrale=c( 45.54755,-73.72085)

mappa <- cbind.data.frame(lat = c( 45.4071, centrale[1], 45.6931),
                          lon = c(-73.9651, centrale[2],-73.4766))

######mappa <- cbind.data.frame(lat = c( centrale[1]-0.15, centrale[1],  centrale[1]+0.15),
######                          lon = c( centrale[2]-0.15, centrale[2],  centrale[2]+0.15))

# save google images

bb <- qbbox(lat = mappa[,'lat'], lon = mappa[,'lon'], 
            margin=list(m=c(0,0,0,0), TYPE='perc'))

GetMap(center = centrale, maptype = 'satellite', 
       format = 'png32', destfile = 'Montreal_google-satellite.png', zoom=10)

Map <- GetMap.bbox(bb$lonR, bb$latR, centrale, destfile = 'Montreal_google-map.png',  maptype = 'satellite',NEWMAP = FALSE, GRAYSCALE  = T)



### Now use this to create our figure

windows()
PlotOnStaticMap(Map, add = FALSE, TrueProj=F,  FUN = points)    # plot the background google image

# and add whatever we like:

# add external boundary

for (nb in 1:100)
{
    points(p[Montreal_Tri_abconstrtot[nb,1:2],1],p[Montreal_Tri_abconstrtot[nb,1:2],2],type="l",lwd=3,col="yellow")
}


# add two holes corresponding to Dorval and Eastend, with the corresponding internal boundaries

#polygon(p[Montreal_Tri_abconstrtot[101:118,1],1],p[Montreal_Tri_abconstrtot[101:118,1],2],col="grey")

#polygon(p[Montreal_Tri_abconstrtot[119:133,1],1],p[Montreal_Tri_abconstrtot[119:133,1],2],col="grey")

for (nb in 101:118)
{
    points(p[Montreal_Tri_abconstrtot[nb,1:2],1],p[Montreal_Tri_abconstrtot[nb,1:2],2],type="l",lwd=2,col="yellow")
}

for (nb in 119:133)
{
    points(p[Montreal_Tri_abconstrtot[nb,1:2],1],p[Montreal_Tri_abconstrtot[nb,1:2],2],type="l",lwd=2,col="yellow")
}


for (nb in which(Montreal_bound_diric[1:132,1]==1 & Montreal_bound_diric[2:133,1]==1))
{
    points(p[Montreal_Tri_abconstrtot[nb,1:2],1],
           p[Montreal_Tri_abconstrtot[nb,1:2],2],type="l",lwd=3,col="blue")
}



########points( p[(493+which(Montreal_bound_diric[1:133,1]==1)),1], p[(493+which(Montreal_bound_diric[1:133,1]==1)),2],col=2,cex=2)


# add tractcenters

points(tractcenter[,1],tractcenter[,2],col="yellow",pch=19,cex=0.5)



windows()
PlotOnStaticMap(Map, add = FALSE, TrueProj=F,  FUN = points)    # plot the background google image


# Now let us plot insted the triangulation

#polygon(p[Montreal_Tri_abconstrtot[101:118,1],1],p[Montreal_Tri_abconstrtot[101:118,1],2],col="grey")
#polygon(p[Montreal_Tri_abconstrtot[119:133,1],1],p[Montreal_Tri_abconstrtot[119:133,1],2],col="grey")
points(p[,1],p[,2],type="n")
polygon(c(p[t[1,1],1],p[t[1,2],1],p[t[1,3],1]),c(p[t[1,1],2],p[t[1,2],2],p[t[1,3],2]),border="yellow")
for (ne in 1:nt)
{polygon(c(p[t[ne,1],1],p[t[ne,2],1],p[t[ne,3],1]),c(p[t[ne,1],2],p[t[ne,2],2],p[t[ne,3],2]),border="yellow")
}




for (nb in which(Montreal_bound_diric[1:132,1]==1 & Montreal_bound_diric[2:133,1]==1))
{
    points(p[Montreal_Tri_abconstrtot[nb,1:2],1],
           p[Montreal_Tri_abconstrtot[nb,1:2],2],type="l",lwd=3,col="red")
}





######################################################################
# To plot one data point at a time

indtract=152



#########################################
######## START SELECTING HERE ###########


graphics.off()
indtract= indtract+1

# ZOOM
centralezoom= c(tractcenter[indtract,2],tractcenter[indtract,1])

# create zoomed-in map
mappazoom <- cbind.data.frame(lat = c(centralezoom[1]-0.02, centralezoom[1], centralezoom[1]+0.02),
                              lon = c(centralezoom[2]-0.02, centralezoom[2], centralezoom[2]+0.02))

bbzoom <- qbbox(lat = mappazoom[,'lat'], lon = mappazoom[,'lon'], 
                margin=list(m=c(0,0,0,0), TYPE='perc'))

GetMap(center = centralezoom, maptype = 'satellite', 
       format = 'png32', destfile = 'Montreal_google-satellite-zoom.png', zoom=10)

Mapzoom <- GetMap.bbox(bbzoom$lonR, bbzoom$latR, centralezoom, destfile = 'Montreal_google-mapzoom.png', NEWMAP = FALSE)

# plot zoomed-in map around specific data point

PlotOnStaticMap(Mapzoom, add = FALSE, TrueProj=F,  FUN = points)
points(tractcenter[indtract,1],tractcenter[indtract,2],col="red",pch=19,cex=2)
print(paste(indtract,", lat= ", centralezoom[1],", lon= ", centralezoom[2]))  # print data index, and latitude and longitude of the data point


######## END SELECTING HERE ###########
#######################################







###################################################################################################################
######### BOUNDARY



indbd=np

#########################################
######## START SELECTING HERE ###########


graphics.off()
indbd= indbd+1

# ZOOM
centralezoom= c(p[indbd,2],p[indbd,1])

# create zoomed-in map
mappazoom <- cbind.data.frame(lat = c(centralezoom[1]-0.03, centralezoom[1], centralezoom[1]+0.03),
                              lon = c(centralezoom[2]-0.03, centralezoom[2], centralezoom[2]+0.03))

bbzoom <- qbbox(lat = mappazoom[,'lat'], lon = mappazoom[,'lon'], 
                margin=list(m=c(0,0,0,0), TYPE='perc'))

GetMap(center = centralezoom, maptype = 'satellite', 
       format = 'png32', destfile = 'Montreal_google-satellite-zoom.png', zoom=10)

Mapzoom <- GetMap.bbox(bbzoom$lonR, bbzoom$latR, centralezoom, destfile = 'Montreal_google-mapzoom.png', NEWMAP = FALSE)

# plot zoomed-in map around specific data point

PlotOnStaticMap(Mapzoom, add = FALSE, TrueProj=F,  FUN = points)
points(p[indbd,1],p[indbd,2],col="red",pch=19,cex=2)
print(paste(indbd,", lat= ", centralezoom[1],", lon= ", centralezoom[2]))  # print data index, and latitude and longitude of the data point


######## END SELECTING HERE ###########
#######################################
















#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################




# ANALYSES


Montreal_indicator_resNOres = as.matrix(read.table('Montreal_indicator_resNOres.txt', header=FALSE))

data = matrix(0,nrow=np,ncol=2)
data[,1] = 1:np 
##########data[,2] = sqrt(Montreal_popden/10000)

data[,2] = Montreal_popden/1000



centrale=c( 45.54755,-73.72085)

mappazoom <- cbind.data.frame(lat = c( centrale[1]-0.10, centrale[1],  centrale[1]+0.10),
                              lon = c( centrale[2]-0.10, centrale[2],  centrale[2]+0.10))

bbzoom <- qbbox(lat = mappazoom[,'lat'], lon = mappazoom[,'lon'], 
                margin=list(m=c(0,0,0,0), TYPE='perc'))

GetMap(center = centrale, maptype = 'satellite', 
       format = 'png32', destfile = 'Montreal_google-satellite-zoom.png', zoom=10)

Mapzoom <- GetMap.bbox(bbzoom$lonR, bbzoom$latR, centrale, destfile = 'Montreal_google-map-zoom.png',NEWMAP = FALSE, GRAYSCALE  = T)

PlotOnStaticMap(Mapzoom, add = FALSE, TrueProj=F,  FUN = points)    # plot the background google image

points(tractcenter[which(Montreal_indicator_resNOres==1),1],tractcenter[which(Montreal_indicator_resNOres==1),2],col="yellow",pch=19,cex=1)
points(tractcenter[which(Montreal_indicator_resNOres==0),1],tractcenter[which(Montreal_indicator_resNOres==0),2],col="red",pch=19,cex=1)
points(tractcenter[c(21,203),1],tractcenter[c(21,203),2],col="blue",pch=19,cex=1)


windows()
PlotOnStaticMap(Mapzoom, add = FALSE, TrueProj=F,  FUN = points)    # plot the background google image
for (k in 1:np)
{points(tractcenter[k,1],p[k,2],pch=19,cex=2*(data[k,2]-min(data[,2]))/(max(data[,2])-min(data[,2])))
}



Montreal_ind = Montreal_indicator_resNOres


windows()
boxplot(Montreal_popden[which(Montreal_ind==1)],
        Montreal_popden[which(Montreal_ind==0)], names=c("res","NO res"))
wilcox.test(Montreal_popden[which(Montreal_ind==1)],
            Montreal_popden[which(Montreal_ind==0)])


max(Montreal_popden)
which.max(Montreal_popden)
max(Montreal_popden[-which.max(Montreal_popden)])
which.max(Montreal_popden[-which.max(Montreal_popden)])

Montreal_ind[which(Montreal_popden >= max(Montreal_popden[-which.max(Montreal_popden)]))]

Montreal_ind[which(Montreal_popden >= max(Montreal_popden[-which.max(Montreal_popden)]))] = 1


# How census tract with 0 population density have been classified


PlotOnStaticMap(Mapzoom, add = FALSE, TrueProj=F,  FUN = points)    # plot the background google image

points(tractcenter[which(Montreal_popden==0),1],tractcenter[which(Montreal_popden==0),2],col="blue",pch=19,cex=1)

Montreal_ind[which(Montreal_popden==0)]

# Since a census tract with 0 population density cannot be a residential tract, I'll set these to NON residential

Montreal_ind[which(Montreal_popden==0)] = 0



boxplot(Montreal_popden[which(Montreal_ind==1)],
        Montreal_popden[which(Montreal_ind==0)], names=c("residential","NON residential"))
wilcox.test(Montreal_popden[which(Montreal_ind==1)],
            Montreal_popden[which(Montreal_ind==0)])



which(Montreal_popden[which(Montreal_ind==1)]<100)

which(Montreal_popden<100)

Montreal_ind[which(Montreal_popden<100)] = 0

graphics.off()
boxplot(Montreal_popden[which(Montreal_ind==1)],
        Montreal_popden[which(Montreal_ind==0)], names=c("residential","NON residential"))
wilcox.test(Montreal_popden[which(Montreal_ind==1)],
            Montreal_popden[which(Montreal_ind==0)])


breaks = (0:24)*1000
par(mfrow=c(2,2))
hist(Montreal_popden[which(Montreal_ind==1)],xlim=c(min(Montreal_popden),max(Montreal_popden)),main="residential",breaks=breaks,xlab="popden")
hist(Montreal_popden[which(Montreal_ind==0)],xlim=c(min(Montreal_popden),max(Montreal_popden)),main="NON residential",xlab="popden")


##graphics.off()
##PlotOnStaticMap(Mapzoom, add = FALSE, TrueProj=F,  FUN = points)    # plot the background google image
##
##points(tractcenter[which(Montreal_ind==1),1],tractcenter[which(Montreal_ind==1),2],col="yellow",pch=19,cex=1)
##points(tractcenter[which(Montreal_ind==0),1],tractcenter[which(Montreal_ind==0),2],col="red",pch=19,cex=1)
##
##
##windows()
##PlotOnStaticMap(Mapzoom, add = FALSE, TrueProj=F,  FUN = points)    # plot the background google image
##title("change")
##points(tractcenter[which((Montreal_indicator_resNOres-Montreal_ind)==1),1],tractcenter[which((Montreal_indicator_resNOres-Montreal_ind)==1),2],col="blue",pch=19,cex=1)



desmat1 = matrix(1,nrow=np,ncol=1)
desmat1[,1] =  Montreal_ind


desmat2 = matrix(1,nrow=np,ncol=2)
desmat2[,1] =  Montreal_ind
desmat2[,2] = sqrt(Montreal_income/10000)


par(mfrow=c(2,2))
plot(desmat2[which(Montreal_ind==1),2],data[which(Montreal_ind==1),2],xlim=c(min(desmat2[,2]),max(desmat2[,2])),ylim=c(min(data[,2]),max(data[,2])),main="residential",xlab="sqrt(income/10000)",ylab="sqrt(popden/1000)")
hist(Montreal_popden[which(Montreal_ind==1)],xlim=c(min(Montreal_popden),max(Montreal_popden)),main="residential",xlab="popden",breaks=breaks)
plot(desmat2[which(Montreal_ind==0),2],data[which(Montreal_ind==0),2],xlim=c(min(desmat2[,2]),max(desmat2[,2])),ylim=c(min(data[,2]),max(data[,2])),main="NON residential",xlab="sqrt(income/10000)",ylab="sqrt(popden/1000)")
hist(Montreal_popden[which(Montreal_ind==0)],xlim=c(min(Montreal_popden),max(Montreal_popden)),main="NON residential",xlab="popden")




#data = matrix(0,nrow=np,ncol=2)
#data[,1] = 1:np 
#data[,2] = log((Montreal_popden+1),10)
#
#
#desmat1 = matrix(1,nrow=np,ncol=1)
#desmat1[,1] =  Montreal_ind
#
#
#desmat2 = matrix(1,nrow=np,ncol=2)
#desmat2[,1] =  Montreal_ind
#desmat2[,2] =  log((Montreal_income+1),10)
#
#
#par(mfrow=c(2,2))
#plot(desmat2[which(Montreal_ind==1),2],data[which(Montreal_ind==1),2],main="residential",xlab="sqrt(income/10000)",ylab="sqrt(popden/1000)")
#hist(Montreal_popden[which(Montreal_ind==1)],xlim=c(min(Montreal_popden),max(Montreal_popden)),main="residential",xlab="popden",breaks=breaks)
#plot(desmat2[which(Montreal_ind==0),2],data[which(Montreal_ind==0),2],main="NON residential",xlab="sqrt(income/10000)",ylab="sqrt(popden/1000)")
#hist(Montreal_popden[which(Montreal_ind==0)],xlim=c(min(Montreal_popden),max(Montreal_popden)),main="NON residential",xlab="popden")




e = NULL

order=2

#  set up the FEM basis object

basisobj = create.FEM.basis(p, e, t, order)

#  set up a dummy FEM functional data object

Montrealfd = fd(numeric(basisobj$nbasis),basisobj)






loglambdaseq=seq(-4,2,by=1)

# Covar1

#GCVres_Covar1 = smooth.FEM.fd.Covar.GCV(data,desmat1,Montrealfd,loglambdaseq,CI_level=0.95)

loglambdaoptim_Covar1 = loglambdaseq[which.min(GCVres_Covar1$GCVseq)]
lambdaoptim_Covar1=10^loglambdaoptim_Covar1
loglambdaoptim_Covar1
min(GCVres_Covar1$GCVseq)

plot(loglambdaseq, GCVres_Covar1$GCVseq,main="Covar1")
abline(v=loglambdaoptim_Covar1)



# Covar2

#GCVres_Covar2 = smooth.FEM.fd.Covar.GCV(data,desmat2,Montrealfd,loglambdaseq,CI_level=0.95)

loglambdaoptim_Covar2 = loglambdaseq[which.min(GCVres_Covar2$GCVseq)]
lambdaoptim_Covar2=10^loglambdaoptim_Covar2
loglambdaoptim_Covar2
min(GCVres_Covar2$GCVseq)

plot(loglambdaseq, GCVres_Covar2$GCVseq,main="Covar2")
abline(v=loglambdaoptim_Covar2)


# NOCovar

#GCVres_NOCovar = smooth.FEM.fd.GCV(data,Montrealfd,loglambdaseq)

loglambdaoptim_NOCovar = loglambdaseq[which.min(GCVres_NOCovar$GCVseq)]
lambdaoptim_NOCovar=10^loglambdaoptim_NOCovar
loglambdaoptim_NOCovar
min(GCVres_NOCovar$GCVseq)

plot(loglambdaseq, GCVres_NOCovar$GCVseq,main="NOCovar")
abline(v=loglambdaoptim_NOCovar)


# # # # # # # # 
# Estimates

lambda=0.0005

# Covar1

#MontrealfelsplobjCovar1 = smooth.FEM.fd.Covar(data,desmat1,Montrealfd,lambdaoptim_Covar1)
MontrealfelsplobjCovar1 = smooth.FEM.fd.Covar(data,desmat1,Montrealfd,lambda)


fdobjCovar1 = MontrealfelsplobjCovar1$felsplobj

fhatCovar1=MontrealfelsplobjCovar1$felsplobj$coef
betahatCovar1= ( solve( t(desmat1) %*% desmat1 ) ) %*% t(desmat1) %*% (data[,2]-fhatCovar1[1:np,])
betahatCovar1

# Covar2

#MontrealfelsplobjCovar2 = smooth.FEM.fd.Covar(data,desmat2,Montrealfd,lambdaoptim_Covar2)
MontrealfelsplobjCovar2 = smooth.FEM.fd.Covar(data,desmat2,Montrealfd,lambda)

fdobjCovar2 = MontrealfelsplobjCovar2$felsplobj

fhatCovar2=MontrealfelsplobjCovar2$felsplobj$coef
betahatCovar2= ( solve( t(desmat2) %*% desmat2 ) ) %*% t(desmat2) %*% (data[,2]-fhatCovar2[1:np,])
betahatCovar2

#NOCovar

#MontrealfelsplobjNOCovar = smooth.FEM.fd(data,Montrealfd,lambdaoptim_NOCovar)
MontrealfelsplobjNOCovar = smooth.FEM.fd(data,Montrealfd,lambda)

fdobjNOCovar = MontrealfelsplobjNOCovar$felsplobj


# # # # # # # 
# PLOT
windows()
par(mfrow=c(2,2))
plot.FEM(fdobjCovar1)
title("Covar1")

plot.FEM(fdobjCovar2)
title("Covar2")

plot.FEM(fdobjNOCovar)
title("NOCovar")





#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################

### ANALYSIS OF RESIDUALS ###


#  evaluate the solution on the points

MontrealDataFitCovar1 = eval.FEM.fd(tractcenter[,1], tractcenter[,2], fdobjCovar1) +  desmat1 %*% betahatCovar1

MontrealDataFitCovar2 = eval.FEM.fd(tractcenter[,1], tractcenter[,2], fdobjCovar2) +  desmat2 %*% betahatCovar2

MontrealDataFitNOCovar = eval.FEM.fd(tractcenter[,1], tractcenter[,2],  fdobjNOCovar)



#  compute the residuals

MontrealResCovar1 = cbind(data[,1], data[,2] - MontrealDataFitCovar1)
min(MontrealResCovar1[,2])
max(MontrealResCovar1[,2])

MontrealResCovar2 = cbind(data[,1], data[,2] - MontrealDataFitCovar2)
min(MontrealResCovar2[,2])
max(MontrealResCovar2[,2])

MontrealResNOCovar = cbind(data[,1], data[,2] - MontrealDataFitNOCovar)
min(MontrealResNOCovar[,2])
max(MontrealResNOCovar[,2])




#  smooth the residuals

# ATTENZIONE: with what value of lambda?
lambda=0.001

ResMontrealfdCovar1 = smooth.FEM.fd(MontrealResCovar1,Montrealfd,lambda)

par(mfrow=c(2,2))
plot.FEM(ResMontrealfdCovar1$felsplobj)
title('Covar1')

ResMontrealfdCovar2 = smooth.FEM.fd(MontrealResCovar2,Montrealfd,lambda)

plot.FEM(ResMontrealfdCovar2$felsplobj)
title('Covar2')


ResMontrealfdNOCovar = smooth.FEM.fd(MontrealResNOCovar,Montrealfd,lambda)

plot.FEM(ResMontrealfdNOCovar$felsplobj)
title('NOCovar')






#  analysis of residuals

#  compute the squared distances and the residual products

dvec = numeric(np*(np-1)/2)
rvecCovar1  = dvec
rvecCovar2  = dvec
rvecNOCovar = dvec

m = 0
for (i in 2:np)
{   
    for (j in 1:(i-1))
    {  m = m + 1
       dvec[m] = sum((p[i,]-p[j,])^2)
       rvecCovar1[m] = MontrealResCovar1[i,2] * MontrealResCovar1[j,2]
       rvecCovar2[m] = MontrealResCovar2[i,2] * MontrealResCovar2[j,2]
       rvecNOCovar[m] = MontrealResNOCovar[i,2] * MontrealResNOCovar[j,2]
    }
}


windows()
par(mfrow=c(2,2))
plot(log10(dvec), rvecCovar1,xlab="log_{10} Squared distance",ylim=c(min(rvecCovar1,rvecCovar2,rvecNOCovar),max(rvecCovar1,rvecCovar2,rvecNOCovar)),main="Covar1")
abline(h=0,col=2,lwd=2,lty="dashed")
plot(log10(dvec), rvecCovar2,xlab="log_{10} Squared distance",ylim=c(min(rvecCovar1,rvecCovar2,rvecNOCovar),max(rvecCovar1,rvecCovar2,rvecNOCovar)),main="Covar2")
abline(h=0,col=2,lwd=2,lty="dashed")
plot(log10(dvec), rvecNOCovar,xlab="log_{10} Squared distance",ylim=c(min(rvecCovar1,rvecCovar2,rvecNOCovar),max(rvecCovar1,rvecCovar2,rvecNOCovar)),main="NO Covar")
abline(h=0,col=2,lwd=2,lty="dashed")



wilcox.test(rvecCovar1,rvecNOCovar,paired=T)

wilcox.test(rvecCovar2,rvecCovar1,paired=T)



#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################



X=NULL
Y=NULL

if (is.null(X))
{
    xmin = min(p[,1])
    xmax = max(p[,1])
    nx   = 201
    X    = matrix(seq(xmin, xmax, len=nx),ncol=1)
} else
{
    xmin = min(X)
    xmax = max(X)
    nx   = length(X)
}

if (is.null(Y))
{
    ymin = min(p[,2])
    ymax = max(p[,2])
    ny   = 201
    Y    = matrix(seq(ymin, ymax, len=ny),ncol=1)    
} else
{
    ymin = min(Y)
    ymax = max(Y)
    ny   = length(Y)
}




Xmat = X %*% matrix(1,nrow=1,ncol=ny)
Ymat = matrix(1,nrow=nx,ncol=1) %*% t(Y)
Xvec = NULL
for (numc in 1:nx)
{Xvec=c(Xvec,Xmat[,numc])}
Yvec = NULL
for (numc in 1:ny)
{Yvec=c(Yvec,Ymat[,numc])}



evalmatCovar1  = eval.FEM.fd(Xvec, Yvec, fdobjCovar1)
evalmatNOCovar = eval.FEM.fd(Xvec, Yvec, fdobjNOCovar)

evalmatiCovar1  = matrix(evalmatCovar1[,1] ,nrow=nx, ncol=ny, byrow=F)
evalmatiNOCovar = matrix(evalmatNOCovar[,1],nrow=nx, ncol=ny, byrow=F)



par(mfrow=c(2,2))
zlimit=c(min(c(evalmatiCovar1[which.min(evalmatiCovar1)],evalmatiNOCovar[which.min(evalmatiNOCovar)])),max(c(evalmatiCovar1[which.max(evalmatiCovar1)],evalmatiNOCovar[which.max(evalmatiNOCovar)])))

image(X,Y,evalmatiCovar1,col=heat.colors(100, alpha=1), zlim=zlimit,main="Covar1")
contour(X,Y,evalmatiCovar1,add=T, zlim=zlimit)


image(X,Y,evalmatiNOCovar,col=heat.colors(100, alpha=1), zlim=zlimit,main="NO Covar")
contour(X,Y,evalmatiNOCovar,add=T, zlim=zlimit)


# Now look at this!


Map <- GetMap.bbox(bb$lonR, bb$latR, centrale, destfile = 'Montreal_google-map.png',NEWMAP = FALSE, GRAYSCALE  = T)


PlotOnStaticMap(Map, add = FALSE, TrueProj=F,  FUN = points)  # background plot

image(X,Y,evalmatiCovar1,col=heat.colors(100, alpha=0.7), add = TRUE)
contour(X,Y,evalmatiCovar1,add=T)



for (nb in 1:100)
{
    points(p[Montreal_Tri_abconstrtot[nb,1:2],1],p[Montreal_Tri_abconstrtot[nb,1:2],2],type="l",lwd=2)
}

for (nb in 101:118)
{
    points(p[Montreal_Tri_abconstrtot[nb,1:2],1],p[Montreal_Tri_abconstrtot[nb,1:2],2],type="l",lwd=2)
}

for (nb in 119:133)
{
    points(p[Montreal_Tri_abconstrtot[nb,1:2],1],p[Montreal_Tri_abconstrtot[nb,1:2],2],type="l",lwd=2)
}




persp3d(as.vector(X),as.vector(Y),evalmati, xlab="", ylab="", zlab="", col='red', alpha=1)



plot(p[,1],p[,2],type="n",xlab="longitude",ylab="latidute")
polygon(c(p[t[1,1],1],p[t[1,2],1],p[t[1,3],1]),c(p[t[1,1],2],p[t[1,2],2],p[t[1,3],2]))
for (ne in 1:nt)
{polygon(c(p[t[ne,1],1],p[t[ne,2],1],p[t[ne,3],1]),c(p[t[ne,1],2],p[t[ne,2],2],p[t[ne,3],2]))
}










#####################################################################################################


# DIRICH


data[,2] = Montreal_popden/1000



source("SFDA_AllFunctions_MOD.R")

inodeconstr=(493+which(Montreal_bound_diric[1:133,1]==1))


order=1

e = NULL

basisobj = create.FEM.basis(p, e, t, order)

Montrealfd = fd(numeric(basisobj$nbasis),basisobj)


numnodes = dim(basisobj$params$nodes)[[1]]


inodefree=1:numnodes
inodefree=inodefree[-inodeconstr]

fconstrvalues=numeric(length(inodeconstr))

fbar=numeric(numnodes)
fbar[inodeconstr]=fconstrvalues


numnodefree=length(inodefree)



lambda=0.0001

# Covar1

MontrealfelsplobjCovar1 = smooth.FEM.fd.Covar(data,desmat1,Montrealfd,lambda)

GCVres_Covar1 = smooth.FEM.fd.Covar.GCV(data,desmat1,Montrealfd,log(lambda,10),CI_level=0.95)


fdobjCovar1 = MontrealfelsplobjCovar1$felsplobj

fhatCovar1=MontrealfelsplobjCovar1$felsplobj$coef
betahatCovar1= ( solve( t(desmat1) %*% desmat1 ) ) %*% t(desmat1) %*% (data[,2]-fhatCovar1[1:np,])
betahatCovar1




X=NULL
Y=NULL

if (is.null(X))
{
    xmin = min(p[,1])
    xmax = max(p[,1])
    nx   = 201
    X    = matrix(seq(xmin, xmax, len=nx),ncol=1)
} else
{
    xmin = min(X)
    xmax = max(X)
    nx   = length(X)
}

if (is.null(Y))
{
    ymin = min(p[,2])
    ymax = max(p[,2])
    ny   = 201
    Y    = matrix(seq(ymin, ymax, len=ny),ncol=1)    
} else
{
    ymin = min(Y)
    ymax = max(Y)
    ny   = length(Y)
}


Xmat = X %*% matrix(1,nrow=1,ncol=ny)
Ymat = matrix(1,nrow=nx,ncol=1) %*% t(Y)
Xvec = NULL
for (numc in 1:nx)
{Xvec=c(Xvec,Xmat[,numc])}
Yvec = NULL
for (numc in 1:ny)
{Yvec=c(Yvec,Ymat[,numc])}

evalmatCovar1  = eval.FEM.fd(Xvec, Yvec, fdobjCovar1)

evalmatiCovar1  = matrix(evalmatCovar1[,1] ,nrow=nx, ncol=ny, byrow=F)



#Map <- GetMap.bbox(bb$lonR, bb$latR, centrale, destfile = 'Montreal_google-map.png',NEWMAP = FALSE, GRAYSCALE  = T)

windows()
PlotOnStaticMap(Map, add = FALSE, TrueProj=F,  FUN = points)  # background plot

image(X,Y,evalmatiCovar1,col=heat.colors(100, alpha=0.7), add = TRUE)
contour(X,Y,evalmatiCovar1,add=T)



for (nb in 1:100)
{
    points(p[Montreal_Tri_abconstrtot[nb,1:2],1],p[Montreal_Tri_abconstrtot[nb,1:2],2],type="l",lwd=2)
}

for (nb in 101:118)
{
    points(p[Montreal_Tri_abconstrtot[nb,1:2],1],p[Montreal_Tri_abconstrtot[nb,1:2],2],type="l",lwd=2)
}

for (nb in 119:133)
{
    points(p[Montreal_Tri_abconstrtot[nb,1:2],1],p[Montreal_Tri_abconstrtot[nb,1:2],2],type="l",lwd=2)
}




persp3d(as.vector(X),as.vector(Y),evalmatiCovar1, xlab="", ylab="", zlab="", col='red', alpha=1)




MontrealDataFitCovar1 = eval.FEM.fd(tractcenter[,1], tractcenter[,2], fdobjCovar1) +  desmat1 %*% betahatCovar1
MontrealResCovar1 = cbind(data[,1], data[,2] - MontrealDataFitCovar1)
min(MontrealResCovar1[,2])
max(MontrealResCovar1[,2])