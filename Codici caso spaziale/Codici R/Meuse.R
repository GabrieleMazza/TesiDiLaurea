### MEUSE data
### Exploratory analysis described in
### Bivand-Pebesma-GómezRubio, Applied Spatial Data Analysis with R, Springer 2008 


# 30 June 2010

library(lattice)
library(sp)



data(meuse)
coordinates(meuse) <- c("x", "y")
spplot(meuse, "zinc", do.log = T)
bubble(meuse, "zinc", do.log = T, key.space = "bottom")

## zinc concentration is larger close to the
## river Meuse banks

xyplot(log(zinc) ~ sqrt(dist), as.data.frame(meuse))
zn.lm <- lm(log(zinc) ~ sqrt(dist), meuse)

meuse$fitted.s <- predict(zn.lm, meuse) - mean(predict(zn.lm, meuse))
meuse$residuals <- residuals(zn.lm)
spplot(meuse, c("fitted.s", "residuals"))



###########################
###### INTERPOLATION ######
###########################

## interpolation on the regular grid given by the meuse.grid data.frame
## then converted into a SpatialPixelsDataFrame

data(meuse.grid)
coordinates(meuse.grid) <- c("x", "y")
plot(coordinates(meuse.grid))
meuse.grid <- as(meuse.grid, "SpatialPixelsDataFrame")


## Inverse Distance Weighted Interpolation
##
## Inverse distance-based weighted interpolation (IDW) computes a weighted average
## where weights w for observations are computed according to their distance to
## the interpolation location
##
## w(s_i) = ||s_i - s_0||^{-p}
##
## default p = 2


library(gstat)
idw.out <- idw(zinc ~ 1, meuse, meuse.grid, idp = 2.5)

spplot(idw.out)



## Linear Regression


zn.lm <- lm(log(zinc) ~ sqrt(dist), meuse)
meuse.grid$pred <- predict(zn.lm, meuse.grid)
meuse.grid$se.fit <- predict(zn.lm, meuse.grid, se.fit = TRUE)$se.fit    # standard error of predicted means

# Alternatively, using function krige in gstat
 
meuse.lm <- krige(log(zinc) ~ sqrt(dist), meuse, meuse.grid)
spplot(meuse.lm)

## !!! in this case we aren't doing any kriging since no variogram is specified
## but we are simply using linear regression.
## Result identical to that of lm:
max(abs(meuse.lm$var1.pred-meuse.grid$pred))


## polynomial regression 
meuse.tr2 <- krige(log(zinc) ~ 1, meuse, meuse.grid, degree = 2)   # second-order polynomial of spatial coordinates
spplot(meuse.tr2)





########################################
#### Estimating Spatial Correlation ####
########################################


## Exploratory Variogram Analysis

## A simple way to acknowledge that spatial correlation is present or not is
## to make scatter plots of pairs Z(s_i) and Z(s_j), grouped according to their
## separation distance h_ij = ||s_i - s_j||

hscat(log(zinc) ~ 1, meuse, (0:9) * 100)

## A second way to explore spatial correlation is by plotting the variogram
## and the variogram cloud. The variogram cloud is obtained by plotting all
## possible squared differences of observation pairs (Z(s_i)-Z(s_j))^2 
## against their separation distance h_ij

plot(variogram(log(zinc) ~ 1, meuse, cloud = TRUE))           # it assumes isotropy

plot(variogram(log(zinc) ~ 1, meuse))

variogram(log(zinc) ~ 1, meuse,cressie=TRUE)                  # uses Cressie's robust variogram estimate

plot(variogram(log(zinc) ~ 1, meuse)$dist,variogram(log(zinc) ~ 1, meuse)$gamma)
points(variogram(log(zinc) ~ 1, meuse,cressie=TRUE)$dist,variogram(log(zinc) ~ 1, meuse,cressie=TRUE)$gamma,col="red") # robust empirical variogram in red





# Is isotropy a reasonable assumption?

plot(variogram(log(zinc) ~ 1, meuse, alpha = c(0, 45, 90, 135))) # Directional sample variogram 
                                                                 # for four directions (0 is North, 90 is East)
                                                                 # It work dividind data paris in direction groups, e.g.,
                                                                 # any point pair between 22.5deg and 67.5deg is used for 
                                                                 # the 45deg panel

# Isotropy does not seem reasonable!!  (see later on possible remedies) 
                                                               
                                                                 
# Changing cut-off and width in the computation of the empirical variogram

plot(variogram(log(zinc) ~ 1, meuse, cutoff = 1000, width =50))  # default cutoff value is one third of the largest diagonal 
                                                                 # of the bounding box
                                                                 # default width is cutoff value divided by 15
                                   

plot(variogram(log(zinc) ~ 1, meuse, boundaries = c(0, 50, 100, seq(250, 1500, 250))))
                                                                 # irregular intervals, e.g. to zoom in on the short 
                                                                 # distance variogram without revealing irrelevant details 
                                                                 # for the longer distances





###################################
####### Variogram Modelling #######
###################################


# Overview of the basic variogram models (gstat)
show.vgms()

show.vgms(model = "Mat", kappa.range = c(0.1, 0.2, 0.5, 1, 2, 5, 10), max = 10) # Matern class

vgm() # list of implemented variogram models

# valid variogram models are constructed by using one or combining two or more basic variogram models
# e.g.

vgm(1, "Sph", 300)

vgm(1, "Sph", 300, 0.5)

v1 <- vgm(1, "Sph", 300, 0.5)

v2 <- vgm(0.8, "Sph", 800, add.to = v1)

v2


## Fitting a variogram model by non-linear least squares optimization, using the empirical variogram

v <- variogram(log(zinc) ~ 1, meuse)
plot(v)                                        # the spherical model looks like a reasonable choice

v.fit=fit.variogram(v, vgm(1, "Sph", 800, 1))  # Initial values for the variogram fit are needed for fit.variogram
v.fit


fit.variogram(v, vgm(1, "Sph", 10, 1))         # bad specification of initial parameter values 


attr(v.fit, "SSErr")


## Fitting a variogram model by REML (restricted maximum likelihood)
#  Gaussian random field
#  Does not use empirical variogram

fit.variogram.reml(log(zinc) ~ 1, meuse, model = vgm(0.6,"Sph", 800, 0.06))




# Modelling anisotropy by defining a range ellipse instead of a circular or spherical range

v.dir <- variogram(log(zinc) ~ 1, meuse, alpha = (0:3) * 45)
v.anis <- vgm(0.6, "Sph", 1600, 0.05, anis = c(45, 0.3))  # anis parameters define an anisotropy ellipse.
                                                          # The first parameter, 45, refers to the main axis direction: 
                                                          # it is the angle for the principal direction of continuity 
                                                          # (measured in degrees, clockwise from positive Y, i.e. North). 
                                                          # The second parameter, 0.3, is the anisotropy ratio, the ratio 
                                                          # of the minor range to the major range (a value between 0 and 1)
plot(v.dir, v.anis)



# Variogram map

plot(variogram(log(zinc) ~ 1, meuse, map = TRUE, cutoff = 1000,width = 100))

                                                          # bins h vectors in square grid cells over x and y





## Residual Variogram Modelling

# Residual variograms are calculated by default when a more complex model for the trend is used.
# Ordinary least squares residuals are used
# for estimating the trend, with observations considered as independent

v=variogram(log(zinc) ~ sqrt(dist), meuse)
plot(v)

# To honour a dependence structure present, generalised least squares residuals can be calculated instead.
# For this, a variogram model to define the covariance structure is needed, e.g.

f <- log(zinc) ~ sqrt(dist)
vt <- variogram(f, meuse)
vt.fit <- fit.variogram(vt, vgm(1, "Exp", 300, 1))
vt.fit


g.wls <- gstat(NULL, "log-zinc", f, meuse, model = vt.fit, set = list(gls = 1))
vg.wls=variogram(g.wls)
plot(vg.wls)

# for this specific example, very little difference
(variogram(g.wls)$gamma - vt$gamma)/mean(vt$gamma)




######################################
######### SPATIAL PREDICTION #########
######################################


## Universal, Ordinary, and Simple Kriging

lz.sk <- krige(log(zinc) ~ 1, meuse, meuse.grid, v.fit, beta = 5.9)   #[using simple kriging]
spplot(lz.sk)


lz.ok <- krige(log(zinc) ~ 1, meuse, meuse.grid, v.fit)               #[using ordinary kriging]
windows()
spplot(lz.ok)

lz.uk <- krige(log(zinc) ~ sqrt(dist), meuse, meuse.grid, vt.fit)     #[using universal kriging]
windows()
spplot(lz.uk)


## Block kriging

# Ordinary block kriging

lz.ok <- krige(log(zinc) ~ 1, meuse, meuse.grid, v.fit, block = c(40, 40))  # Ordinary block kriging 
spplot(lz.ok)                                                               # for blocks of size 40 × 40
                                                                            # [using ordinary kriging]


 
xy <- expand.grid(x = seq(-20, 20, 4), y = seq(-20, 20, 4))                 # For a circular shape with radius 20, 
xy <- xy[(xy$x^2 + xy$y^2) <= 20^2, ]                                       # centred on the points of meuse.grid, one
                                                                            # could select points on a regular grid within 
                                                                            # a circle
lz.ok <- krige(log(zinc) ~ 1, meuse, meuse.grid, v.fit, block = xy)
spplot(lz.ok) 





####################################
####### Model Diagnostics ##########
####################################



# Cross Validation Residuals
 
# We split the dataset in 100 observations for modelling and 55 for testing

sel100 <- sample(1:155, 100)   
m.model <- meuse[sel100, ]
m.valid <- meuse[-sel100, ]

v100.fit <- fit.variogram(variogram(log(zinc) ~ 1, m.model), vgm(1, "Sph", 800, 1))

m.valid.pr <- krige(log(zinc) ~ 1, m.model, m.valid, v100.fit)

resid.kr <- log(m.valid$zinc) - m.valid.pr$var1.pred

summary(resid.kr)

resid.mean <- log(m.valid$zinc) - mean(log(m.valid$zinc))

R2 <- 1 - sum(resid.kr^2)/sum(resid.mean^2) # kriging prediction is a better predictor than the mean
                                            # with an indicative R2 of
R2 


m.valid.pr$res <- resid.kr
bubble(m.valid.pr, "res")                   # map of cross validation residuals

# more on cross-validation on the book








####################
## Spatial CLASSES

library(sp)
getClass("Spatial")
getClass("CRS")

m <- matrix(c(0, 0, 1, 1), ncol = 2, dimnames = list(NULL, c("min", "max")))
crs <- CRS(projargs = as.character(NA))

crs

S <- Spatial(bbox = m, proj4string = crs)


getClass("SpatialPoints")








##############################################
##############################################
##############################################
##############################################
##############################################
##############################################


library(lattice)
library(sp)
library(gstat)



data(meuse)
coordinates(meuse) <- c("x", "y")


## Residual Variogram Modelling

# Residual variograms are calculated by default when a more complex model for the trend is used.
# Ordinary least squares residuals are used
# for estimating the trend, with observations considered as independent


f=log(zinc) ~ sqrt(dist)+elev

v=variogram(f, meuse)
plot(v)

# To honour a dependence structure present, generalised least squares residuals can be calculated instead.
# For this, a variogram model to define the covariance structure is needed, e.g.

vt.fit <- fit.variogram(v, vgm(1, "Exp", 300, 1))
vt.fit


g.wls <- gstat(NULL, "log-zinc", f, meuse, model = vt.fit, set = list(gls = 1))
vg.wls=variogram(g.wls)
plot(vg.wls)

# the difference is
(vg.wls$gamma - v$gamma)/mean(v$gamma)




######################################
######### SPATIAL PREDICTION #########
######################################


## Universal Kriging

data(meuse.grid)
coordinates(meuse.grid) <- c("x", "y")
plot(coordinates(meuse.grid))
meuse.grid <- as(meuse.grid, "SpatialPixelsDataFrame")



####
coord.augmented=rbind(coordinates(meuse),coordinates(meuse.grid))
dist.augmented=c(meuse$dist,meuse.grid$dist)
meuse.grid.augmented=data.frame(x=coord.augmented[,1], y=coord.augmented[,2],dist=dist.augmented)
coordinates(meuse.grid.augmented) <- c("x", "y")
plot(coordinates(meuse.grid.augmented))
meuse.grid.aug <- as(meuse.grid.augmented, "SpatialPixelsDataFrame")     
#### it does not like it....





lz.uk <- krige(f, meuse, meuse, vt.fit)     #[using universal kriging]
windows()
spplot(lz.uk)
resid.kr <- log(meuse$zinc) - lz.uk$var1.pred
resid.kr 





####################################
####### Model Diagnostics ##########
####################################



# Cross Validation Residuals
 
# We split the dataset in 100 observations for modelling and 55 for testing



np=155
ressave=numeric(np)
    for (idrop in 1:np)
     {  m.model <- meuse[-idrop, ]
        m.valid <- meuse[idrop, ]
        v154.fit <- fit.variogram(variogram(f, m.model), vgm(1, "Sph", 800, 1))
#g.wls <- gstat(NULL, "log-zinc", f, meuse, model = vt.fit, set = list(gls = 1))
#vg.wls=variogram(g.wls)
        m.valid.pr <- krige(f, m.model, m.valid, v154.fit)
        ressave[idrop] <- log(m.valid$zinc) - m.valid.pr$var1.pred
     }
    SE = sqrt(mean(ressave^2))
    
    
    SE












fitmeuse=Krig(x=meuse[,1:2],Y=log(meuse$zinc),Z=cbind(sqrt(meuse$dist),meuse$elev),m=1)

predict(fitmeuse)-log(meuse$zinc)

max(abs(fitmeuse$fitted.values.null- fitmeuse$d[1]- fitmeuse$d[2]*sqrt(meuse$dist)- fitmeuse$d[3]*meuse$elev))
 
fitmeuse$d

np=155
ressave=numeric(np)
    for (idrop in 1:np)
     {  m.model <- meuse[-idrop, ]
        m.valid <- meuse[idrop, ]
        m.valid.pr <- krige(f, m.model, m.valid, v154.fit)
        ressave[idrop] <- log(m.valid$zinc) - m.valid.pr$var1.pred
     }
    SE = sqrt(mean(ressave^2))
    
    
    SE
