# Supporting material - Example R Code for models presented in:
#
# Nicole H. Augustin, Verena M. Trenkel, Simon N. Wood and Pascal Lorance.
# Space-time modelling of blue ling for fisheries stock management.
# Environmetrics.
# The data is confidential and cannot be shared.
# preliminaries
# attach required libraries
.libPaths("/home/nha20/Rlocal")
load("Dat.rda") ### blue ling data
load("conBL.rda") ### data frame with polygons for boundary,
### variables longitude and latitude
knots<-read.table("Knots.csv",header=T,sep=";")
## file with variables longitude and latitude of 74 knot positions
knots[12,2] <- 6390 ## one knot is on boundary hence it is moved
library(mgcv)
library(soap)
library(MASS)
### function to print model summaries
printresults<-function(gobject){
    edf=sum(gobject$edf)
    nobs=length(gobject$y)
    BIC<--2*logLik(gobject)+edf*log(nobs)
    results<-c(sum(gobject$edf), gobject$deviance,summary(gobject)$r.sq,
               gobject$aic,BIC)
    results<-round(results,2)
    names(results)<-c("edf","deviance","adj.R-sq","AIC","BIC")
    results
}
dknots<- nrow(knots)
#Model 1a
#### full model with vessel as random effect
Dat$vessel<- factor(Dat$vessel)
mod1a <-gam(blue_ling~s(vessel,bs="re")+s(duration)+s(depth)+s(month,bs="cc")+
                te(depth,month,bs=c("tp","cc"),k=c(10,10))
            +te(depth,year,bs=c("tp","cr"),k=c(10,7))
            +te(longitude,latitude,year,d=c(2,1),bs=c("sf","cr"),k=c(30,5),
                xt=list(list(bnd=list(conBL)),NULL))+te(longitude,latitude,year,d=c(2,1),
                                                        bs=c("sw","cr"),k=c(dknots,5),xt=list(list(bnd=list(conBL)),NULL)),
            knots=knots,data=Dat,family=Tweedie(p=1.5,link="log"),method="REML")
printresults(mod1a)
#Model 2a
#FULL model with depth (*year+month) + lat*long*year , fishing power
mod2a<-gam(blue_ling~s(duration)+s(depth)+s(month,bs="cc")+
               te(depth,month,bs=c("tp","cc"),k=c(10,10))+te(depth,year,bs=c("tp","cr"),
                                                             k=c(10,7))
           +te(longitude,latitude,year,d=c(2,1),bs=c("sf","cr"),k=c(30,5),
               xt=list(list(bnd=list(conBL)),NULL))+
               te(longitude,latitude,year,d=c(2,1),bs=c("sw","cr"),k=c(dknots,5),
                  xt=list(list(bnd=list(conBL)),NULL))+power,
           knots=knots,data=Dat,family=Tweedie(p=1.5,link="log"),method="REML")
#Model 3a
#### additive model with vessel as random effect
mod3a <-gam(blue_ling~s(vessel,bs="re")+s(duration)+s(depth)+s(month,bs="cc")+
                te(depth,month,bs=c("tp","cc"),k=c(10,10))+te(depth,year,bs=c("tp","cr"),
                                                              k=c(10,7))
            +s(longitude,latitude,bs="so", k=30,xt=list(bnd=list(conBL)))+
                s(year,bs="cr",k=5),knots=knots,data=Dat,family=Tweedie(p=1.5,link="log"),
            method="REML")
#Model 4a
###soap additive with power
mod4a<-gam(blue_ling~s(duration)+s(depth)+s(month,bs="cc")+
               te(depth,month,bs=c("tp","cc"),k=c(10,10))+te(depth,year,bs=c("tp","cr"),
                                                             k=c(10,7))+ s(year,bs="cr",k=5)+
               s(longitude,latitude,bs="so",k=dknots,xt=list(bnd=list(conBL)))+power,
           knots=knots,data=Dat,family=Tweedie(p=1.5,link="log"),method="REML")
#Model 1b
#### full model with vessel as random effect - NOSOAP
mod1b <-gam(blue_ling~s(vessel,bs="re")+s(duration)+s(depth)+s(month,bs="cc")+
                te(depth,month,bs=c("tp","cc"),k=c(10,10))+te(depth,year,bs=c("tp","cr"),
                                                              k=c(10,7))
            +te(longitude,latitude,year,d=c(2,1),bs=c("tp","cr"),k=c(dknots+30,5)),
            data=Dat,family=Tweedie(p=1.5,link="log"),method="REML")
#Model 2b
#### FULL: in comparison without soap smooth
mod2b<-gam(blue_ling~s(duration)+s(depth)+s(month,bs="cc")+
               te(depth,month,bs=c("tp","cc"),k=c(10,10))+te(depth,year,bs=c("tp","cr"),
                                                             k=c(10,7))
           +te(longitude,latitude,year,d=c(2,1),bs=c("tp","cr"),k=c(dknots+30,5))
           +power,
           data=Dat,family=Tweedie(p=1.5,link="log"),method="REML")
#Model 3b
#### additive model with vessel as random effect - NOSOAP
mod3b<-gam(blue_ling~s(vessel,bs="re")+s(duration)+s(depth)+s(month,bs="cc")+
               te(depth,month,bs=c("tp","cc"),k=c(10,10))+te(depth,year,bs=c("tp","cr"),
                                                             k=c(10,7))
           +te(longitude,latitude,d=c(2),bs=c("tp"),k=c(dknots+30))
           +s(year,bs="cr",k=5),
           data=Dat,family=Tweedie(p=1.5,link="log"),method="REML")
#Model 4b
#### additive model with vessel power - NOSOAP
mod4b<-gam(blue_ling~power+s(duration)+s(depth)+s(month,bs="cc")+
               te(depth,month,bs=c("tp","cc"),k=c(10,10))+te(depth,year,bs=c("tp","cr"),
                                                             k=c(10,7))
           +te(longitude,latitude,d=c(2),bs=c("tp"),k=c(dknots+30))
           +s(year,bs="cr",k=5),
           data=Dat,family=Tweedie(p=1.5,link="log"),method="REML")
########### Computing time trends with Bayesian credible intervals
gobject<-mod1a
### obtain design/prediction matrix
ndat<-Dat
M<-predict(gobject,newdata=ndat,type="lpmatrix")
### ndat is the data with the spatial grid we want to predict for
### here we have set it to data as observed, for prediction need to
### need to set values in ndat as required
#### simulate 1000 fitted values from posterior distribution of parameters
simcoef <-mvrnorm(n=1000, coef(gobject), gobject$Vp)
simfit<-as.matrix(M)%*% t(simcoef)
simfit <- aggregate(simfit, by=list(ndat$year),mean,na.rm=TRUE)
years<-simfit[,1]
simfit<-simfit[,-1] ### exclude group index
## obtain quantiles for trendplot
simquant<-apply(simfit,1,quantile,p=c(0.025,0.5,0.975),na.rm=TRUE)
## basic time trend plot
plot(years,simquant[2,],type="o",pch=1,lty=1,col=1)
lines(years,simquant[1,],lty=2,col=1)
lines(years,simquant[3,],lty=2,col=1)
#### plot of spatial predictions
years<- sort(unique(Dat$year))
zlimmi<-range(gobject$linear.predictors)
par(mfrow=c(3,3))
par(pty="s")
oldpar<-par()
par(mar=oldpar$mar-c(4,3,3,2))
for (i in 1:length(years)) {
    print(years[i])
    vis.gam(gobject,view=c("longitude","latitude"),zlim=zlimmi,
            cond=list(year=years[i]),plot.type="contour",type="link",
            too.far=0.05, color="topo", main=" ", xlab="",ylab="")
    text(x=600,y=6200,paste(years[i]))
    points(Dat$longitude[Dat$year==years[i]],Dat$latitude[Dat$year==years[i]]
           ,pch=14,cex=0.3)
    lines(conBL$longitude,conBL$latitude,type="l",xlim=range(Dat$longitude),
          ylim=range(Dat$latitude)) ## add points where sample locations are
}