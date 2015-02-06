load("FrontieraC.RData")
library(mgcv)
library(SDMTools)
library(MASS)

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


# Devo aggiungere del rumore al dato?
noise<-TRUE

# Per quanti istanti di tempo eseguire la stima?
TimePoints<-0:5

xknot<-x[TypePoint]
yknot<-y[TypePoint]
xbound<-x[!TypePoint]
ybound<-y[!TypePoint]

fun=function(x,y,t,xbound,ybound)
{
    Bound<-cbind(xbound,ybound)
    Result<-numeric(length(x))
    
    for(i in 1:length(x))
    {
        if(pnt.in.poly(cbind(x[i],y[i]),Bound)$pip==0)
        {
            Result[i]=NA
        } else
        {
            Result[i]=(x[i]+y[i]+t[i])^2 + log(t[i]+1)
        }
    }
    
    return(Result)
}

Data<-NULL
TimeData<-NULL
xData<-NULL
yData<-NULL

for(i in TimePoints)
{
    Time<-rep(i,length(x[TypePoint]))
    if(noise)
    {
        Data<-c(Data,(fun(xknot,yknot,Time,xbound,ybound)+rnorm(length(xknot),0,0.005)))
    } else
    {
        Data<-c(Data,fun(xknot,yknot,Time,xbound,ybound))
    }
    TimeData<-c(TimeData,Time)
    xData<-c(xData,xknot)
    yData<-c(yData,yknot)
}

#### full model with vessel as random effect
mod <-gam(Data~s(TimeData,bs="cr")+te(xData,yData,year,d=c(2,1),bs=c("so","so"),k=c(30,5)),
            knots=cbind(xData,yData),data=Data,family=Tweedie(p=0,link="log"),method="REML")
printresults(mod1a)