load("FrontieraC.RData")
library(mgcv)

xbound=x[!TypePoint]
ybound=y[!TypePoint]
NX=400
NY=200
# Genero i dati per il grafico
xGen<-seq(-1,3.5,length.out=NX)
yGen<-seq(-1,1,length.out=NY)

xPlot<-NULL
yPlot<-NULL
PosMatrix<-NULL
RealValuesPlot<-NULL

for(i in 1:NX)
{
    for(j in 1:NY)
    {
        if(!is.na(fs.test(xGen[i],yGen[j])))
        {
            xPlot<-c(xPlot,xGen[i])
            yPlot<-c(yPlot,yGen[j])
            PosMatrix<-rbind(PosMatrix,c(i,j))
            RealValuesPlot<-c(RealValuesPlot,fs.test(xGen[i],yGen[j]))
        }
    }
}

PlotMatrix=matrix(nrow=NX,ncol=NY)
for(i in 1:(dim(PosMatrix)[1]))
{
    PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=RealValuesPlot[i]
}
png(filename="DomC_fstest")
image(xGen,yGen,PlotMatrix,col=heat.colors(100),main="Funzione spaziale",xlab="",ylab="")
lines(xbound,ybound,lwd=3)
lines(c(xbound[1],xbound[length(xbound)]),c(ybound[1],ybound[length(ybound)]),lwd=3)
contour(xGen,yGen,PlotMatrix,nlevels=10,add=TRUE)
dev.off()
