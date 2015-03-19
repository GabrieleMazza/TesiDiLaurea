library(rgl)
library(mgcv)
library(SDMTools)
library(ggplot2)
source("SpazioTempo.R")
load("FrontieraC.RData")
load("Validazione.RData")

set.seed(10)
# Devo aggiungere del rumore al dato?
noise<-TRUE

SigmaNoise<-0.5

# Per quanti istanti di tempo eseguire la stima?
TimePoints<-seq(0,2*pi,length.out=9)

# Distinguo i punti interni dai punti di bordo
xknot<-x[TypePoint]
yknot<-y[TypePoint]
xbound<-x[!TypePoint]
ybound<-y[!TypePoint]

# Function for perturbation
fun=function(x,y,t)
{
    y=fs.test(x,y)*cos(t)
    return(y)
}

# Tutti i punti della validazione sono contenuti all'interno del bordo?
Bound<-cbind(xbound,ybound)
sum(pnt.in.poly(cbind(xValid,yValid),Bound)$pip)==length(xValid)




##### OGGETTI UTILI #####

TimeBasisObj<-Create.Bspline.Time.Basis(TimePoints,TimeOrder=4,DerivativeOrder=2,PlotIt=F)
SpaceBasisObj<-Create.FEM.Space.Basis(cbind(x,y),Triang,TypePoint,1)



##### DATI E BUBBLEPLOT DEI DATI #####

# Genero una volta i dati allo stesso modo di come farò dopo nel ciclo,
# per vederne i bubbleplot...

# Ho bisogno dei un dataframe
Knot=data.frame(xknot=c(xknot,10,10),yknot=c(yknot,1,-1))
Boundary<-data.frame(xbound,ybound)
# Genero i dati
DataMatrix<-NULL
for(t in TimePoints)
{
    if(noise)
    {
        DataMatrix<-cbind(DataMatrix,fun(xknot,yknot,rep(t,length(xknot)))+rnorm(length(xknot),0,SigmaNoise))
    } else
    {
        DataMatrix<-cbind(DataMatrix,fun(xknot,yknot,rep(t,length(xknot))))
    }
}

for(j in 1: length(TimePoints))
{
    Gen<-c(DataMatrix[,j],min(DataMatrix),max(DataMatrix))
    
    ggplot(data=Knot, aes(x=xknot, y=yknot)) +
        geom_point(aes(size=(Gen),color=Gen)) +
        scale_color_gradient(low="yellow", high="red") +
        geom_polygon(data=Boundary,aes(x=xbound,y=ybound),alpha=1,colour="black", fill=NA, size=1.1) +
        theme_bw() +
        labs(color="Dato")+
        scale_x_continuous(name="x", limits=c(-1,3.5)) +
        scale_y_continuous(name="y",limits=c(-1,1)) +
        labs(title=paste("Dati Generati al tempo ",round(TimePoints[j],2),sep="")) + 
        theme(plot.title = element_text(size = rel(1.5), face="bold")) +
        guides(size=FALSE)
    ggsave(file=paste("Dati Generati al tempo ",round(TimePoints[j],2),".png",sep=""))
}

rm(Knot,Boundary)

dev.off()

Data<-NULL
for(i in 1:(dim(DataMatrix)[1]))
{
    Data<-c(Data,DataMatrix[i,])
}

##### GCV #####

# LogS = seq(-1,0,by=0.125)
# LogT = seq(-3.75,-2.75,by=0.125)
# GCVResult<-ST.GCV(Data,SpaceBasisObj,TimeBasisObj,LogS,LogT)
# png(filename="GCV Matrix.png")
# image(LogS,LogT,GCVResult$GCVMatrix,col=heat.colors(100),main="GCV Matrix",xlab="logLambdaS",ylab="logLambdaT")
# dev.off()
# 
# save(file="GCVResult.RData",GCVResult,LogS,LogT)
# 
# LambdaS=10^GCVResult$Best[1]
# LambdaT=10^GCVResult$Best[2]

LambdaS=10^-0.375
LambdaT=10^-3.250




##### SOLUZIONE #####

##### RISOLUZIONE DEL SISTEMA #####

SolutionObj<-ST.Smooth(Data,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)

# Ora salvo i risultati
write.table(SolutionObj$C,file="VettoreC.txt",row.names=FALSE,col.names=FALSE)


##### GRAFICI ANNO PER ANNO SOLO DELLA FUNZIONE SENZA COVARIATE #####

PlotMatrix=matrix(ncol=(length(yvec)+1),nrow=length(xvec))
# E' stata aggiunta una colonna, che sarà utile per avere massimi e minimi uniformi
# tra tutti i plot

# Voglio ricreare perfettamente le posizioni della image matrix tru

# RICERCA DEL MASSIMO
# PER AVERE UNA SCALA DI COLORI PER IL GRAFICO UNIFORME IN TUTTI GLI ISTANTI DI TEMPO

ResultFitted<-NULL

for(j in TimePoints)
{
    Time<-rep(j,length(xValid))
    Result<-ST.Eval(xValid,yValid,Time,SolutionObj)
    
    ResultFitted<-cbind(ResultFitted,Result)
}

ResultReal<-NULL

for(j in TimePoints)
{
    Time<-rep(j,length(xValid))
    Result<-fun(xValid,yValid,Time)
    ResultReal<-cbind(ResultReal,Result)
}

save(file="ResultFitted.RData", ResultFitted)
Max=max(ResultFitted,ResultReal,na.rm=TRUE)
Min=min(ResultFitted,ResultReal,na.rm=TRUE)
zlim<-c(Min,Max)

#AGGIUSTO LA COLONNA IN PIU'... QUI AGGIUNGO IL MASSIMO E IL MINIMO
PlotMatrix[1,(length(yvec)+1)]=Max
PlotMatrix[2,(length(yvec)+1)]=Min

# GRAFICI

# Ora salvo tutti i grafici
# Voglio ricreare perfettamente le posizioni della image matrix tru

nx<-length(xvec)
ny<-length(yvec)
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))

for(j in 1:length(TimePoints))
{
    
    # IMP! PER UNIFORMARE LA SCALA SPAZIALE, TRACCIO LONTANO LA COLONNA CON IL MASSIMO
    # E IL MINIMO, E POI LA TAGLIO CON ylim
    
    
    # Funzione Reale
    for(i in 1:(dim(PosMatrix)[1]))
    {
        PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=ResultReal[i,j]
    }
    # Plot
    png(filename=paste("Tempo ",round(TimePoints[j],2)," reale.png",sep=""))
    image(xvec,c(yvec,4),PlotMatrix,col=heat.colors(100),main=paste("Funzione reale tempo ",round(TimePoints[j],2),sep=""),ylim=c(-1,+1),xlab=" ",ylab=" ")
    lines(xbound,ybound,lwd=3)
    lines(c(xbound[1],xbound[length(xbound)]),c(ybound[1],ybound[length(ybound)]),lwd=3)
    contour(xvec,c(yvec,4),PlotMatrix,nlevels=10,add=TRUE,lwd=2,labcex=1.1)
    dev.off()
    
    # Matrice con la stimata
    for(i in 1:(dim(PosMatrix)[1]))
    {
        PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=ResultFitted[i,j]
    }
    # Plot
    png(filename=paste("Tempo ",round(TimePoints[j],2)," stimata.png",sep=""))
    image(xvec,c(yvec,4),PlotMatrix,col=heat.colors(100),main=paste("Funzione stimata tempo ",round(TimePoints[j],2),sep=""),ylim=c(-1,+1),xlab=" ",ylab=" ")
    lines(xbound,ybound,lwd=3)
    lines(c(xbound[1],xbound[length(xbound)]),c(ybound[1],ybound[length(ybound)]),lwd=3)
    contour(xvec,c(yvec,4),PlotMatrix,nlevels=10,add=TRUE,lwd=2,labcex=0.9)
    dev.off()
}

i<-98
xP<-xknot[i]
yP<-yknot[i]
png(filename=paste("Plot per un punto fissato",fs.test(xP,yP),".png",sep=" "))
FixedPointPlot(xP,yP,SolutionObj,lwd=2,NameLocation = paste("(",round(xP,2),",",round(yP,2),")",sep=""),ylim=c(-1,1))
points(TimePoints,DataMatrix[i,],pch=16,col="red")
points(seq(min(TimePoints),max(TimePoints),length.out=100),fs.test(xP,yP)*cos(seq(min(TimePoints),max(TimePoints),length.out=100)),type='l',col="blue",lwd=2)
legend("bottomleft",c("reale", "stimata"), lty = c(1,1),col=c("blue","black"),lwd=2)
dev.off()

i<-108
xP<-xknot[i]
yP<-yknot[i]
png(filename=paste("Plot per un punto fissato",fs.test(xP,yP),".png",sep=" "))
FixedPointPlot(xP,yP,SolutionObj,lwd=2,NameLocation = paste("(",round(xP,2),",",round(yP,2),")",sep=""),ylim=c(-2,2))
points(TimePoints,DataMatrix[i,],pch=16,col="red")
points(seq(min(TimePoints),max(TimePoints),length.out=100),fs.test(xP,yP)*cos(seq(min(TimePoints),max(TimePoints),length.out=100)),type='l',col="blue",lwd=2)
legend("bottomleft",c("reale", "stimata"), lty = c(1,1),col=c("blue","black"),lwd=2)
dev.off()

i<-38
xP<-xknot[i]
yP<-yknot[i]
png(filename=paste("Plot per un punto fissato",fs.test(xP,yP),".png",sep=" "))
FixedPointPlot(xP,yP,SolutionObj,lwd=2,NameLocation = paste("(",round(xP,2),",",round(yP,2),")",sep=""),ylim=c(-4.5,4.5))
points(TimePoints,DataMatrix[i,],pch=16,col="red")
points(seq(min(TimePoints),max(TimePoints),length.out=100),fs.test(xP,yP)*cos(seq(min(TimePoints),max(TimePoints),length.out=100)),type='l',col="blue",lwd=2)
legend("bottomleft",c("reale", "stimata"), lty = c(1,1),col=c("blue","black"),lwd=2)
dev.off()



##### PLOT DEI RESIDUI #####

#Ricavo il vettore con tutte le stime
zHat<-NULL
for(i in 1:length(xknot))
{
    zHat<-c(zHat,ST.Eval(rep(xknot[i],length(TimePoints)),rep(yknot[i],length(TimePoints)),TimePoints,SolutionObj))
}


Residuals=Data-zHat

# Faccio qualche plot
png(filename="Scatterplot1.png")
plot(Data,Residuals,main="Scatterplot Residui",xlab="Dati",ylab="Residui")
abline(h=0,col="blue",lwd=2)
dev.off()
# Faccio qualche plot
png(filename="Scatterplot2.png")
plot(zHat,Residuals,main="Scatterplot Residui",xlab="Valori stimati",ylab="Residui")
abline(h=0,col="blue",lwd=2)
dev.off()

png(filename="QQplot.png")
qqnorm(Residuals)
qqline(Residuals,lwd=2,col="blue")
dev.off()

sink(file="Shapiro Test.txt")
print(shapiro.test(Residuals))
sink()