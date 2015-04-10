library(rgl)
library(mgcv)
library(SDMTools)
library(ggplot2)
source("SpazioTempo.R")
load("FrontieraC.RData")
load("Validazione.RData")
     
library(animation)
ani.options(convert = 'C:\\Program Files\\ImageMagick-6.9.0-Q16\\convert.exe')

set.seed(10)
# Devo aggiungere del rumore al dato?
noise<-TRUE

# DEVIAZIONE STANDARD DEL RUMORE
SigmaNoise=0.5

# Per la covariata
Mu<-0
Sigma<-1
Beta<-1

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

SolutionObj=ReadSolutionObjCovar(SpaceBasisObj,TimeBasisObj,"VettoreC.txt","BetaHat.txt")


TimeValid=seq(0,2*pi,length.out=100)

ResultFittedAni<-NULL
for(j in TimeValid)
{
    Time<-rep(j,length(xValid))
    Result<-ST.Eval(xValid,yValid,Time,SolutionObj)
    
    ResultFittedAni<-cbind(ResultFittedAni,Result)
}

ResultRealAni<-NULL

for(j in TimeValid)
{
    Time<-rep(j,length(xValid))
    Result<-fun(xValid,yValid,Time)
    ResultRealAni<-cbind(ResultRealAni,Result)
}

save(file="Animation.RData", ResultFittedAni,ResultRealAni,TimeValid)

load(file="Animation.RData")
Max=max(ResultFittedAni,ResultRealAni,na.rm=TRUE)
Min=min(ResultFittedAni,ResultRealAni,na.rm=TRUE)
zlim<-c(Min,Max)

#AGGIUSTO LA COLONNA IN PIU'... QUI AGGIUNGO IL MASSIMO E IL MINIMO
PlotMatrix=matrix(ncol=(length(yvec)+1),nrow=length(xvec))
PlotMatrix[1,(length(yvec)+1)]=Max
PlotMatrix[2,(length(yvec)+1)]=Min

# GRAFICI

# Ora salvo tutti i grafici
# Voglio ricreare perfettamente le posizioni della image matrix tru

nx<-length(xvec)
ny<-length(yvec)
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))


ani.options(interval=.5)
# Begin animation loop
# Note the brackets within the parentheses
saveGIF({
    
    # For the most part, it's safest to start with graphical settings in 
    # the animation loop, as the loop adds a layer of complexity to 
    # manipulating the graphs. For example, the layout specification needs to 
    # be within animation loop to work properly.
    layout(matrix(c(1, rep(2, 5)), 6, 1))
    
    # Adjust the margins a little
    par(mar=c(4,4,2,1) + 0.1)
    
    # Begin the loop that creates the 150 individual graphs
    for (j in 1:length(TimeValid))
    {
        par(oma=c(4,0,4,0),mar=c(2,2,2,2))
        layout(matrix(c(1,1,
                        2,3
        ), nrow=2,ncol=2,byrow = TRUE),height=c(1,4),width=c(1,1))
        
        # Set up the top chart that keeps track of the current frame/iteration
        # Dress it up a little just for fun
        plot(-5, xlim = c(0,2*pi), ylim = c(0, .3), xlab = "", ylab = "", main = "Time",axes=F)
        abline(v=TimeValid[j], lwd=5, col = rgb(0, 0, 255, 255, maxColorValue=255))
        
        # Bring back the X axis
        xticks <- round(seq(0,2*pi,length.out=9),2)
        axis(side=1, at=xticks, labels=xticks)
        
        
        # Funzione Reale
        for(i in 1:(dim(PosMatrix)[1]))
        {
            PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=ResultRealAni[i,j]
        }
        # Plot
        image(xvec,c(yvec,4),PlotMatrix,col=heat.colors(100),ylim=c(-1,+1),xlab=" ",ylab=" ",main="Funzione reale")
        lines(xbound,ybound,lwd=3)
        lines(c(xbound[1],xbound[length(xbound)]),c(ybound[1],ybound[length(ybound)]),lwd=3)
        contour(xvec,c(yvec,4),PlotMatrix,nlevels=10,add=TRUE,lwd=2,labcex=1.1)
        
        # Matrice con la stimata
        for(i in 1:(dim(PosMatrix)[1]))
        {
            PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=ResultFittedAni[i,j]
        }
        # Plot
        image(xvec,c(yvec,4),PlotMatrix,col=heat.colors(100),ylim=c(-1,+1),xlab=" ",ylab=" ",main="Funzione stimata")
        lines(xbound,ybound,lwd=3)
        lines(c(xbound[1],xbound[length(xbound)]),c(ybound[1],ybound[length(ybound)]),lwd=3)
        contour(xvec,c(yvec,4),PlotMatrix,nlevels=10,add=TRUE,lwd=2,labcex=1.1)
    }
})


#STAMPA PER ANIMAZIONE A SCORRIMENTO
for (j in 1:length(TimeValid))
{
    png(file=paste("image",j,".png",sep=""))
    par(oma=c(4,0,4,0),mar=c(2,2,2,2))
    layout(matrix(c(1,1,
                    2,3
    ), nrow=2,ncol=2,byrow = TRUE),height=c(1,4),width=c(1,1))
    
    # Set up the top chart that keeps track of the current frame/iteration
    # Dress it up a little just for fun
    plot(-5, xlim = c(0,2*pi), ylim = c(0, .3), xlab = "", ylab = "", main = "Time",axes=F)
    abline(v=TimeValid[j], lwd=5, col = rgb(0, 0, 255, 255, maxColorValue=255))
    
    # Bring back the X axis
    xticks <- round(seq(0,2*pi,length.out=9),2)
    axis(side=1, at=xticks, labels=xticks)
    
    
    # Funzione Reale
    for(i in 1:(dim(PosMatrix)[1]))
    {
        PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=ResultRealAni[i,j]
    }
    # Plot
    image(xvec,c(yvec,4),PlotMatrix,col=heat.colors(100),ylim=c(-1,+1),xlab=" ",ylab=" ",main="Funzione reale")
    lines(xbound,ybound,lwd=3)
    lines(c(xbound[1],xbound[length(xbound)]),c(ybound[1],ybound[length(ybound)]),lwd=3)
    contour(xvec,c(yvec,4),PlotMatrix,nlevels=10,add=TRUE,lwd=2,labcex=1.1)
    
    # Matrice con la stimata
    for(i in 1:(dim(PosMatrix)[1]))
    {
        PlotMatrix[PosMatrix[i,1],PosMatrix[i,2]]=ResultFittedAni[i,j]
    }
    # Plot
    image(xvec,c(yvec,4),PlotMatrix,col=heat.colors(100),ylim=c(-1,+1),xlab=" ",ylab=" ",main="Funzione stimata")
    lines(xbound,ybound,lwd=3)
    lines(c(xbound[1],xbound[length(xbound)]),c(ybound[1],ybound[length(ybound)]),lwd=3)
    contour(xvec,c(yvec,4),PlotMatrix,nlevels=10,add=TRUE,lwd=2,labcex=1.1)
    dev.off()
}