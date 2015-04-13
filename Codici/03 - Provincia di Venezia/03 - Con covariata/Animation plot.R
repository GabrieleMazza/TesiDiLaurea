library(rgl)
library(mgcv)
library(SDMTools)
library(ggplot2)
source("SpazioTempo.R")
load("Territorio.RData")
     
library(animation)
ani.options(convert = 'C:\\Program Files\\ImageMagick-6.9.0-Q16\\convert.exe')

set.seed(10)
# Devo aggiungere del rumore al dato?
noise<-TRUE

SigmaNoise<-0.5

# Leggo i dati
CoordinateCovariate<-read.table("CoordinateCovariate.txt",header=T)
Risposta<-read.table("Risposta.txt",header=T)

# Punti di tempo
TimePoints<-1997:2011
InternalPoints<-!is.na(Codici)
nint<-sum(InternalPoints)

# Basi
TimeBasisObj<-Create.Bspline.Time.Basis(TimePoints,TimeOrder=4,DerivativeOrder=2,PlotIt=F)
SpaceBasisObj<-Create.FEM.Space.Basis(cbind(x,y),Triang,InternalPoints,1)

SolutionObj=ReadSolutionObjCovar(SpaceBasisObj,TimeBasisObj,"VettoreC.txt","BetaHat.txt")


nx<-200
ny<-200
xvec<-seq(min(x),max(x),length.out=nx)
yvec<-seq(min(y),max(y),length.out=ny)
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))




TimeValid=seq(1997,2011,length.out=100)

ResultFittedAni<-NULL
for(j in TimeValid)
{
    Time<-rep(j,length(xx))
    Result<-ST.Eval(xx,yy,Time,SolutionObj)
    
    ResultFittedAni<-cbind(ResultFittedAni,Result)
}

save(file="Animation.RData", ResultFittedAni,TimeValid)

load(file="Animation.RData")
Max=max(ResultFittedAni,na.rm=TRUE)
Min=min(ResultFittedAni,na.rm=TRUE)

zlim<-c(Min,Max)

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
        layout(matrix(c(1,
                        2
        ), nrow=2,ncol=1,byrow = TRUE),height=c(1,4),width=c(1))
        
        # Set up the top chart that keeps track of the current frame/iteration
        # Dress it up a little just for fun
        plot(-5, xlim = c(1997,2011), ylim = c(0, .3), xlab = "", ylab = "", main = "Anno",axes=F)
        abline(v=TimeValid[j], lwd=5, col = rgb(0, 0, 255, 255, maxColorValue=255))
        
        # Bring back the X axis
        xticks <- seq(1997,2011)
        axis(side=1, at=xticks, labels=xticks)
        
        
        # Matrice
        Mat <- matrix(ResultFittedAni[,j],nrow=nx,ncol=ny,byrow=F)
        Mat<-cbind(Mat,numeric(dim(Mat)[1]))
        Mat[1,ny+1]=zlim[1]
        Mat[2,ny+1]=zlim[2]
        
        # Plot
        image(xvec,c(yvec,50),Mat,col=heat.colors(100),main="Funzione stimata",xlab="",ylab="",ylim=c(yvec[1],yvec[ny]))
        lines(xbound,ybound,lwd=1)
        lines(c(xbound[1],xbound[length(xbound)]),c(ybound[1],ybound[length(ybound)]),lwd=1)
        contour(xvec,c(yvec,50),Mat,nlevels=10,add=TRUE)
           }
})


#STAMPA PER ANIMAZIONE A SCORRIMENTO
for (j in 1:length(TimeValid))
{
    png(file=paste("image",j,".png",sep=""))
    par(oma=c(4,0,4,0),mar=c(2,2,2,2))
    layout(matrix(c(1,
                    2
    ), nrow=2,ncol=1,byrow = TRUE),height=c(1,4),width=c(1))
    
    # Set up the top chart that keeps track of the current frame/iteration
    # Dress it up a little just for fun
    plot(-5, xlim = c(1997,2011), ylim = c(0, .3), xlab = "", ylab = "", main = "Anno",axes=F)
    abline(v=TimeValid[j], lwd=5, col = rgb(0, 0, 255, 255, maxColorValue=255))
    
    # Bring back the X axis
    xticks <- seq(1997,2011)
    axis(side=1, at=xticks, labels=xticks)
    
    
    # Matrice
    Mat <- matrix(ResultFittedAni[,j],nrow=nx,ncol=ny,byrow=F)
    Mat<-cbind(Mat,numeric(dim(Mat)[1]))
    Mat[1,ny+1]=zlim[1]
    Mat[2,ny+1]=zlim[2]
    
    # Plot
    image(xvec,c(yvec,50),Mat,col=heat.colors(100),main="Funzione stimata",xlab="",ylab="",ylim=c(yvec[1],yvec[ny]))
    lines(xbound,ybound,lwd=1)
    lines(c(xbound[1],xbound[length(xbound)]),c(ybound[1],ybound[length(ybound)]),lwd=1)
    contour(xvec,c(yvec,50),Mat,nlevels=10,add=TRUE)
    dev.off()
}