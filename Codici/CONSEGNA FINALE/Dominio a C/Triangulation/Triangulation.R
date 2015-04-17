# Creation of boundary and triangulation for C shaped domain
# INPUT:
#   DomC_Boundary.txt -> Boundary definition
# OUTPUT:
#   DomC_InternalPoints -> Internal Points
#   DomC_Triangulation -> Triangulation
#   DomC_TriangulationPoints -> Triangulation Points

library(rgl)
library(mgcv)
library(SDMTools)
library(RTriangle)
source("Functions.R")

# Read boundary
Bound=read.table(file="DomC_Boundary.txt",header=T)

# There are 108 points for boundary. How many points for generation?
# Only those that will be in C shaped domai will be valid
set.seed(1)
n=150
# Bigger n means more triangles

fsb <- list(fs.boundary())
# Simulation of internal points
x <- runif(n)*5-1
y <-runif(n)*2-1
muvec <- fs.test(x,y,b=1)
ind <- inSide(fsb,x=x,y=y) ## remove outsiders
muvec <- muvec[ind]
x <- x[ind]
y <- y[ind] 

# How many points?
nint<-length(x)
nint

Points<-data.frame(xInt=x,yInt=y)
write.table(Points,file="DomC_InternalPoints.txt",row.names=F)


# Objects for triangulation
xtot<-c(x,Bound$xBound)
ytot<-c(y,Bound$yBound)

Boundaries<-NULL
for(i in (nint+1):(length(xtot)-1))
{
    Boundaries<-rbind(Boundaries, c(i,i+1))
}
Boundaries<-rbind(Boundaries, c(length(xtot),nint+1))

# Now I triangulate
pslg_obj<-pslg(cbind(xtot,ytot),S=Boundaries)
mesh<-triangulate(pslg_obj,Y=TRUE,D=TRUE)
Triang<-mesh$T

TPoints=data.frame(xTriang=xtot,yTriang=ytot)
write.table(Triang,file="DomC_Triangulation.txt",row.names=F,col.names=F)
write.table(TPoints,file="DomC_TriangulationPoints.txt",row.names=F)

png(filename = paste("Triangulation.png",sep=""))
plot(xtot,ytot,type="n")
for (ne in 1:dim(Triang)[1])
{
    polygon(c(xtot[Triang[ne,1]],xtot[Triang[ne,2]],xtot[Triang[ne,3]]),c(ytot[Triang[ne,1]],ytot[Triang[ne,2]],ytot[Triang[ne,3]]))
}
points(Bound$xBound,Bound$yBound,type='l')
dev.off()

