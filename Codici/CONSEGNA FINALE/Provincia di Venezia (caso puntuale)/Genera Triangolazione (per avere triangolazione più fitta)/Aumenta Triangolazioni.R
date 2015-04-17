# This script is used to create the riangulation

# INPUT FILES
#       Boundary.txt -> Boundary of Venice Province
#       CoordinateCovariate.txt -> Points of municipalities

# OUTPUT FILES
#       TriangulationPoints.txt -> Spatial points of triangulation
#       Triangulation.txt -> Index of TriangulationPoints for each triangle

library(RTriangle)

Boundary<-read.table(file="Boundary.txt",header=T)
CoordinateCovariate<-read.table(file="CoordinateCovariate.txt",header=T)

#Need total vectors of x and y, for triangulation
x<-c(CoordinateCovariate$Longitudine,Boundary$xBound)
y<-c(CoordinateCovariate$Latitudine,Boundary$yBound)
#Boundaries object for RTriangle
Boundaries<-NULL
for(i in (length(CoordinateCovariate$Longitudine)+1):(length(x)-1))
{
    Boundaries<-rbind(Boundaries, c(i,i+1))
}
Boundaries<-rbind(Boundaries, c(length(x),length(CoordinateCovariate$Longitudine)+1))

# pslg object
pslg_obj<-pslg(cbind(x,y),S=Boundaries)
#Maximum value for area
amax=0.003
mesh<-triangulate(pslg_obj,Y=FALSE,D=TRUE,a=amax)
Triang<-mesh$T
dimnames(Triang)[[2]]<-c("Vertex1","Vertex2","Vertex3")
xnew<-mesh$P[,1]
ynew<-mesh$P[,2]
Points<-data.frame(xTriang=xnew,yTriang=ynew)

# Plot
png(filename="Triangulation.png")
plot(xnew,ynew,type="n",xlab=" ",ylab=" ",main=paste("Triangulation (",length(xnew), " points, ",dim(Triang)[1]," triangles)",sep=""))
for (ne in 1:dim(Triang)[1])
{
    polygon(c(xnew[Triang[ne,1]],xnew[Triang[ne,2]],xnew[Triang[ne,3]]),c(ynew[Triang[ne,1]],ynew[Triang[ne,2]],ynew[Triang[ne,3]]))
}
dev.off()

#Save new files
write.table(file="Triangulation.txt",Triang,row.names=F)
write.table(file="TriangulationPoints.txt",Points,row.names=F)