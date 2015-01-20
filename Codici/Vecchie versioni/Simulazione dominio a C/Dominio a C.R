# Simulate spatial data on a C-shaped domain
# Last modified 24 July 2010 by laura

library(rgl)
library(mgcv)
library(fda)
library(RTriangle)
library(SDMTools)
source("SpazioTempo.R")

######################################################################
######################   DATA GENERATION   ###########################
######################################################################


# Plot the function, and its boundary
fsb <- list(fs.boundary())
nx<-30
ny<-10
#Sequenza delle x e delle y
xvec <- seq(-1,4,length=nx)
yvec<-seq(-1,1,length=ny)

#Li replico e formo una specie di dominio quadrato
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))
#Replica ogni valore di yvec 250 volte


iY=1:ny
#fstest fa tutto da solo, se nx e ny stanno dentro il domino, gli da un valore
#se no na
tru <- matrix(fs.test(xx,yy),nx,ny) ## truth
#Corrispondenti ai punti di xvec e yvec


image(xvec,yvec,tru,col=heat.colors(100),xlab="x",ylab="y",asp=1)
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=3)
contour(xvec,yvec,tru,levels=seq(-5,5,by=.25),add=TRUE)


#Ora studio tru, creo il vettore delle x e delle y con i dati e con la soluzione

x<-NULL
y<-NULL
data<-NULL
for(i in 1:(dim(tru)[1]))
{
    for(j in 1:(dim(tru)[2]))
    {
        if(!is.na(tru[i,j]))
        {
            x<-c(x,xvec[i])
            y<-c(y,yvec[j])
            data<-c(data,tru[i,j])
        }
    }
}

xbound<-fsb[[1]]$x
ybound<-fsb[[1]]$y

#Ci sono punti esterni alla frontiera?
PolyPoints<-cbind(xbound,ybound)
if(sum(pnt.in.poly(cbind(x,y),PolyPoints)$pip)==length(x))
{
    print("Tutti i punti stanno dentro")
} else
{
    n<-length(x)-sum(pnt.in.poly(cbind(x,y),PolyPoints)$pip)
    paste("Esistono",n,"punti esterni alla frontiera",sep="")
}

ID<-pnt.in.poly(cbind(x,y),PolyPoints)$pip==1
#Li tolgo
x<-x[ID]
y<-y[ID]
data<-data[ID]

#Ora riprovo
PolyPoints<-cbind(xbound,ybound)
if(sum(pnt.in.poly(cbind(x,y),PolyPoints)$pip)==length(x))
{
    print("Tutti i punti stanno dentro")
} else
{
    print("Esistono punti esterni alla frontiera")
}

#Creo la triangolazione spaziale
#Oggetti completi
xtot<-c(x,xbound)
ytot<-c(y,ybound)
#Creo i Boundaries
Boundaries<-NULL
for(i in (length(x)+1):(length(xtot)-1))
{
    Boundaries<-rbind(Boundaries, c(i,i+1))
}
Boundaries<-rbind(Boundaries, c(length(xtot),length(x)+1))
#Ora triangolazione
#Oggetto pslg
pslg_obj<-pslg(cbind(xtot,ytot),S=Boundaries)
#Creo la mesh
#Y dice di non aggiungere Steiner Points
#D dice di triangolare con Delaunay
mesh<-triangulate(pslg_obj,Y=TRUE,D=TRUE)
#Estrazione dei triangoli
Triang<-mesh$T
#Plot della triangolazione
plot(xtot,ytot,col="white")
for (ne in 1:dim(Triang)[1])
{
    polygon(c(xtot[Triang[ne,1]],xtot[Triang[ne,2]],xtot[Triang[ne,3]]),c(ytot[Triang[ne,1]],ytot[Triang[ne,2]],ytot[Triang[ne,3]]))
}


### PERTURBAZIONE TEMPORALE DEI DATI ###
#Perturbo secondo la funzione cos(x)
#Ora introduco una variazione temporale. Creo 5 istanti di tempo
TimePoints<-0:5
Data<-NULL
for(i in TimePoints)
{
    Data<-cbind(Data,data*cos(i))
}

#Creo le basi in spazio e tempo
TimeBasisObj<-Create.Bspline.Time.Basis(TimePoints,3,T)
SpaceBasisObj<-Create.FEM.Space.Basis(cbind(xtot,ytot),Triang,1)

