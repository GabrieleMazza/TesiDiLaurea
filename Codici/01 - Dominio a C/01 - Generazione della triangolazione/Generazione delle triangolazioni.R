# CREAZIONE DEL DOMINIO A C E DELLA TRIANGOLAZIONE ASSOCIATA

library(rgl)
library(mgcv)
library(SDMTools)
library(RTriangle)
source("Functions.R")

# Ho già di base 108 punti di frontiera..
# Quanti punti simulare?
n=1600

# Inanzitutto considero la funzione reale
fsb <- list(fs.boundary())
nx<-250
ny<-100 
xvec <- seq(-1,3.5,length=nx)
yvec<-seq(-1,1,length=ny)
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))
iY=1:ny
tru <- matrix(fs.test(xx,yy),nx,ny) ## truth
png(filename="Funzione reale.png")
image(xvec,yvec,tru,col=heat.colors(100),xlab="x",ylab="y",asp=1,ylim=c(-1,1))
lines(fsb[[1]]$x,fsb[[1]]$y,lwd=3)
contour(xvec,yvec,tru,levels=seq(-5,5,by=.25),add=TRUE)
dev.off()

# Carico la frontiera da file
 
Bound<-read.table(file = "Frontiera.txt", header=T) 

attach(Bound)

# Ci sono 108 punti di frontiera secondo questa definizione
# Ora simulo i dati interni

x <- runif(n)*5-1
y <-runif(n)*2-1
muvec <- fs.test(x,y,b=1)
ind <- inSide(fsb,x=x,y=y) ## remove outsiders
muvec <- muvec[ind]
x <- x[ind]
y <- y[ind] 

# Quanti punti sono rimasti?
nint<-length(x)
nint
nbound<-length(xbound)

# Creo gli oggetti totali
x<-c(x,xbound)
y<-c(y,ybound)

Boundaries<-NULL
for(i in (nint+1):(length(x)-1))
{
    Boundaries<-rbind(Boundaries, c(i,i+1))
}
Boundaries<-rbind(Boundaries, c(length(x),nint+1))
#Ora triangolazione
#Oggetto pslg
pslg_obj<-pslg(cbind(x,y),S=Boundaries)
#Creo la mesh
#Y dice di non aggiungere Steiner Points
#D dice di triangolare con Delaunay
mesh<-triangulate(pslg_obj,Y=TRUE,D=TRUE)
#Estrazione dei triangoli
Triang<-mesh$T
#Plot della triangolazione
png(filename = paste("Triangolazione.png",sep=""))
plot(x,y,type="n")
for (ne in 1:dim(Triang)[1])
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
}
dev.off()
#Controllo se la triangolazione ha triangoli con solo punti di bordo
BorderTR<-BorderTriangles(mesh$T,Boundaries)
BorderTR
#Li coloro
for (ne in BorderTR)
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]),col="green")
}
detach(Bound)

# Ora mi occupo di generare i dati per i punti interni al dominio
Data <- fs.test(x[1:nint],y[1:nint])
Data<-as.vector(Data)

# Creo un vettore che mi dica se il punto è interno o di frontiera
TypePoint<-rep(FALSE,length(x))
TypePoint[1:nint]=TRUE

save(file="FrontieraC.RData",x,y,TypePoint,Data,Triang)



### SALVATAGGIO DEI PUNTI PER LA VALIDAZIONE DEL RISULTATO ###

# Ho a disposizione xvec, yvec e tru
xValid<-NULL
yValid<-NULL
PosMatrix<-NULL

for(i in 1:(dim(tru)[1]))
{
    for(j in 1:(dim(tru)[2]))
    {
        if(!is.na(tru[i,j]))
        {
            xValid<-c(xValid,xvec[i])
            yValid<-c(yValid,yvec[j])
            PosMatrix<-rbind(PosMatrix,c(i,j))
        }
    }
}

save(file="Validazione.RData",xvec,yvec,xValid,yValid,PosMatrix)
