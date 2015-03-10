# Aumento le triangolazioni
# Voglio rendere il numero di triangoli più alto, diminuendo le aree

library(RTriangle)
load("TerritorioOld.RData")

# Creo l'oggetto boundaries
nint<-length(Codici[!is.na(Codici)])
Boundaries<-NULL
for(i in (nint+1):(length(x)-1))
{
    Boundaries<-rbind(Boundaries, c(i,i+1))
}
Boundaries<-rbind(Boundaries, c(length(x),nint+1))

# Ora pslg object
pslg_obj<-pslg(cbind(x,y),S=Boundaries)

#Creo la mesh
#Y dice di non aggiungere Steiner Points
#D dice di triangolare con Delaunay
mesh<-triangulate(pslg_obj,Y=FALSE,D=TRUE,a=0.001)

Triang<-mesh$T
xnew<-mesh$P[,1]
ynew<-mesh$P[,2]
# Inserisco i punti da leggere
png(filename="Triangolazione.png")
plot(xnew,ynew,type="n",xlab=" ",ylab=" ",main=paste("Triangolazione (da ",length(x)," a ",length(xnew), " punti)",sep=""))
for (ne in 1:dim(Triang)[1])
{
    polygon(c(xnew[Triang[ne,1]],xnew[Triang[ne,2]],xnew[Triang[ne,3]]),c(ynew[Triang[ne,1]],ynew[Triang[ne,2]],ynew[Triang[ne,3]]))
}
dev.off()

#Controllo
sum(abs(x-xnew[1:length(x)]))

Codici<-c(Codici,rep(NA,(length(xnew)-length(x))))
x<-xnew
y<-ynew

save(file="TerritorioNew.RData",x,y,Codici,Triang,xbound,ybound)
