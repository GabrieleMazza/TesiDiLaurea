# Plot di triangolazioni

# Inserisco il file da leggere
load("FrontieraC.RData")

# Inserisco i punti da leggere
x<-x
y<-y
Internal<-TypePoint

nint<-length(x[Internal])
png(filename = "Triangolazione.png")
plot(x,y,type="n",xlab=" ",ylab=" ",main=paste("Triangolazione (",nint,"/",length(x),").png",sep=""))
for (ne in 1:dim(Triang)[1])
{
    polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
}
dev.off()
