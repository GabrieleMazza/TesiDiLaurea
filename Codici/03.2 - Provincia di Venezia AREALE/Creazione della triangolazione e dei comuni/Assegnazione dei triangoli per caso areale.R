# Creazione della triangolazione e assegnamento ai comuni per il caso areale
 
library(RTriangle)
library(SDMTools)
library(geosphere)
load("TerritorioOld.RData")
# Qui leggo le informazioni usate dal caso con dati puntuali
# Quindi ci sono anche i punti dei dati replicati!

xbound=x[is.na(Codici)]
ybound=y[is.na(Codici)]
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

# SELEZIONA L'AREA MASSIMA
amax=0.001
mesh<-triangulate(pslg_obj,Y=FALSE,D=TRUE,a=amax)

Triang<-mesh$T
xnew<-mesh$P[,1]
ynew<-mesh$P[,2]
# Inserisco i punti da leggere
png(filename=paste("Triangolazione",amax,".png"))
plot(xnew,ynew,type="n",xlab=" ",ylab=" ",main=paste("Triangolazione (da ",length(x)," a ",length(xnew), " punti, Amax=",amax,", ", dim(Triang)[1], " triangoli)",sep=""))
for (ne in 1:dim(Triang)[1])
{
    polygon(c(xnew[Triang[ne,1]],xnew[Triang[ne,2]],xnew[Triang[ne,3]]),c(ynew[Triang[ne,1]],ynew[Triang[ne,2]],ynew[Triang[ne,3]]))
}
dev.off()






##### ASSEGNAZIONE DEI TRIANGOLI AI COMUNI #####
# Scarico la struttura dei comuni
require(raster)
ProvinciaVZ =  subset(getData('GADM', country='ITA', level=3), NAME_2=="Venezia")
#plot(ProvinciaVZ)

#save(file="ProvinciaVZ.RData",ProvinciaVZ)
load(file="ProvinciaVZ.RData")

Nomi<-slot(ProvinciaVZ, 'data')[,c(1,10)]

# La frontiera è memorizzata come una unione di 131 poligoni (a causa delle isole della)
# laguna veneta)

# Calcolo i baricentri
# Un triangolo sarà asseganto ad un comune se il baricentro è in quel territorio
xG<-NULL
yG<-NULL
for (i in 1:(dim(Triang)[1]))
{
    xG<-c(xG,(xnew[Triang[i,1]]+xnew[Triang[i,2]]+xnew[Triang[i,3]])/3)
    yG<-c(yG,(ynew[Triang[i,1]]+ynew[Triang[i,2]]+ynew[Triang[i,3]])/3)
}

##### CICLO DI ASSEGNAZIONE #####
tmp<- slot(ProvinciaVZ, 'polygons')
sub.tmp <- slot(tmp[[2]],'Polygons') 
order<-slot(sub.tmp[[17]],"coords")
labels<-numeric(length(xG))

# siccome alcuni comuni sono costituiti da più di un poligono (tipo venezia)
# devo controllare più poligoni per comune
for(i in 1:length(xG))
{
    xT<-xG[i]
    yT<-yG[i]
    
    done=FALSE
    for(j in 1:44)
    {
        tmp<- slot(ProvinciaVZ, 'polygons')
        sub.tmp <- slot(tmp[[j]],'Polygons') 
        for(k in 1:length(sub.tmp))
        {
            PolyPoints<-slot(sub.tmp[[k]],"coords")
            if(pnt.in.poly(cbind(xT,yT),PolyPoints)$pip==1)
            {
                if(done)
                {
                    print("DOPPIO COMUNE! ERRORE!")
                }
                labels[i]=j
                done=TRUE
            }
        }
        
    }
}

# ORA LABELS CONTIENE PER OGNI TRIANGOLO IL COMUNE IN CUI RICADE IL SUO BARICCENTRO
# ALCUNI PERO' NON SONO IDENTIFICATI (codificati con zero)

print(paste("Assegnati ",sum(labels!=0), " triangoli su ", length(xG)))
print(paste(sum(labels==0), "mancanti"))

Colors<-rainbow(44)
png(filename="Comuni da assegnare.png")
plot(xnew,ynew,type="n",xlab=" ",ylab=" ",main=paste("Triangolazione (da ",length(x)," a ",length(xnew), " punti, Amax=",amax,", ", dim(Triang)[1], " triangoli)",sep=""))
for (ne in 1:dim(Triang)[1])
{
    if(labels[ne]!=0)
    {
        polygon(c(xnew[Triang[ne,1]],xnew[Triang[ne,2]],xnew[Triang[ne,3]]),c(ynew[Triang[ne,1]],ynew[Triang[ne,2]],ynew[Triang[ne,3]]),col=Colors[labels[ne]])
    }
    if(labels[ne]==0)
    {
        polygon(c(xnew[Triang[ne,1]],xnew[Triang[ne,2]],xnew[Triang[ne,3]]),c(ynew[Triang[ne,1]],ynew[Triang[ne,2]],ynew[Triang[ne,3]]),col="white")
    }
}
dev.off()

# Assegno manualmente i mancanti.. No troppi

# Decido di creare un metodo per definire come assegnare la distanza
# E questo metodo è: controllo la distanza del baricentro dal poligono
# e scelgo quello a distanza minore
# Funzione per la distanza punto-poligono
p2poly <- function(pt, poly){
    # Closing the polygon
    if(!identical(poly[1,],poly[nrow(poly),])){poly<-rbind(poly,poly[1,])}
    # A simple distance function
    dis <- function(x0,x1,y0,y1){sqrt((x0-x1)^2 +(y0-y1)^2)}
    d <- c()   # Your distance vector
    for(i in 1:(nrow(poly)-1)){
        ba <- c((pt[1]-poly[i,1]),(pt[2]-poly[i,2])) #Vector BA
        bc <- c((poly[i+1,1]-poly[i,1]),(poly[i+1,2]-poly[i,2])) #Vector BC
        dbc <- dis(poly[i+1,1],poly[i,1],poly[i+1,2],poly[i,2]) #Distance BC
        dp <- (ba[1]*bc[1]+ba[2]*bc[2])/dbc          #Projection of A on BC
        if(dp<=0){ #If projection is outside of BC on B side
            d[i] <- dis(pt[1],poly[i,1],pt[2],poly[i,2])
        }else if(dp>=dbc){ #If projection is outside of BC on C side
            d[i] <- dis(poly[i+1,1],pt[1],poly[i+1,2],pt[2])
        }else{ #If projection is inside of BC
            d[i] <- sqrt(abs((ba[1]^2 +ba[2]^2)-dp^2))
        }
    }
    min(d)
}

for(i in 1:length(xG))
{
    xT<-xG[i]
    yT<-yG[i]
    # Se non ho ancora processato il corrispondente triangolo...
    if(labels[i]==0)
    {
        d= Inf
        for(j in 1:44)
        {
            tmp<- slot(ProvinciaVZ, 'polygons')
            sub.tmp <- slot(tmp[[j]],'Polygons') 
            for(k in 1:length(sub.tmp))
            {
                # Se un comune ha più poligoni, scelgo quello a distanza minore
                PolyPoints<-slot(sub.tmp[[k]],"coords")
                tmp=p2poly(c(xT,yT),PolyPoints)
                if(tmp<d)
                {
                    labels[i]=j
                    d=tmp
                }
            }
            
        }
        
    }
}

print(paste("Assegnati ",sum(labels!=0), " triangoli su ", length(xG)))
print(paste(sum(labels==0), "mancanti"))

# Creo assegnamenti
Assegnamenti<-NULL
for(k in 1:(dim(Triang)[1]))
{
    Assegnamenti<-rbind(Assegnamenti,c(k,paste(Nomi[labels[k],2]),NA))
    Assegnamenti[k,3] <- which(Nomi$NAME_3==Assegnamenti[k,2])  
}

# Ora voglio fare una mappa colorata
Colors<-rainbow(44)
png(filename="Comuni.png")
plot(xnew,ynew,type="n",xlab=" ",ylab=" ",main=paste("Triangolazione (da ",length(x)," a ",length(xnew), " punti, Amax=",amax,", ", dim(Triang)[1], " triangoli)",sep=""))
for (ne in 1:dim(Triang)[1])
{
    polygon(c(xnew[Triang[ne,1]],xnew[Triang[ne,2]],xnew[Triang[ne,3]]),c(ynew[Triang[ne,1]],ynew[Triang[ne,2]],ynew[Triang[ne,3]]),col=Colors[labels[ne]])
}
dev.off()

## Ora plot di un comune alla volta...
for(i in 1:44)
{
    ID=which(Assegnamenti[,3]==i)
    if (length(ID)==0)
    {
        print(paste("COMUNE ",i," CON ZERO TRIANGOLI"))
    }
    # Ora ho i triangoli da stampare
    name=Nomi$NAME_3[i]
    png(filename=paste(name,".png"))
    plot(xnew,ynew,type="n",xlab=" ",ylab=" ",main=paste(name,", ",length(ID)," triangoli"))
    for (ne in 1:dim(Triang)[1])
    {
        polygon(c(xnew[Triang[ne,1]],xnew[Triang[ne,2]],xnew[Triang[ne,3]]),c(ynew[Triang[ne,1]],ynew[Triang[ne,2]],ynew[Triang[ne,3]]))
    }
    for (ne in ID)
    {
        polygon(c(xnew[Triang[ne,1]],xnew[Triang[ne,2]],xnew[Triang[ne,3]]),c(ynew[Triang[ne,1]],ynew[Triang[ne,2]],ynew[Triang[ne,3]]),col="red")
    }
    dev.off()
    
}



##### ERRORI DA CORREGGERE #####

# Correzione di labels per l'isola di pellestrina (sensibile a chioggia ma
# in realtà è dentro venezia)
# Definisco un rettangolo
xmin=12.25
xmax=12.35
ymin=45.235
ymax=45.4
for(i in 1:length(xG))
{
    xT<-xG[i]
    yT<-yG[i]
    
    if(xT>xmin & xT<xmax)
    {
        if(yT>ymin & yT<ymax)
        {
            labels[i]=43
        }
            
    }    
}


# CONTROLLO
# 9 chioggia
# 43 venezia
# TOGLIERE I SEGUENTI COMMENTI PER FARE IL PLOT INGRANDITO

# i=43
# ID=which(Assegnamenti[,3]==i)
# if (length(ID)==0)
# {
#     print(paste("COMUNE ",i," CON ZERO TRIANGOLI"))
# }
# # Ora ho i triangoli da stampare
# name=Nomi$NAME_3[i]
# png(filename=paste(name,".png"))
# plot(xnew,ynew,type="n",xlab=" ",ylab=" ",main=paste(name,", ",length(ID)," triangoli"),xlim=c(12.25,12.4),ylim=c(45.225,45.28))
# for (ne in 1:dim(Triang)[1])
# {
#     polygon(c(xnew[Triang[ne,1]],xnew[Triang[ne,2]],xnew[Triang[ne,3]]),c(ynew[Triang[ne,1]],ynew[Triang[ne,2]],ynew[Triang[ne,3]]))
# }
# for (ne in ID)
# {
#     polygon(c(xnew[Triang[ne,1]],xnew[Triang[ne,2]],xnew[Triang[ne,3]]),c(ynew[Triang[ne,1]],ynew[Triang[ne,2]],ynew[Triang[ne,3]]),col="red")
# }
# dev.off()


# ERRORI NELLA PARTE SUPERIORE DELLA LAGUNA
plot(xnew,ynew,type="n",xlab=" ",ylab=" ",xlim=c(12.5,12.62),ylim=c(45.5,45.7))
for (ne in 1:dim(Triang)[1])
{
    if(labels[ne]!=0)
    {
        polygon(c(xnew[Triang[ne,1]],xnew[Triang[ne,2]],xnew[Triang[ne,3]]),c(ynew[Triang[ne,1]],ynew[Triang[ne,2]],ynew[Triang[ne,3]]),col=Colors[labels[ne]])
    }
    if(labels[ne]==0)
    {
        polygon(c(xnew[Triang[ne,1]],xnew[Triang[ne,2]],xnew[Triang[ne,3]]),c(ynew[Triang[ne,1]],ynew[Triang[ne,2]],ynew[Triang[ne,3]]),col="white")
    }
}
# Qui è possibile che siano dati dei comuni a venezia poichè vicini alle isole interne
# della laguna
# Venezia è rossa -> CODICE 43
# San donà di piave viola (ha uno sbocco sul mare in teoria) -> 34
# musile di piave azzurro -> 26

# Quindi li correggo identificando i baricentri
points(xG,yG,pch=16,col="white")
identify(xG,yG)

# Modifico manualmente quelli che devono essere cambiati
labels[537]=34


plot(xnew,ynew,type="n",xlab=" ",ylab=" ",xlim=c(12.5,12.62),ylim=c(45.5,45.7))
for (ne in 1:dim(Triang)[1])
{
    if(labels[ne]!=0)
    {
        polygon(c(xnew[Triang[ne,1]],xnew[Triang[ne,2]],xnew[Triang[ne,3]]),c(ynew[Triang[ne,1]],ynew[Triang[ne,2]],ynew[Triang[ne,3]]),col=Colors[labels[ne]])
    }
    if(labels[ne]==0)
    {
        polygon(c(xnew[Triang[ne,1]],xnew[Triang[ne,2]],xnew[Triang[ne,3]]),c(ynew[Triang[ne,1]],ynew[Triang[ne,2]],ynew[Triang[ne,3]]),col="white")
    }
}



##### CONCLUSIONE #####

Assegnamenti<-NULL
for(k in 1:(dim(Triang)[1]))
{
    Assegnamenti<-rbind(Assegnamenti,c(k,paste(Nomi[labels[k],2]),NA))
    Assegnamenti[k,3] <- which(Nomi$NAME_3==Assegnamenti[k,2])  
}

## Ora plot di un comune alla volta...
for(i in 1:44)
{
    ID=which(Assegnamenti[,3]==i)
    if (length(ID)==0)
    {
        print(paste("COMUNE ",i," CON ZERO TRIANGOLI"))
    }
    # Ora ho i triangoli da stampare
    name=Nomi$NAME_3[i]
    png(filename=paste(name,".png"))
    plot(xnew,ynew,type="n",xlab=" ",ylab=" ",main=paste(name,", ",length(ID)," triangoli"))
    for (ne in 1:dim(Triang)[1])
    {
        polygon(c(xnew[Triang[ne,1]],xnew[Triang[ne,2]],xnew[Triang[ne,3]]),c(ynew[Triang[ne,1]],ynew[Triang[ne,2]],ynew[Triang[ne,3]]))
    }
    for (ne in ID)
    {
        polygon(c(xnew[Triang[ne,1]],xnew[Triang[ne,2]],xnew[Triang[ne,3]]),c(ynew[Triang[ne,1]],ynew[Triang[ne,2]],ynew[Triang[ne,3]]),col="red")
    }
    dev.off()
    
}


# Salvo i file degli assegnamenti
write.table(Assegnamenti,file="Assegnamenti.txt",row.names=FALSE,col.names=FALSE)

# I punti
Points<-data.frame(xTriang=xnew,yTriang=ynew)
write.table(Points,file="TriangulationPoints.txt",row.names=FALSE,col.names=FALSE)

# La triangolazione
write.table(Triang,file="Triangulation.txt",row.names=FALSE,col.names=FALSE)
