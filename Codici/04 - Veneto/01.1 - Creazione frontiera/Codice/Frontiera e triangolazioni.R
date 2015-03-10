# Creazione della frontiera in poligono unico usamdo una tecnica di smoothing
# La tecnica scelta è regression splines


# Funzioni di appoggio
source("Functions.R")

library(fda)
library(KernSmooth)
library(geosphere)
library(fda)
library(RTriangle)
library(SDMTools)




##### DOWNLOAD DELLA FRONTIERA #####

# Scarico la frontiera
# require(raster)
# veneto =  subset(getData('GADM', country='ITA', level=1), NAME_1=="Veneto")
# plot(veneto)
# save(file="Veneto.RData",veneto)
load("Veneto.RData")

# La frontiera è memorizzata come una unione di 131 poligoni (a causa delle isole della)
# laguna veneta)

tmp<- slot(veneto, 'polygons')
sub.tmp <- slot(tmp[[1]],'Polygons') 
order<-slot(tmp[[1]],'plotOrder')

xbound<-NULL
ybound<-NULL
labels<-NULL
for(i in 1:131)
{
    xbound<-c(xbound,sub.tmp[[order[i]]]@coords[,1])
    ybound<-c(ybound,sub.tmp[[order[i]]]@coords[,2])
    labels<-c(labels,rep(i,length(sub.tmp[[order[i]]]@coords[,1])))
}

# In totale ci sono 34690 punti per la provincia di venezia, in 131 poligoni

# Leggo anche i comuni
Comuni<-read.table(file="Coordinate.txt",header=T)





##### VENETO #####

# Ora devo studiare il caso del Veneto
# Il Veneto è nel poligono 1
# Studio tutti i dati in dipendenza dall'ascissa curvilinea


xVeneto_tot=xbound[labels==1]
yVeneto_tot=ybound[labels==1]
s<-0
previous<-0

for(i in 2:length(xVeneto_tot))
{
    # Aggiorno la distanza
    previous<-previous+sqrt((xVeneto_tot[i]-xVeneto_tot[i-1])^2+(yVeneto_tot[i]-yVeneto_tot[i-1])^2)
    s<-c(s,previous)  
}

# Ora devo fare smoothing di queste due funzioni

# Quanti punti fisso per la definizione del confine?
N<-150
# Nel range di ascissa curvilinea fisso N punti
Val<-seq(0,max(s),by=max(s)/N)
# Decido di usare delle splines cubiche
m<-4




##### REGRESSION SPLINES VENETO #####

do_cycles=FALSE




##### ANALISI DEL MIGLIOR nbasis #####

# Scelgo come upper bound il numero di punti della valutazione per il numero di basi con
# cui comporre la frontiera
# Come scelgo l'ottimo del numero di basi?
# Mi interessa che sia minimizzata la distanza totale quadratica euclidea
# poligono reale del Veneto
if(do_cycles==TRUE)
{
    bestsum<-9999999999999
    best<-0
    sum_vect<-NULL
    best_vect<-NULL
    #Qui mi salvo un po' di risultati numerici
    sink(file = "ROutput.txt", append = FALSE)
    print("nbasis sum best Intersections Trinagoli Comuni")
    for(nbasis in (m+1):N)
    {
        basis = create.bspline.basis(c(0,max(s)), nbasis, m)
        xsmooth = smooth.basis(argvals=s, y=xVeneto_tot, basis)
        ysmooth = smooth.basis(argvals=s, y=yVeneto_tot, basis)
        xVeneto  = eval.fd(Val, xsmooth$fd)
        yVeneto  = eval.fd(Val, ysmooth$fd)
        png(filename = paste(nbasis,".png",sep=""))
        plot(xVeneto_tot,yVeneto_tot,type='l',main=paste(nbasis,sep=""))
        points(xVeneto,yVeneto,col="blue",type='l')
        dev.off()
        sum=0
        for(i in 1:N)
        {
            sum<-sum+dist2Line(c(xVeneto[i],yVeneto[i]), cbind(xVeneto_tot,yVeneto_tot), distfun=SquareEuclideanDistance)[1]
        }
        if(sum<bestsum)
        {
            best<-nbasis
            bestsum=sum
        }
        Intersect<-Intersections(xVeneto,yVeneto)
        
        # Oggetti completi
        x<-c(Comuni$Longitudine,xVeneto)
        y<-c(Comuni$Latitudine,yVeneto)
        # Creo i Boundaries
        Boundaries<-NULL
        for(i in (length(Comuni$Longitudine)+1):(length(x)-1))
        {
            Boundaries<-rbind(Boundaries, c(i,i+1))
        }
        Boundaries<-rbind(Boundaries, c(length(x),length(Comuni$Longitudine)+1))
        # Ora triangolazione
        # Oggetto pslg
        pslg_obj<-pslg(cbind(x,y),S=Boundaries)
        # Creo la mesh
        # Y dice di non aggiungere Steiner Points
        # D dice di triangolare con Delaunay
        mesh<-triangulate(pslg_obj,Y=TRUE,D=TRUE)
        # Estrazione dei triangoli
        Triang<-mesh$T
        # Plot della triangolazione
        png(filename = paste("Triangolazione nella laguna ", nbasis, ".png", sep=""))
        plot(x,y,col="white")
        for (ne in 1:dim(Triang)[1])
        {
            polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
        }
        # Controllo se la triangolazione ha triangoli con solo punti di bordo
        BorderTR<-BorderTriangles(mesh$T,Boundaries)
        BorderTR
        # Li coloro
        for (ne in BorderTR)
        {
            polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]),col="green")
        }
        # Ci sono comuni esterni alla frontiera?
        PolyPoints<-cbind(xVeneto,yVeneto)
        points(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],Comuni$Latitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],col="red",pch=16)
        dev.off()
        
        length(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0)
        comuniout=580-length(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0])
        print(paste(
            nbasis,
            sum,
            best,
            dim(Intersect)[1],
            length(BorderTR),
            comuniout,
            sep=" "))
        sum_vect<-c(sum_vect,sum)
        best_vect<-c(best_vect,best)
    }
    sink()
    png(filename = "Sum.png")
    plot(sum_vect,type='l')
    dev.off()
    png(filename = "Best.png")
    plot(best_vect,type='l')
    dev.off()
}




##### RISULTATI E SMOOTHING FINALE #####

# Se io scelgo di avere 150 punti, il miglior nbasis è 136 con la distanza,
# ma a livello di quanti comuni restano fuori e di quanti triangoli vengono tolti, è meglio
# nbasis 145
nbasis<-145
basis = create.bspline.basis(c(0,max(s)), nbasis, m)
xsmooth = smooth.basis(argvals=s, y=xVeneto_tot, basis)
ysmooth = smooth.basis(argvals=s, y=yVeneto_tot, basis)
xVeneto  = eval.fd(Val, xsmooth$fd)
yVeneto  = eval.fd(Val, ysmooth$fd)

# Traccio ora i risultati
png(filename = "Veneto - Smoothing x (regression splines).png")
plot(s,xVeneto_tot,main="Smoothing in x",xlab="Ascissa curvilinea",ylab="xVeneto")
points(Val,xVeneto,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Veneto - Smoothing y (regression splines).png")
plot(s,yVeneto_tot,main="Smoothing in y",xlab="Ascissa curvilinea",ylab="yVeneto")
points(Val,yVeneto,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Veneto - Smoothing laguna (regression splines).png")
plot(xVeneto_tot,yVeneto_tot,type='l',main="Nuova definizione della regione",xlab="xVeneto",ylab="yVeneto",xlim=c(12.1,12.6),ylim=c(45,45.6))
points(xVeneto,yVeneto,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Veneto - Smoothing regione (regression splines).png")
plot(xVeneto_tot,yVeneto_tot,type='l',main="Nuova definizione della regione",xlab="xVeneto",ylab="yVeneto")
points(xVeneto,yVeneto,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()

# Controllo di intersezioni
Intersect<-Intersections(xVeneto,yVeneto)
Intersect

#Controllo se ci sono dei comuni che rimangono fuori dalla frontiera
PolyPoints<-cbind(xVeneto,yVeneto)
if(sum(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip)==length(Comuni$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}
#Quali comuni restano fuori?
Comuni$Comune[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0]




##### COMUNI RIMASTI FUORI #####

# Ingrandisco nelle zone interessate dai comuni rimasti fuori

# Inizio da Lestebasse
plot(xVeneto_tot,yVeneto_tot,col="black",xlim=c(11.2,11.3),ylim=c(45.85,45.95),type='l',main="Comune di Lestebasse",xlab="xVeneto",ylab="yVeneto")
points(xVeneto,yVeneto,type='l',col="blue")
points(xVeneto,yVeneto,pch=16,col="blue")
points(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],Comuni$Latitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],col="red",pch=16)
# Se allargo un pò la regione, ci sta dentro anche questo comune
# identify(xVeneto,yVeneto,pos=T,plot=T)
xVeneto[127]=11.26
yVeneto[127]=45.92
points(xVeneto[127],yVeneto[127],col="blue",pch=16)
points(xVeneto,yVeneto,type='l',col="blue",lwd=2)

# Ora il comune di Melara
plot(xVeneto_tot,yVeneto_tot,col="black",xlim=c(11.10,11.30),ylim=c(44.95,45.15),type='l',main="Comune di Melara",xlab="xVeneto",ylab="yVeneto")
points(xVeneto,yVeneto,type='l',col="blue")
points(xVeneto,yVeneto,pch=16,col="blue")
points(Comuni$Longitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],Comuni$Latitudine[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0],col="red",pch=16)
# Se allargo un pò la regione, ci sta dentro anche questo comune
# identify(xVeneto,yVeneto,pos=T,plot=T)
xVeneto[105]=11.19
yVeneto[105]=45.05
points(xVeneto[105],yVeneto[105],col="blue",pch=16)
xVeneto[106]=11.145
yVeneto[106]=45.10
points(xVeneto[106],yVeneto[106],col="blue",pch=16)
points(xVeneto,yVeneto,type='l',col="blue",lwd=2)

# Riprovo tutto con queste modifiche
# Controllo di intersezioni
Intersect<-Intersections(xVeneto,yVeneto)
Intersect




##### SMOOTHING SPLINES VENEZIA #####

# Ora studio l'isola di venezia.
# Anche qui devo cercare di ridurre al minimo i triangoli con vertici tutti di frontiera
xVenezia_tot=xbound[labels==6]
yVenezia_tot=ybound[labels==6]
s<-0
previous<-0

for(i in 2:length(xVenezia_tot))
{
    # Aggiorno la distanza
    previous<-previous+sqrt((xVenezia_tot[i]-xVenezia_tot[i-1])^2+(yVenezia_tot[i]-yVenezia_tot[i-1])^2)
    s<-c(s,previous)  
}
N<-15
# Nel range di ascissa curvilinea fisso N punti
Val<-seq(0,max(s),by=max(s)/N)
# Decido di usare delle splines cubiche
m<-3

# Scelgo come upper bound il numero di punti della valutazione per il numero di basi con
# cui comporre la frontiera
# Come scelgo l'ottimo del numero di basi?
# Mi interessa che sia minimizzata la distanza totale quadratica euclidea
# poligono reale di Venezia
if(do_cycles==TRUE)
{
    bestsum<-9999999999999
    best<-0
    sum_vect<-NULL
    best_vect<-NULL
    # Qui mi salvo un po' di risultati numerici
    sink(file = "ROutputVenezia.txt", append = FALSE)
    print("nbasis sum best Intersections Trinagoli")
    for(nbasis in (m+1):N)
    {
        basis = create.bspline.basis(c(0,max(s)), nbasis, m)
        xsmooth = smooth.basis(argvals=s, y=xVenezia_tot, basis)
        ysmooth = smooth.basis(argvals=s, y=yVenezia_tot, basis)
        xVenezia  = eval.fd(Val, xsmooth$fd)
        yVenezia  = eval.fd(Val, ysmooth$fd)
        png(filename = paste(nbasis,".png",sep=""))
        plot(xVenezia_tot,yVenezia_tot,type='l',main=paste(nbasis,sep=""))
        points(xVenezia,yVenezia,col="blue",type='l')
        dev.off()
        sum=0
        for(i in 1:N)
        {
            sum<-sum+dist2Line(c(xVenezia[i],yVenezia[i]), cbind(xVenezia_tot,yVenezia_tot), distfun=SquareEuclideanDistance)[1]
        }
        if(sum<bestsum)
        {
            best<-nbasis
            bestsum=sum
        }
        Intersect<-Intersections(xVenezia,yVenezia)
        
        # Oggetti completi
        x<-c(Comuni$Longitudine[Comuni$Codice==425],xVenezia)
        y<-c(Comuni$Latitudine[Comuni$Codice==425],yVenezia)
        # Creo i Boundaries
        Boundaries<-NULL
        for(i in 2:(length(x)-1))
        {
            Boundaries<-rbind(Boundaries, c(i,i+1))
        }
        Boundaries<-rbind(Boundaries, c(length(x),2))
        # Ora triangolazione
        # Oggetto pslg
        pslg_obj<-pslg(cbind(x,y),S=Boundaries)
        # Creo la mesh
        # Y dice di non aggiungere Steiner Points
        # D dice di triangolare con Delaunay
        mesh<-triangulate(pslg_obj,Y=TRUE,D=TRUE)
        # Estrazione dei triangoli
        Triang<-mesh$T
        # Plot della triangolazione
        png(filename = paste("Triangolazione a Venezia", nbasis, ".png", sep=""))
        plot(x,y,col="white")
        for (ne in 1:dim(Triang)[1])
        {
            polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
        }
        # Controllo se la triangolazione ha triangoli con solo punti di bordo
        BorderTR<-BorderTriangles(mesh$T,Boundaries)
        BorderTR
        # Li coloro
        for (ne in BorderTR)
        {
            polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]),col="green")
        }
        # Ci sono comuni esterni alla frontiera?
        points(Comuni$Longitudine[Comuni$Codice==425],Comuni$Latitudine[Comuni$Codice==425],col="red",pch=16)
        
        print(paste(
            nbasis,
            sum,
            best,
            dim(Intersect)[1],
            length(BorderTR),
            sep=" "))
        sum_vect<-c(sum_vect,sum)
        best_vect<-c(best_vect,best)
        dev.off()
    }
    sink()
    png(filename = "Sum.png")
    plot(sum_vect,type='l')
    dev.off()
    png(filename = "Best.png")
    plot(best_vect,type='l')
    dev.off()
}

# Troppi punti fanno in modo che si crei sempre un elevato numero di triangoli da eliminare
# Quindi devo fare qualcosa di più semmplice
# Scelgo nbasis 11
nbasis<-11
basis = create.bspline.basis(c(0,max(s)), nbasis, m)
xsmooth = smooth.basis(argvals=s, y=xVenezia_tot, basis)
ysmooth = smooth.basis(argvals=s, y=yVenezia_tot, basis)
xVenezia  = eval.fd(Val, xsmooth$fd)
yVenezia  = eval.fd(Val, ysmooth$fd)
# Traccio ora i risultati
png(filename = "Venezia - Smoothing x (regression splines).png")
plot(s,xVenezia_tot,main="Smoothing in x",xlab="Ascissa curvilinea",ylab="xVenezia")
points(Val,xVenezia,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Venezia - Smoothing y (regression splines).png")
plot(s,yVenezia_tot,main="Smoothing in y",xlab="Ascissa curvilinea",ylab="yVenezia")
points(Val,yVenezia,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Venezia - Smoothing regione (regression splines).png")
plot(xVenezia_tot,yVenezia_tot,type='l',main="Nuova definizione della regione",xlab="xVenezia",ylab="yVenezia")
points(xVenezia,yVenezia,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()

# Controllo di intersezioni
Intersect<-Intersections(xVenezia,yVenezia)
Intersect

# Devo modificare questa visione dell'isola
plot(xVenezia,yVenezia,type='l')
# identify(xVenezia,yVenezia)
# Alcun li devo modificare
xVenezia[3]=12.329
yVenezia[3]=45.450

xVenezia[4]=12.35
yVenezia[4]=45.439

xVenezia[5]=12.354
yVenezia[5]=45.438

xVenezia[6]=12.362
yVenezia[6]=45.441

xVenezia[7]=12.367
yVenezia[7]=45.425

xVenezia[8]=12.363
yVenezia[8]=45.4245

xVenezia[9]=12.344
yVenezia[9]=45.4342

xVenezia[10]=12.333
yVenezia[10]=45.428

xVenezia[11]=12.317
yVenezia[11]=45.431

xVenezia[12]=12.312
yVenezia[12]=45.431

xVenezia[13]=12.304
yVenezia[13]=45.43462

# Preparo il raccordo
# Troppo largo il settimo punto di Venezia
xVenezia[2]=12.3134
yVenezia[2]=45.4460

xVenezia[16]=12.3124
yVenezia[16]=45.4452

# Tolgo il primo punto
xVenezia<-xVenezia[-c(1,15)]
yVenezia<-yVenezia[-c(1,15)]
plot(xVenezia,yVenezia,type='l')

# Unisco venezia alla regione totale
plot(xVeneto,yVeneto,type='l',xlim=c(12.27,12.33),ylim=c(45.42,45.49))
points(xVenezia,yVenezia,type='l')

points(12.280237,45.463851)
points(12.280913,45.464272)
# Così va bene

# Devo però avere più punti all'interno del ponte, in modo da avere più triangoli
m=(45.464272-yVenezia[1])/(12.280913-xVenezia[1])
q=(yVenezia[1]-m*xVenezia[1])
provax1<-12.29
provay1<-m*provax1+q
provax2<-12.30
provay2<-m*provax2+q
points(provax1,provay1)
points(provax2,provay2)
m=(45.463851-yVenezia[length(xVenezia)])/(12.280237-xVenezia[length(xVenezia)])
q=(yVenezia[length(xVenezia)]-m*xVenezia[length(xVenezia)])
provax3<-12.29
provay3<-m*provax3+q
provax4<-12.30
provay4<-m*provax4+q
points(provax3,provay3)
points(provax4,provay4)
xVeneto=c(xVeneto[1:54],12.280913,provax1,provax2,xVenezia,provax4,provax3,12.280237,xVeneto[55:length(xVeneto)])
yVeneto=c(yVeneto[1:54],45.464272,provay1,provay2,yVenezia,provay4,provay3,45.463851,yVeneto[55:length(yVeneto)])

plot(xVeneto,yVeneto,type='l',xlim=c(12.29,12.33),ylim=c(45.42,45.49))

# Perfettamente unite





##### CHIOGGIA #####

# Ora studio l'isola di chioggia. Non la descrivo troppo bene però.
# In questo caso mi basta "inglobarla" alla regione
# Anche qui devo cercare di ridurre al minimo i triangoli con vertici tutti di frontiera

plot(xVeneto,yVeneto,type='l',xlim=c(12.24,12.32),ylim=c(45.1,45.3))

xChioggia1=xbound[labels==16]
yChioggia1=ybound[labels==16]

points(xChioggia1,yChioggia1,type='l')

xChioggia2=xbound[labels==4]
yChioggia2=ybound[labels==4]

points(xChioggia2,yChioggia2,type='l')

#In pratica occorre modificare il confine tra i punt 83 e 85 perchè
#il punto 80 è finito dentro chioggia

points(12.26,45.19)
points(12.27,45.23)
points(12.31,45.235)
points(12.31,45.16)

xVeneto=c(xVeneto[1:90],12.26,12.275,12.305,12.31,xVeneto[92:length(xVeneto)])
yVeneto=c(yVeneto[1:90],45.19,45.225,45.232,45.16,yVeneto[92:length(yVeneto)])

plot(xVeneto,yVeneto,type='l',xlim=c(12.24,12.32),ylim=c(45.1,45.3))
# Perfetto



##### REVISIONE DELLA ZONA DI CAVALLINO-TREPORTI #####

# Non è per nulla aderente alla realtà
plot(xVeneto,yVeneto,type='l',xlim=c(12.4,12.55),ylim=c(45.4,45.5))
points(xVeneto_tot,yVeneto_tot,type='l')

yVeneto[35]=45.42
xVeneto[36]=12.41892
yVeneto[36]=45.44426
xVeneto[37]=12.45
yVeneto[37]=45.47453

points(Comuni$Longitudine,Comuni$Latitudine,col="red",pch=16)
# Devo spostare un pò anche cavallino-treporti
Comuni$Longitudine[581]=12.465
Comuni$Latitudine[581]=45.455




##### REVISIONE DELLA LAGUNA #####


# Non è per nulla aderente alla realtà di GoogleMaps
plot(xVeneto,yVeneto,type='l',xlim=c(12.3,12.6),ylim=c(45.46,45.56))
points(xVeneto_tot,yVeneto_tot,type='l')
points(Comuni$Longitudine,Comuni$Latitudine,col="red",pch=16)



# Devo andare a modificare il punto 35
xVeneto[48]=xVeneto[49]
yVeneto[48]=45.54

xVeneto=xVeneto[-c(49:52)]
yVeneto=yVeneto[-c(49:52)]

plot(xVeneto,yVeneto,type='l',xlim=c(12.3,12.6),ylim=c(45.46,45.56))
points(xVeneto_tot,yVeneto_tot,type='l')
points(Comuni$Longitudine,Comuni$Latitudine,col="red",pch=16)

# Tolgo anche altri due punti, che non sono importanti per la regione
# identify(xVeneto,yVeneto)

xVeneto=xVeneto[-c(43:44)]
yVeneto=yVeneto[-c(43:44)]

plot(xVeneto,yVeneto,type='l',xlim=c(12.3,12.6),ylim=c(45.46,45.56))
points(xVeneto_tot,yVeneto_tot,type='l')
points(Comuni$Longitudine,Comuni$Latitudine,col="red",pch=16)

# Fatto




##### AGGIUNTA DELL'ISOLA DI LIDO #####

# Ora studio l'isola di lido.
# Anche qui devo cercare di ridurre al minimo i triangoli con vertici tutti di frontiera
xLido_tot=xbound[labels==5]
yLido_tot=ybound[labels==5]
s<-0
previous<-0

for(i in 2:length(xLido_tot))
{
    # Aggiorno la distanza
    previous<-previous+sqrt((xLido_tot[i]-xLido_tot[i-1])^2+(yLido_tot[i]-yLido_tot[i-1])^2)
    s<-c(s,previous)  
}
N<-15
# Nel range di ascissa curvilinea fisso N punti
Val<-seq(0,max(s),by=max(s)/N)
# Decido di usare delle splines cubiche
m<-4

# Scelgo come upper bound il numero di punti della valutazione per il numero di basi con
# cui comporre la frontiera
# Come scelgo l'ottimo del numero di basi?
# Mi interessa che sia minimizzata la distanza totale quadratica euclidea
# poligono reale di Lido
if(do_cycles==TRUE)
{
    bestsum<-9999999999999
    best<-0
    sum_vect<-NULL
    best_vect<-NULL
    # Qui mi salvo un po' di risultati numerici
    sink(file = "ROutputLido.txt", append = FALSE)
    print("nbasis sum best Intersections Trinagoli")
    for(nbasis in (m+1):N)
    {
        basis = create.bspline.basis(c(0,max(s)), nbasis, m)
        xsmooth = smooth.basis(argvals=s, y=xLido_tot, basis)
        ysmooth = smooth.basis(argvals=s, y=yLido_tot, basis)
        xLido  = eval.fd(Val, xsmooth$fd)
        yLido  = eval.fd(Val, ysmooth$fd)
        png(filename = paste(nbasis,".png",sep=""))
        plot(xLido_tot,yLido_tot,type='l',main=paste(nbasis,sep=""))
        points(xLido,yLido,col="blue",type='l')
        dev.off()
        sum=0
        for(i in 1:N)
        {
            sum<-sum+dist2Line(c(xLido[i],yLido[i]), cbind(xLido_tot,yLido_tot), distfun=SquareEuclideanDistance)[1]
        }
        if(sum<bestsum)
        {
            best<-nbasis
            bestsum=sum
        }
        Intersect<-Intersections(xLido,yLido)
        
        # Oggetti completi
        x<-c(Comuni$Longitudine[Comuni$Codice==425],xLido)
        y<-c(Comuni$Latitudine[Comuni$Codice==425],yLido)
        # Creo i Boundaries
        Boundaries<-NULL
        for(i in 2:(length(x)-1))
        {
            Boundaries<-rbind(Boundaries, c(i,i+1))
        }
        Boundaries<-rbind(Boundaries, c(length(x),2))
        # Ora triangolazione
        # Oggetto pslg
        pslg_obj<-pslg(cbind(x,y),S=Boundaries)
        # Creo la mesh
        # Y dice di non aggiungere Steiner Points
        # D dice di triangolare con Delaunay
        mesh<-triangulate(pslg_obj,Y=TRUE,D=TRUE)
        # Estrazione dei triangoli
        Triang<-mesh$T
        # Plot della triangolazione
        png(filename = paste("Triangolazione a Lido", nbasis, ".png", sep=""))
        plot(x,y,col="white")
        for (ne in 1:dim(Triang)[1])
        {
            polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
        }
        # Controllo se la triangolazione ha triangoli con solo punti di bordo
        BorderTR<-BorderTriangles(mesh$T,Boundaries)
        BorderTR
        # Li coloro
        for (ne in BorderTR)
        {
            polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]),col="green")
        }
        # Ci sono comuni esterni alla frontiera?
        points(Comuni$Longitudine[Comuni$Codice==425],Comuni$Latitudine[Comuni$Codice==425],col="red",pch=16)
        
        print(paste(
            nbasis,
            sum,
            best,
            dim(Intersect)[1],
            length(BorderTR),
            sep=" "))
        sum_vect<-c(sum_vect,sum)
        best_vect<-c(best_vect,best)
        dev.off()
    }
    sink()
    png(filename = "Sum.png")
    plot(sum_vect,type='l')
    dev.off()
    png(filename = "Best.png")
    plot(best_vect,type='l')
    dev.off()
}

# Scelgo nbasis 8
nbasis<-8
basis = create.bspline.basis(c(0,max(s)), nbasis, m)
xsmooth = smooth.basis(argvals=s, y=xLido_tot, basis)
ysmooth = smooth.basis(argvals=s, y=yLido_tot, basis)
xLido  = eval.fd(Val, xsmooth$fd)
yLido  = eval.fd(Val, ysmooth$fd)
# Traccio ora i risultati
png(filename = "Lido - Smoothing x (regression splines).png")
plot(s,xLido_tot,main="Smoothing in x",xlab="Ascissa curvilinea",ylab="xLido")
points(Val,xLido,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Lido - Smoothing y (regression splines).png")
plot(s,yLido_tot,main="Smoothing in y",xlab="Ascissa curvilinea",ylab="yLido")
points(Val,yLido,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Lido - Smoothing regione (regression splines).png")
plot(xLido_tot,yLido_tot,type='l',main="Nuova definizione della regione",xlab="xLido",ylab="yLido")
points(xLido,yLido,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()


# Controllo di intersezioni
Intersect<-Intersections(xLido,yLido)
Intersect

# Devo modificare questa visione dell'isola
plot(xLido,yLido,type='l')
# identify(xLido,yLido)
# Alcun li devo modificare



xLido[9]=12.331
yLido[9]=45.35

xLido[10]=12.33
yLido[10]=45.337

xLido[11]=12.308
yLido[11]=45.340

xLidoo[4]=12.40166
yLido[4]=45.43

xLido<-xLido[-c(1:3)]
yLido<-yLido[-c(1:3)]

# Unisco ora all'isola di Venezia

plot(xLido,yLido,type='l')
points(xVeneto,yVeneto,type='l')

#Unisco in questo modo
xVeneto<-c(xVeneto[1:56],12.3665,12.3805,xLido[length(xLido)],xLido[1:(length(xLido)-1)],12.38,xVeneto[57:length(xVeneto)])
yVeneto<-c(yVeneto[1:56],45.427,45.429,yLido[length(yLido)],yLido[1:(length(yLido)-1)],45.427,yVeneto[57:length(yVeneto)])



##### AGGIUNTA DELL'ISOLA DI PELLESTRINA #####

# Ora studio l'isola di Pellestrina.
# Anche qui devo cercare di ridurre al minimo i triangoli con vertici tutti di frontiera
xPellestrina_tot=xbound[labels==11]
yPellestrina_tot=ybound[labels==11]
s<-0
previous<-0

for(i in 2:length(xPellestrina_tot))
{
    # Aggiorno la distanza
    previous<-previous+sqrt((xPellestrina_tot[i]-xPellestrina_tot[i-1])^2+(yPellestrina_tot[i]-yPellestrina_tot[i-1])^2)
    s<-c(s,previous)  
}
N<-15
# Nel range di ascissa curvilinea fisso N punti
Val<-seq(0,max(s),by=max(s)/N)
# Decido di usare delle splines cubiche
m<-4

# Scelgo come upper bound il numero di punti della valutazione per il numero di basi con
# cui comporre la frontiera
# Come scelgo l'ottimo del numero di basi?
# Mi interessa che sia minimizzata la distanza totale quadratica euclidea
# poligono reale di Pellestrina
if(do_cycles==TRUE)
{
    bestsum<-9999999999999
    best<-0
    sum_vect<-NULL
    best_vect<-NULL
    # Qui mi salvo un po' di risultati numerici
    sink(file = "ROutputPellestrina.txt", append = FALSE)
    print("nbasis sum best Intersections Trinagoli")
    for(nbasis in (m+1):N)
    {
        basis = create.bspline.basis(c(0,max(s)), nbasis, m)
        xsmooth = smooth.basis(argvals=s, y=xPellestrina_tot, basis)
        ysmooth = smooth.basis(argvals=s, y=yPellestrina_tot, basis)
        xPellestrina  = eval.fd(Val, xsmooth$fd)
        yPellestrina  = eval.fd(Val, ysmooth$fd)
        png(filename = paste(nbasis,".png",sep=""))
        plot(xPellestrina_tot,yPellestrina_tot,type='l',main=paste(nbasis,sep=""))
        points(xPellestrina,yPellestrina,col="blue",type='l')
        dev.off()
        sum=0
        for(i in 1:N)
        {
            sum<-sum+dist2Line(c(xPellestrina[i],yPellestrina[i]), cbind(xPellestrina_tot,yPellestrina_tot), distfun=SquareEuclideanDistance)[1]
        }
        if(sum<bestsum)
        {
            best<-nbasis
            bestsum=sum
        }
        Intersect<-Intersections(xPellestrina,yPellestrina)
        
        # Oggetti completi
        x<-c(Comuni$Longitudine[Comuni$Codice==425],xPellestrina)
        y<-c(Comuni$Latitudine[Comuni$Codice==425],yPellestrina)
        # Creo i Boundaries
        Boundaries<-NULL
        for(i in 2:(length(x)-1))
        {
            Boundaries<-rbind(Boundaries, c(i,i+1))
        }
        Boundaries<-rbind(Boundaries, c(length(x),2))
        # Ora triangolazione
        # Oggetto pslg
        pslg_obj<-pslg(cbind(x,y),S=Boundaries)
        # Creo la mesh
        # Y dice di non aggiungere Steiner Points
        # D dice di triangolare con Delaunay
        mesh<-triangulate(pslg_obj,Y=TRUE,D=TRUE)
        # Estrazione dei triangoli
        Triang<-mesh$T
        # Plot della triangolazione
        png(filename = paste("Triangolazione a Pellestrina", nbasis, ".png", sep=""))
        plot(x,y,col="white")
        for (ne in 1:dim(Triang)[1])
        {
            polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
        }
        # Controllo se la triangolazione ha triangoli con solo punti di bordo
        BorderTR<-BorderTriangles(mesh$T,Boundaries)
        BorderTR
        # Li coloro
        for (ne in BorderTR)
        {
            polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]),col="green")
        }
        # Ci sono comuni esterni alla frontiera?
        points(Comuni$Longitudine[Comuni$Codice==425],Comuni$Latitudine[Comuni$Codice==425],col="red",pch=16)
        
        print(paste(
            nbasis,
            sum,
            best,
            dim(Intersect)[1],
            length(BorderTR),
            sep=" "))
        sum_vect<-c(sum_vect,sum)
        best_vect<-c(best_vect,best)
        dev.off()
    }
    sink()
    png(filename = "Sum.png")
    plot(sum_vect,type='l')
    dev.off()
    png(filename = "Best.png")
    plot(best_vect,type='l')
    dev.off()
}

# Scelgo nbasis 9
nbasis<-9
basis = create.bspline.basis(c(0,max(s)), nbasis, m)
xsmooth = smooth.basis(argvals=s, y=xPellestrina_tot, basis)
ysmooth = smooth.basis(argvals=s, y=yPellestrina_tot, basis)
xPellestrina  = eval.fd(Val, xsmooth$fd)
yPellestrina  = eval.fd(Val, ysmooth$fd)
# Traccio ora i risultati
png(filename = "Pellestrina - Smoothing x (regression splines).png")
plot(s,xPellestrina_tot,main="Smoothing in x",xlab="Ascissa curvilinea",ylab="xPellestrina")
points(Val,xPellestrina,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Pellestrina - Smoothing y (regression splines).png")
plot(s,yPellestrina_tot,main="Smoothing in y",xlab="Ascissa curvilinea",ylab="yPellestrina")
points(Val,yPellestrina,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Pellestrina - Smoothing regione (regression splines).png")
plot(xPellestrina_tot,yPellestrina_tot,type='l',main="Nuova definizione della regione",xlab="xPellestrina",ylab="yPellestrina")
points(xPellestrina,yPellestrina,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()


# Controllo di intersezioni
Intersect<-Intersections(xPellestrina,yPellestrina)
Intersect

# Devo modificare questa visione dell'isola
# identify(xPellestrina,yPellestrina)
# Alcun li devo modificare

xPellestrina[1]=12.331
yPellestrina[1]=45.332

xPellestrina[15]=12.32117
yPellestrina[15]=45.330

xPellestrina[16]=12.31594
yPellestrina[16]=45.331

xPellestrina[8]=12.30018
yPellestrina[8]=45.25

xPellestrina[9]=12.2995
yPellestrina[9]=45.255

xPellestrina[9]=12.2995
yPellestrina[9]=45.255

xPellestrina=xPellestrina[-11]
yPellestrina=yPellestrina[-11]


#Cerco due punti sul lato basso di Lido
m=(yVeneto[67]-yVeneto[66])/(xVeneto[67]-xVeneto[66])
q=yVeneto[66]-m*xVeneto[66]

xp2<-12.31
yp2<-m*xp2+q
xp1<-12.3105
yp1<-m*xp1+q


xVeneto<-c(xVeneto[1:66],xp1,12.3161,xPellestrina,xp2,xVeneto[67:length(xVeneto)])
yVeneto<-c(yVeneto[1:66],yp1,45.332,yPellestrina,yp2,yVeneto[67:length(yVeneto)])
plot(xPellestrina,yPellestrina,type='l',xlim=c(12.30,12.33),ylim=c(45.23,45.36))
points(xVeneto,yVeneto,type='l')

xVeneto<-xVeneto[-74]
yVeneto<-yVeneto[-74]

# Controllo se ci sono dei comuni che rimangono fuori dalla frontiera
PolyPoints<-cbind(xVeneto,yVeneto)
if(sum(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip)==length(Comuni$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}
# Quali comuni restano fuori?
Comuni$Comune[pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip==0]
# Infatti devo ancora fare murano




##### AGGIUNTA DELL'ISOLA DI MURANO #####

# Ora studio l'isola di Murano.
# Anche qui devo cercare di ridurre al minimo i triangoli con vertici tutti di frontiera
xMurano_tot=xbound[labels==12]
yMurano_tot=ybound[labels==12]
s<-0
previous<-0

for(i in 2:length(xMurano_tot))
{
    # Aggiorno la distanza
    previous<-previous+sqrt((xMurano_tot[i]-xMurano_tot[i-1])^2+(yMurano_tot[i]-yMurano_tot[i-1])^2)
    s<-c(s,previous)  
}
N<-10
# Nel range di ascissa curvilinea fisso N punti
Val<-seq(0,max(s),by=max(s)/N)
# Decido di usare delle splines cubiche
m<-4

# Scelgo come upper bound il numero di punti della valutazione per il numero di basi con
# cui comporre la frontiera
# Come scelgo l'ottimo del numero di basi?
# Mi interessa che sia minimizzata la distanza totale quadratica euclidea
# poligono reale di Murano
if(do_cycles==TRUE)
{
    bestsum<-9999999999999
    best<-0
    sum_vect<-NULL
    best_vect<-NULL
    # Qui mi salvo un po' di risultati numerici
    sink(file = "ROutputMurano.txt", append = FALSE)
    print("nbasis sum best Intersections Trinagoli")
    for(nbasis in (m+1):N)
    {
        basis = create.bspline.basis(c(0,max(s)), nbasis, m)
        xsmooth = smooth.basis(argvals=s, y=xMurano_tot, basis)
        ysmooth = smooth.basis(argvals=s, y=yMurano_tot, basis)
        xMurano  = eval.fd(Val, xsmooth$fd)
        yMurano  = eval.fd(Val, ysmooth$fd)
        png(filename = paste(nbasis,".png",sep=""))
        plot(xMurano_tot,yMurano_tot,type='l',main=paste(nbasis,sep=""))
        points(xMurano,yMurano,col="blue",type='l')
        dev.off()
        sum=0
        for(i in 1:N)
        {
            sum<-sum+dist2Line(c(xMurano[i],yMurano[i]), cbind(xMurano_tot,yMurano_tot), distfun=SquareEuclideanDistance)[1]
        }
        if(sum<bestsum)
        {
            best<-nbasis
            bestsum=sum
        }
        Intersect<-Intersections(xMurano,yMurano)
        
        # Oggetti completi
        x<-c(Comuni$Longitudine[Comuni$Codice==425],xMurano)
        y<-c(Comuni$Latitudine[Comuni$Codice==425],yMurano)
        # Creo i Boundaries
        Boundaries<-NULL
        for(i in 2:(length(x)-1))
        {
            Boundaries<-rbind(Boundaries, c(i,i+1))
        }
        Boundaries<-rbind(Boundaries, c(length(x),2))
        # Ora triangolazione
        # Oggetto pslg
        pslg_obj<-pslg(cbind(x,y),S=Boundaries)
        # Creo la mesh
        # Y dice di non aggiungere Steiner Points
        # D dice di triangolare con Delaunay
        mesh<-triangulate(pslg_obj,Y=TRUE,D=TRUE)
        # Estrazione dei triangoli
        Triang<-mesh$T
        # Plot della triangolazione
        png(filename = paste("Triangolazione a Murano", nbasis, ".png", sep=""))
        plot(x,y,col="white")
        for (ne in 1:dim(Triang)[1])
        {
            polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
        }
        # Controllo se la triangolazione ha triangoli con solo punti di bordo
        BorderTR<-BorderTriangles(mesh$T,Boundaries)
        BorderTR
        # Li coloro
        for (ne in BorderTR)
        {
            polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]),col="green")
        }
        # Ci sono comuni esterni alla frontiera?
        points(Comuni$Longitudine[Comuni$Codice==425],Comuni$Latitudine[Comuni$Codice==425],col="red",pch=16)
        
        print(paste(
            nbasis,
            sum,
            best,
            dim(Intersect)[1],
            length(BorderTR),
            sep=" "))
        sum_vect<-c(sum_vect,sum)
        best_vect<-c(best_vect,best)
        dev.off()
    }
    sink()
    png(filename = "Sum.png")
    plot(sum_vect,type='l')
    dev.off()
    png(filename = "Best.png")
    plot(best_vect,type='l')
    dev.off()
}

# Scelgo nbasis 9
nbasis<-7
basis = create.bspline.basis(c(0,max(s)), nbasis, m)
xsmooth = smooth.basis(argvals=s, y=xMurano_tot, basis)
ysmooth = smooth.basis(argvals=s, y=yMurano_tot, basis)
xMurano  = eval.fd(Val, xsmooth$fd)
yMurano  = eval.fd(Val, ysmooth$fd)
# Traccio ora i risultati
png(filename = "Murano - Smoothing x (regression splines).png")
plot(s,xMurano_tot,main="Smoothing in x",xlab="Ascissa curvilinea",ylab="xMurano")
points(Val,xMurano,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Murano - Smoothing y (regression splines).png")
plot(s,yMurano_tot,main="Smoothing in y",xlab="Ascissa curvilinea",ylab="yMurano")
points(Val,yMurano,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()
png(filename = "Murano - Smoothing regione (regression splines).png")
plot(xMurano_tot,yMurano_tot,type='l',main="Nuova definizione della regione",xlab="xMurano",ylab="yMurano")
points(xMurano,yMurano,col="blue",type='l')
legend("bottomleft", legend=c("Reale","Smmothing"), col=c("black","blue"), lty=1)
dev.off()


# Controllo di intersezioni
Intersect<-Intersections(xMurano,yMurano)
Intersect

# Devo modificare questa visione dell'isola
plot(xMurano,yMurano,type='l')
# identify(xMurano,yMurano)
# Alcun li devo modificare

xMurano[7]=12.348
yMurano[7]=45.451

yMurano[8]=45.454

xMurano[9]=12.345
yMurano[9]=45.4635

# Unisco al resto della regione
xVeneto<-c(xVeneto[1:53],12.3405,12.3477,xMurano[8:length(xMurano)],xMurano[1:7],12.3408,xVeneto[54:length(xVeneto)])
yVeneto<-c(yVeneto[1:53],45.4437,45.4513,yMurano[8:length(yMurano)],yMurano[1:7],45.4435,yVeneto[54:length(yVeneto)])

# Controllo se ci sono dei comuni che rimangono fuori dalla frontiera
PolyPoints<-cbind(xVeneto,yVeneto)
if(sum(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip)==length(Comuni$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}

# Salvo tutto, per la triangolazione della provincia di venezia
save(file="PerPV.RData",xVeneto,yVeneto)



##### PLOT DEI TRIANGOLI CON PUNTI DI BORDO #####

x<-c(Comuni$Longitudine,xVeneto)
y<-c(Comuni$Latitudine,yVeneto)
#Creo i Boundaries
Boundaries<-NULL
for(i in (length(Comuni$Longitudine)+1):(length(x)-1))
{
    Boundaries<-rbind(Boundaries, c(i,i+1))
}
Boundaries<-rbind(Boundaries, c(length(x),length(Comuni$Longitudine)+1))
#Ora triangolazione
#Oggetto pslg
pslg_obj<-pslg(cbind(x,y),S=Boundaries)
#Creo la mesh
#Y dice di non aggiungere Steiner Points
#D dice di triangolare con Delaunay
mesh<-triangulate(pslg_obj,Y=TRUE,D=TRUE)
#Estrazione dei triangoli
Triang<-mesh$T

# Triangoli con soli punti di Bordo
BorderTR<-BorderTriangles(mesh$T,Boundaries)
BorderTR

#Controllo i baricentri dei triangoli
xG<-NULL
yG<-NULL
for (i in 1:length(BorderTR))
{
    xG<-c(xG,(x[Triang[BorderTR[i],1]]+x[Triang[BorderTR[i],2]]+x[Triang[BorderTR[i],3]])/3)
    yG<-c(yG,(y[Triang[BorderTR[i],1]]+y[Triang[BorderTR[i],2]]+y[Triang[BorderTR[i],3]])/3)
}

for(k in 1:length(xG))
{
    xlim1<-c(xG[k]-0.02,xG[k]+0.02)
    xlim2<-c(xG[k]-0.08,xG[k]+0.08)
    
    ylim1<-c(yG[k]-0.1,yG[k]+0.1)
    ylim2<-c(yG[k]-0.3,yG[k]+0.3)
    
    
    png(filename=paste("BorderTR",k," 1.png",sep=""))
    plot(x,y,col="white",xlim=xlim1,ylim=ylim1,main=paste(k))
    for (ne in 1:dim(Triang)[1])
    {
        polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
    }
    #Controllo se la triangolazione ha triangoli con solo punti di bordo
    #Li coloro
    for (ne in BorderTR)
    {
        polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]),col="green")
    }
    points(xG[k],yG[k],pch=16,col="red")
    dev.off()
    
    png(filename=paste("BorderTR",k," 2.png",sep=""))
    plot(x,y,col="white",xlim=xlim2,ylim=ylim2,main=paste(k))
    for (ne in 1:dim(Triang)[1])
    {
        polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
    }
    #Controllo se la triangolazione ha triangoli con solo punti di bordo
    #Li coloro
    for (ne in BorderTR)
    {
        polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]),col="green")
    }
    points(xG[k],yG[k],pch=16,col="red")
    dev.off()
}

# Dall'analisi dei grafici si ha che:
IDDelete=c(52,78)
IDChoose=c(14,21,39,40,44,45,46,56,77)
IDKeep=1:78
IDKeep=IDKeep[-c(IDDelete,IDChoose)]
#Ma questi sono secondo gli indici di BorderTR, ora devo trovare gli ID di Triang
IDDelete<-BorderTR[IDDelete]
IDChoose<-BorderTR[IDChoose]
IDKeep<-BorderTR[IDKeep]



##### SALVATAGGIO DEI RISULTATI SENZA ELIMINAZIONE DEI TRIANGOLI #####

x<-c(Comuni$Longitudine,xVeneto)
y<-c(Comuni$Latitudine,yVeneto)
#Creo i Boundaries
Boundaries<-NULL
for(i in (length(Comuni$Longitudine)+1):(length(x)-1))
{
    Boundaries<-rbind(Boundaries, c(i,i+1))
}
Boundaries<-rbind(Boundaries, c(length(x),length(Comuni$Longitudine)+1))
#Ora triangolazione
#Oggetto pslg
pslg_obj<-pslg(cbind(x,y),S=Boundaries)
#Creo la mesh
#Y dice di non aggiungere Steiner Points
#D dice di triangolare con Delaunay
mesh<-triangulate(pslg_obj,Y=TRUE,D=TRUE)
#Estrazione dei triangoli
Triang<-mesh$T

#Ci sono comuni esterni alla frontiera?
PolyPoints<-cbind(xVeneto,yVeneto)
if(sum(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip)==length(Comuni$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}


Codici<-Comuni$Codice
#Creo anche un vettore di codici di comune
for (i in (length(Codici)+1):length(x))
{
    Codici<-c(Codici,NA)
}

save(file="ConTriangoli.RData",x,y,Codici,Triang)









##### SALVATAGGIO DEI RISULTATI TOGLIENDO SOLO ALCUNI TRIANGOLI #####

#Creo i Boundaries
x<-c(Comuni$Longitudine,xVeneto)
y<-c(Comuni$Latitudine,yVeneto)
Boundaries<-NULL
for(i in (length(Comuni$Longitudine)+1):(length(x)-1))
{
    Boundaries<-rbind(Boundaries, c(i,i+1))
}
Boundaries<-rbind(Boundaries, c(length(x),length(Comuni$Longitudine)+1))
#Ora triangolazione
#Oggetto pslg
pslg_obj<-pslg(cbind(x,y),S=Boundaries)
#Creo la mesh
#Y dice di non aggiungere Steiner Points
#D dice di triangolare con Delaunay
mesh<-triangulate(pslg_obj,Y=TRUE,D=TRUE)
#Estrazione dei triangoli
Triang<-mesh$T

L<-CleanPoints(Triang,IDDelete,x,y)
x<-L[[1]]
y<-L[[2]]
Triang<-L[[3]]

Codici<-Comuni$Codice
#Controllo che il nuovo bordo non abbia intersezioni
Intersect<-Intersections(x[(length(Codici)+1):length(x)],y[(length(Codici)+1):length(x)])
Intersect
#Creo anche un vettore di codici di comune
for (i in (length(Codici)+1):length(x))
{
    Codici<-c(Codici,NA)
}


#Controllo che tutti i comuni siano dentro

PolyPoints<-cbind(x[is.na(Codici)],y[is.na(Codici)])
if(sum(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip)==length(Comuni$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}
#Salvo tutto
save(file="ConMenoTriangoli.RData",x,y,Codici,Triang)




##### SALVATAGGIO DEI RISULTATI TOGLIENDO TUTTI TRIANGOLI POSSIBILI #####

#Creo i Boundaries
x<-c(Comuni$Longitudine,xVeneto)
y<-c(Comuni$Latitudine,yVeneto)
Boundaries<-NULL
for(i in (length(Comuni$Longitudine)+1):(length(x)-1))
{
    Boundaries<-rbind(Boundaries, c(i,i+1))
}
Boundaries<-rbind(Boundaries, c(length(x),length(Comuni$Longitudine)+1))
#Ora triangolazione
#Oggetto pslg
pslg_obj<-pslg(cbind(x,y),S=Boundaries)
#Creo la mesh
#Y dice di non aggiungere Steiner Points
#D dice di triangolare con Delaunay
mesh<-triangulate(pslg_obj,Y=TRUE,D=TRUE)
#Estrazione dei triangoli
Triang<-mesh$T

L<-CleanPoints(Triang,c(IDDelete,IDChoose),x,y)
x<-L[[1]]
y<-L[[2]]
Triang<-L[[3]]

Codici<-Comuni$Codice
#Controllo che il nuovo bordo non abbia intersezioni
Intersect<-Intersections(x[(length(Codici)+1):length(x)],y[(length(Codici)+1):length(x)])
Intersect
#Creo anche un vettore di codici di comune
for (i in (length(Codici)+1):length(x))
{
    Codici<-c(Codici,NA)
}

#Controllo che tutti i comuni siano dentro

PolyPoints<-cbind(x[is.na(Codici)],y[is.na(Codici)])
if(sum(pnt.in.poly(cbind(Comuni$Longitudine,Comuni$Latitudine),PolyPoints)$pip)==length(Comuni$Longitudine))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}
#Salvo tutto
save(file="SenzaTriangoli.RData",x,y,Codici,Triang)