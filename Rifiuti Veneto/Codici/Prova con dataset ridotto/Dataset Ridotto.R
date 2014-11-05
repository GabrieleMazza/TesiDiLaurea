#Applicazione al dataset ridotti

#Punti di confine e di comune
load("Ridotto.RData")

#library(RTriangle)
library(deldir)

##### TRIANGOLAZIONE CON I PUNTI DI BORDO #####

ncom<-length(xcom_rid)
nbound<-length(xbound_rid)

#Identifico per la funzione di RTriangle pslg quali sono di confine e quali interni
# 0 --> interni
# 1 --> confine

IDcom<-rep(0,ncom)
IDbound<-rep(1,nbound)

#Creazione degli oggetti definitivi

x_rid<-c(xcom_rid,xbound_rid)
y_rid<-c(ycom_rid,ybound_rid)
coord_rid<-cbind(x_rid,y_rid)

#Devo creare la matrice di confine per pslg
Boundaries<-NULL
for(i in (ncom+1):(ncom+nbound-1))
{
    Boundaries<-rbind(Boundaries,c(i,i+1))
}
Boundaries<-rbind(Boundaries,c(ncom+nbound,ncom+1))


#CON RTRIANGLE
#Creo l'oggetto
#pslgobj<-pslg(coord_rid[,1:2],S=Boundaries)
#Warning sul controllo degli NA, ma non c'è problema...
#Ora triangolazione
#mesh_rid<-triangulate(pslgobj, Y=TRUE,D=TRUE)
#plot(mesh_rid)
#points(coord_rid[,1:2],col=coord_rid[,3])
#T_rid<-mesh_rid$T
#dim(T_rid)

#CON DELDIR
tr<-deldir(x_rid,y_rid,plotit=TRUE)
T_rid<-triMat(tr)
source("Functions.R")
T_rid<-CleanTriangulation(x_rid, y_rid, T_rid, Boundaries)

#Plot
plot(coord_rid[,1],coord_rid[,2],type="n")
polygon(c(coord_rid[T_rid[1,1],1],coord_rid[T_rid[1,2],1],coord_rid[T_rid[1,3],1]),c(coord_rid[T_rid[1,1],2],coord_rid[T_rid[1,2],2],coord_rid[T_rid[1,3],2]))
for (ne in 1:dim(T_rid)[1])
{
    polygon(c(coord_rid[T_rid[ne,1],1],coord_rid[T_rid[ne,2],1],coord_rid[T_rid[ne,3],1]),c(coord_rid[T_rid[ne,1],2],coord_rid[T_rid[ne,2],2],coord_rid[T_rid[ne,3],2]))
}

save(file="T_RidBoundaries.RData",T_rid)


##### TRIANGOLAZIONE SENZA I PUNTI DI BORDO #####

#Ora uso solo i punti di comune

#CON DELDIR
tr<-deldir(xcom_rid,ycom_rid,plotit=TRUE)
T_rid<-triMat(tr)
#Triangolazione dell'inviluppo convesso

#CON RTRIANGLE
#coord_rid<-cbind(xcom_rid,ycom_rid)
#pslgobj<-pslg(coord_rid)
#mesh_rid<-triangulate(pslgobj,Y=TRUE,D=TRUE)
#plot(mesh_rid)
#T_rid<-mesh_rid$T

#Plot
plot(coord_rid[,1],coord_rid[,2],type="n")
polygon(c(coord_rid[T_rid[1,1],1],coord_rid[T_rid[1,2],1],coord_rid[T_rid[1,3],1]),c(coord_rid[T_rid[1,1],2],coord_rid[T_rid[1,2],2],coord_rid[T_rid[1,3],2]))
for (ne in 1:dim(T_rid)[1])
{
    polygon(c(coord_rid[T_rid[ne,1],1],coord_rid[T_rid[ne,2],1],coord_rid[T_rid[ne,3],1]),c(coord_rid[T_rid[ne,1],2],coord_rid[T_rid[ne,2],2],coord_rid[T_rid[ne,3],2]))
}

save(file="T_RidComuni.RData",T_rid)


##### APPLICAZIONE CON I PUNTI DI BORDO #####

library(fda)
library(rgl)
source("2013_SSR_AllFunctions.R")
load("T_RidBoundaries.RData")

#Traccio i punti

plot(coord_rid[,1],coord_rid[,2])

#Traccio la mesh

plot(coord_rid[,1],coord_rid[,2],type="n")
polygon(c(coord_rid[T_rid[1,1],1],coord_rid[T_rid[1,2],1],coord_rid[T_rid[1,3],1]),c(coord_rid[T_rid[1,1],2],coord_rid[T_rid[1,2],2],coord_rid[T_rid[1,3],2]))
for (ne in 1:dim(T_rid)[1])
{
    polygon(c(coord_rid[T_rid[ne,1],1],coord_rid[T_rid[ne,2],1],coord_rid[T_rid[ne,3],1]),c(coord_rid[T_rid[ne,1],2],coord_rid[T_rid[ne,2],2],coord_rid[T_rid[ne,3],2]))
}

#No alla edge matrix

e = NULL

#Ordine

order=1

basisobj = create.FEM.basis(coord_rid, e, T_rid, order)

#Creo l'oggetto fd

Rifiutifd = fd(numeric(basisobj$nbasis),basisobj)

#Risposta

data = matrix(0,nrow=dim(coord_rid)[1],ncol=2)
data[,1] = 1:dim(coord_rid)[1]
data[1:length(TotC_rid),2] = TotC_rid
for (i in length(TotC_rid):dim(coord_rid)[1])
{
    data[i,2] = 400
}


#Applico il modello

lambda = 10^(2)
Rifiuti = smooth.FEM.fd(data,Rifiutifd,lambda)

#Stima nei punti

fhat=Rifiuti$felsplobj$coef
fhat

#Plot

plot.FEM(Rifiuti$felsplobj)


##### APPLICAZIONE SENZA I PUNTI DI BORDO #####

library(fda)
library(rgl)
source("2013_SSR_AllFunctions.R")
load("T_RidComuni.RData")

#Le coordinate sono disponibili, devo unirle
coord_rid<-cbind(xcom_rid,ycom_rid)

#Traccio i punti

plot(coord_rid[,1],coord_rid[,2])

#Traccio la mesh

plot(coord_rid[,1],coord_rid[,2],type="n")
polygon(c(coord_rid[T_rid[1,1],1],coord_rid[T_rid[1,2],1],coord_rid[T_rid[1,3],1]),c(coord_rid[T_rid[1,1],2],coord_rid[T_rid[1,2],2],coord_rid[T_rid[1,3],2]))
for (ne in 1:dim(T_rid)[1])
{
    polygon(c(coord_rid[T_rid[ne,1],1],coord_rid[T_rid[ne,2],1],coord_rid[T_rid[ne,3],1]),c(coord_rid[T_rid[ne,1],2],coord_rid[T_rid[ne,2],2],coord_rid[T_rid[ne,3],2]))
}

#No alla edge matrix

e = NULL

#Ordine

order=1

basisobj = create.FEM.basis(coord_rid, e, T_rid, order)

#Creo l'oggetto fd

Rifiutifd = fd(numeric(basisobj$nbasis),basisobj)

#Risposta

data = matrix(0,nrow=length(TotC_rid),ncol=2)
data[,1] = 1:length(TotC_rid)
data[,2] = TotC_rid


#Applico il modello

lambda = 10^(-3)
Rifiuti = smooth.FEM.fd(data,Rifiutifd,lambda)

#Stima nei punti

fhat=Rifiuti$felsplobj$coef
#fhat
mean((fhat-TotC_rid)^2)

#Plot

plot.FEM(Rifiuti$felsplobj)

#MOLTO LENTO!!!!
#Cross-validazione

loglam = seq(1,5,by=0.5)
nloglam = length(loglam)
SEsaveCovar = numeric(nloglam)
np<-dim(coord_rid)[1]
for (ilam in 1:nloglam)
{
    print(c(ilam,loglam[ilam]))
    lami = 10^loglam[ilam]
    ressave = numeric(np)
    for (idrop in 1:np)
    {
        datai = data[-idrop,]
        Rifiuti_i   = smooth.FEM.fd(datai,Rifiutifd,lami)
        fhati = Rifiuti_i$felsplobj$coef
        RifiutiDataFiti = eval.FEM.fd(coord_rid[idrop,1], coord_rid[idrop,2], Rifiuti_i$felsplobj)
        ressave[idrop] = data[idrop,2] - RifiutiDataFiti
    }
    SEsaveCovar[ilam] = sqrt(mean(ressave^2))
}

print(rbind(loglam, SEsaveCovar))

plot(loglam,SEsaveCovar)
points(loglam,SEsaveCovar,type="l")
