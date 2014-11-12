#Applicazione al dataset ridotti

#Librerie
library(fda)
library(rgl)
library(deldir)
#library(RTriangle)
#Punti di confine e di comune
load("Ridotto.RData")
#Per le funzioni di modello
source("2013_SSR_AllFunctions.R")
#Per le funzioni di appoggio
source("Functions.R")



##### TRIANGOLAZIONE CON I PUNTI DI BORDO #####

ncom<-length(xcom_rid)
nbound<-length(xbound_rid)

#Creazione degli oggetti definitivi

x_rid<-c(xcom_rid,xbound_rid)
y_rid<-c(ycom_rid,ybound_rid)

#Devo creare la matrice di confine per pslg
Boundaries<-NULL
for(i in (ncom+1):(ncom+nbound-1))
{
    Boundaries<-rbind(Boundaries,c(i,i+1))
}
Boundaries<-rbind(Boundaries,c(ncom+nbound,ncom+1))


#CON RTRIANGLE
#Creo l'oggetto
#pslgobj<-pslg(cbind(x_rid,y_rid),S=Boundaries)
#Warning sul controllo degli NA, ma non c'è problema...
#Ora triangolazione
#mesh_rid<-triangulate(pslgobj, Y=TRUE,D=TRUE)
#plot(mesh_rid)
#T_rid<-mesh_rid$T
#dim(T_rid)

#CON DELDIR
tr<-deldir(x_rid,y_rid)
T_rid<-triMat(tr)
dim(T_rid)[1]
T_rid<-CleanTriangulation(x_rid, y_rid, T_rid, Boundaries)
dim(T_rid)[1]

#Plot
plot(x_rid,y_rid,type="n")
for (ne in 1:dim(T_rid)[1])
{
    polygon(c(x_rid[T_rid[ne,1]],x_rid[T_rid[ne,2]],x_rid[T_rid[ne,3]]),c(y_rid[T_rid[ne,1]],y_rid[T_rid[ne,2]],y_rid[T_rid[ne,3]]))
}
#Controllo dei triangoli duplicati 
if(sum(DuplicatedTriangulation(T_rid)==0))
{
    print("Nessun triangolo duplicato")
} else
{
    print("Triangoli duplicati!!!")
}
PolyPoints<-NULL
for(i in 1:dim(Boundaries)[1])
{
    PolyPoints<-rbind(PolyPoints,c(xPoints[Boundaries[i,1]],yPoints[Boundaries[i,1]]))
}
points(PolyPoints,type='l',col="red")


save(file="T_RidBoundaries.RData",T_rid)





##### TRIANGOLAZIONE SENZA I PUNTI DI BORDO #####

#Ora uso solo i punti di comune
#CON DELDIR
tr<-deldir(xcom_rid,ycom_rid)
T_rid<-triMat(tr)
dim(T_rid)[1]
#Triangolazione dell'inviluppo convesso

#CON RTRIANGLE
#pslgobj<-pslg(cbind(xcom_rid,ycom_rid))
#mesh_rid<-triangulate(pslgobj,Y=TRUE,D=TRUE)
#plot(mesh_rid)
#T_rid<-mesh_rid$T

#Plot
plot(xcom_rid,ycom_rid,type="n")
for (ne in 1:dim(T_rid)[1])
{
    polygon(c(xcom_rid[T_rid[ne,1]],xcom_rid[T_rid[ne,2]],xcom_rid[T_rid[ne,3]]),c(ycom_rid[T_rid[ne,1]],ycom_rid[T_rid[ne,2]],ycom_rid[T_rid[ne,3]]))
}

#Controllo dei triangoli duplicati 
if(sum(DuplicatedTriangulation(T_rid)==0))
{
    print("Nessun triangolo duplicato")
} else
{
    print("Triangoli duplicati!!!")
}

save(file="T_RidComuni.RData",T_rid)





##### APPLICAZIONE CON I PUNTI DI BORDO #####


load("T_RidBoundaries.RData")

#Traccio i punti
x_rid<-c(xcom_rid,xbound_rid)
y_rid<-c(ycom_rid,ybound_rid)

#Traccio la mesh

plot(x_rid,y_rid,type="n")
for (ne in 1:dim(T_rid)[1])
{
    polygon(c(x_rid[T_rid[ne,1]],x_rid[T_rid[ne,2]],x_rid[T_rid[ne,3]]),c(y_rid[T_rid[ne,1]],y_rid[T_rid[ne,2]],y_rid[T_rid[ne,3]]))
}

#No alla edge matrix

e = NULL

#Matrice disegno
load("Covar.RData")
desmat = matrix(1,nrow=length(Alb_rid),ncol=2)
desmat[,1] = Alb_rid
desmat[,2] = Comp_rid


#Ordine

order=1

basisobj = create.FEM.basis(cbind(x_rid,y_rid),  e, T_rid, order)

#Creo l'oggetto fd

Rifiutifd = fd(numeric(basisobj$nbasis),basisobj)

#Risposta

data = matrix(0,nrow=length(TotC_rid),ncol=2)
data[,1] = 1:dim(data)[1]
data[,2] = TotC_rid


#Applico il modello

lambda = 10^(5)
Rifiuti = smooth.FEM.fd.Covar(data,desmat,Rifiutifd,lambda)

#Stima nei punti

fhat=Rifiuti$felsplobj$coef
fhat

#Plot

plot.FEM(Rifiuti$felsplobj)


##### APPLICAZIONE SENZA I PUNTI DI BORDO #####

load("T_RidComuni.RData")

#Traccio la mesh

plot(xcom_rid,ycom_rid,type="n")
for (ne in 1:dim(T_rid)[1])
{
    polygon(c(xcom_rid[T_rid[ne,1]],xcom_rid[T_rid[ne,2]],xcom_rid[T_rid[ne,3]]),c(ycom_rid[T_rid[ne,1]],ycom_rid[T_rid[ne,2]],ycom_rid[T_rid[ne,3]]))
}

#No alla edge matrix

e = NULL

#Ordine

order=1

basisobj = create.FEM.basis(cbind(xcom_rid,ycom_rid), e, T_rid, order)

#Creo l'oggetto fd

Rifiutifd = fd(numeric(basisobj$nbasis),basisobj)

#Risposta

data = matrix(0,nrow=length(TotC_rid),ncol=2)
data[,1] = 1:dim(data)[1]
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
np<-dim(x_rid)[1]
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
        RifiutiDataFiti = eval.FEM.fd(x_rid[idrop], y_rid[idrop], Rifiuti_i$felsplobj)
        ressave[idrop] = data[idrop,2] - RifiutiDataFiti
    }
    SEsaveCovar[ilam] = sqrt(mean(ressave^2))
}

print(rbind(loglam, SEsaveCovar))

plot(loglam,SEsaveCovar)
points(loglam,SEsaveCovar,type="l")
