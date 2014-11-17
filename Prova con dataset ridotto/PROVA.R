#PROVA CON DATASET RIDOTTO
source("Functions.R")

##################################
### RIDUZIONE DELLA FRONTIERA ####
##################################

#La frontiera è memorizzata come una unione di 130 poligoni (a causa delle isole della)
#laguna veneta)
#Prendo solo la regione, ed è formata da un poligono di 20000 vertici.. lo riduco in base
#alla tolleranza indicata
#pIù bassa è la tolleranza, più punti sono tenuti per il poligono di confine

#Punti di confine
load("Boundaries.RData")

#OGGETTI IN CUI METTERO' I PUNTI
xbound_rid<-NULL
ybound_rid<-NULL

#Scelgo un criterio per fissare solo alcuni punti.
#Inanzitutto per ora studio solo la label 1
#Ci sono 34690 punti

dist=0
xbound_rid<-c(xbound_rid,xbound[1])
ybound_rid<-c(ybound_rid,ybound[1])
tol<-0.005

for(i in 2:34690)
{
    if(labels[i]==1)
    {
        #Aggiorno la distanza
        dist<-dist+sqrt((xbound[i]-xbound[i-1])^2+(ybound[i]-ybound[i-1])^2)
        if(dist>tol)
        {
            xbound_rid<-c(xbound_rid,xbound[i])
            ybound_rid<-c(ybound_rid,ybound[i])
            dist<-0
        }
    }  
}
#NUMERO DI PUNTI CON CUI DESCRIVO LA FRONTIERA
length(xbound_rid)

#Controllo se ci sono intersezioni
#Intersect<-Intersections(xbound_rid,ybound_rid)
#Intersect
#Se sono segnalate intersezioni, allora devo intervenire, togliendo alcuni punti

plot(xbound_rid,ybound_rid,type='l',xlim=c(12.415,12.425),ylim=c(44.80,44.82))
identify(xbound_rid,ybound_rid,pos=T,plot=T)

#Casi analizzati
#TOLLERANZA 0.05
IDdelete<-c(182,185,186,192)
#TOLLERANZA 0.005
IDdelete<-c(1714,1715)

xbound_rid=xbound_rid[-IDdelete]
ybound_rid=ybound_rid[-IDdelete]

###########################################
###### RIDUZIONE DEI PUNTI DI COMUNI ######
###########################################

#Punti interni
load("Comuni.RData")
#Leggo la risposta (rifiuti prodotto).. anche questa va ridotta
Data<-read.table("Rifiuti.txt",header=T)

#Devo eliminare venezia e chioggia dal dataset, perchè sono isole e non sono nel
#poligono conservato. Ecco gli ID
#Chioggia--->391
#Venezia---->425

#Devo affiancare ai punti di comune gli id e la risposta

TotC<-NULL

for (i in 1:length(IDcom))
{
    control<-FALSE
    j<-1
    while (control==FALSE)
    {
        if(Data$IDans[j]==IDcom[i])
        {
            TotC<-c(TotC,Data$TotC[j])
            control<-TRUE
        }
        j<-j+1;
    }
}


#Ora devo trovare gli indici di riga corrispondenti agli ID da eliminare, che sono
#Chioggia--->391
#Venezia---->425

IDdelete<-c(391,425)
ROWdelete<-NULL

for (i in 1:length(IDcom))
{
    for (j in 1:length(IDdelete))
    {
        if (IDcom[i]==IDdelete[j])
        {
            ROWdelete<-c(ROWdelete,i)
        }
    }
}

xcom_rid<-xcom[-ROWdelete]
ycom_rid<-ycom[-ROWdelete]
IDcom_rid<-IDcom[-ROWdelete]
TotC_rid<-TotC[-ROWdelete]

#IMPORTANTE!!
#LA FRONTIERA RIDOTTA DESCRIVE UN POLIGONO CHE RACCHIUDE TUTTI I COMUNI?

library(SDMTools)
PolyPoints<-cbind(xbound_rid,ybound_rid)
if(sum(pnt.in.poly(cbind(xcom_rid,ycom_rid),PolyPoints)$pip)==length(xcom_rid))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}

#SE NO, SI ABBASSA LA TOLLERANZA...

save(file="Ridotto.RData",xbound_rid,ybound_rid,xcom_rid,ycom_rid,IDcom_rid,TotC_rid)

##############################################
##### TRIANGOLAZIONE E PROVA DEL MODELLO #####
##############################################

#Librerie
library(fda)
library(rgl)
library(RTriangle)
#Punti di confine e di comune
load("Ridotto.RData")
#Per le funzioni di modello
source("2013_SSR_AllFunctions.R")

ncom<-length(xcom_rid)
nbound<-length(xbound_rid)

#Creazione degli oggetti definitivi

x_rid<-c(xcom_rid,xbound_rid)
y_rid<-c(ycom_rid,ybound_rid)

#Devo creare la matrice di confine
Boundaries<-NULL
for(i in (ncom+1):(ncom+nbound-1))
{
    Boundaries<-rbind(Boundaries,c(i,i+1))
}
Boundaries<-rbind(Boundaries,c(ncom+nbound,ncom+1))
PolyPoints<-cbind(xbound_rid,ybound_rid)

#TRIANGOLAZIONE CON RTRIANGLE
#Prima devo creare l'oggetti pslg. Per poterlo fare mi serve
coord<-cbind(x_rid,y_rid)
#Oggetto pslg
pslg_obj<-pslg(coord,S=Boundaries)
#Creo la mesh
#Y dice di non aggiungere Steiner Points
#D dice di triangolare con Delaunay
mesh<-triangulate(pslg_obj,Y=TRUE,D=TRUE)
#Estrazione dei triangoli
T_rid<-mesh$T

#CONTROLLO: SONO STATI AGGIUNTI DEI PUNTI?
dim(mesh$P)[1]-dim(coord)[1]

#Plot
plot(x_rid,y_rid,type="n")
for (ne in 1:dim(T_rid)[1])
{
    polygon(c(x_rid[T_rid[ne,1]],x_rid[T_rid[ne,2]],x_rid[T_rid[ne,3]]),c(y_rid[T_rid[ne,1]],y_rid[T_rid[ne,2]],y_rid[T_rid[ne,3]]))
}
points(PolyPoints,type='l',col="red")


#No edge matrix
e = NULL

#Ordine
order=1
basisobj = create.FEM.basis(cbind(x_rid,y_rid), e, T_rid, order)

#Creo l'oggetto fd
Rifiutifd = fd(numeric(basisobj$nbasis),basisobj)

#Risposta
data = matrix(0,nrow=length(TotC_rid),ncol=2)
data[,1] = 1:dim(data)[1]
data[,2] = TotC_rid


#Applico il modello
lambda = 10^(3)
Rifiuti = smooth.FEM.fd(data,Rifiutifd,lambda)

#Stima nei punti
fhat=Rifiuti$felsplobj$coef
fhat

#Plot
plot.FEM(Rifiuti$felsplobj)