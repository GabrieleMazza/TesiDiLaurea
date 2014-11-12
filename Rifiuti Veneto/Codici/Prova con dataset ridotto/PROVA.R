#PROVA CON DATASET RIDOTTO

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
tol<-0.05

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


plot(xbound_rid,ybound_rid,type='l',xlim=c(12.4,12.8),ylim=c(45.4,45.6))
identify(xbound_rid,ybound_rid,pos=T,plot=T)


#Tolgo alcuni punti per semplificare..
IDdelete<-c(71,74,78,80,83,85,89,90,94,95,114,115)
xbound_rid=xbound_rid[-IDdelete]
ybound_rid=ybound_rid[-IDdelete]

###########################################
###### RIDUZIONE DEI PUNTI DI COMUNI ######
###########################################

#Punti interni
load("Comuni.RData")
#Devo eliminare venezia e chioggia dal dataset, perchè sono isole e non sono nel
#poligono conservato. Ecco gli ID
#Chioggia--->391
#Venezia---->425

#leggo la risposta (rifiuti prodotto).. anche questa va ridotta
Data<-read.table("Rifiuti.txt",header=T)

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
#Ariano del Polesine (che tolgo in quanto esce dal dominio) -->533
IDdelete<-c(391,425,533)
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
library(deldir)
#library(RTriangle)
#Punti di confine e di comune
load("Ridotto.RData")
#Per le funzioni di modello
source("2013_SSR_AllFunctions.R")
#Per le funzioni di appoggio
source("Functions.R")

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


#CON DELDIR
tr<-deldir(x_rid,y_rid)
T_rid<-triMat(tr)
dim(T_rid)[1]
#ELIMINO I TRIANGOLI CHE NON SONO ALL'INTERNO DELLA REGIONE
T_rid<-CleanTriangulation(x_rid, y_rid, T_rid, Boundaries)
dim(T_rid)[1]

#Plot
plot(x_rid,y_rid,type="n",xlim=c(12.4,12.8),ylim=c(45.4,45.6))
for (ne in 1:dim(T_rid)[1])
{
    polygon(c(x_rid[T_rid[ne,1]],x_rid[T_rid[ne,2]],x_rid[T_rid[ne,3]]),c(y_rid[T_rid[ne,1]],y_rid[T_rid[ne,2]],y_rid[T_rid[ne,3]]))
}


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

