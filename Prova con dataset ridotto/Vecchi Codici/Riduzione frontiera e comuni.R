#Riduzione dei punti di frontiera del veneto

#Punti di confine
load("Boundaries.RData")
#Punti interni
load("Comuni.RData")

##### RIDUZIONE DELLA FRONTIERA #####

xbound_rid<-NULL
ybound_rid<-NULL

#Scelgo un criterio per fissare solo alcuni punti.
#Inanzitutto per ora studio solo la label 1
#Ci sono 34690 punti

dist=0
xbound_rid<-c(xbound_rid,xbound[1])
ybound_rid<-c(ybound_rid,ybound[1])
tol<-0.01

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
length(xbound_rid)


##### RIDUZIONE DEI COMUNI E RISPOSTA #####

#Devo eliminare venezia e chioggia dal dataset. Ecco gli ID
#Chioggia--->391
#Venezia---->425

Data<-read.table("Rifiuti.txt",header=T)
Covar<-read.table("Covariate.txt",header=T)

#Devo affiancare ai punti di comune gli id e la risposta

TotC<-NULL
Alb<-NULL
Comp<-NULL

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
for (i in 1:length(IDcom))
{
    control<-FALSE
    j<-1
    while (control==FALSE)
    {
        if(Data$IDans[j]==IDcom[i])
        {
            Alb<-c(Alb,Covar$Alb[j])
            Comp<-c(Comp,Covar$Comp[j])
            control<-TRUE
        }
        j<-j+1;
    }
}

#Ora devo trovare gli indici di riga corrispondenti agli ID da eliminare, che sono
#Chioggia--->391
#Venezia---->425

IDdelete<-c(391,425)

#Comuni
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

#Covariate

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

Alb_rid<-Alb[-ROWdelete]
Comp_rid<-Comp[-ROWdelete]

#Ci sono dei punti che si replicano?
library(deldir)
if(sum(duplicatedxy(c(xcom_rid,xbound_rid),c(ycom_rid,ybound_rid)))==0)
{
    print("Nessun punto duplicato tra frontiera e comuni")
} else
{
    print("Punti duplicati!!!")
}

library(SDMTools)
PolyPoints<-cbind(xbound_rid,ybound_rid)
if(sum(pnt.in.poly(cbind(xcom_rid,ycom_rid),PolyPoints)$pip)==length(xcom_rid))
{
    print("Tutti i comuni stanno dentro")
} else
{
    print("Esistono comuni esterni alla frontiera")
}

save(file="Ridotto.RData",xbound_rid,ybound_rid,xcom_rid,ycom_rid,IDcom_rid,TotC_rid)
save(file="Covar.RData",Alb_rid,Comp_rid)