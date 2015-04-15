# Creazione della frontiera in poligono unico usamdo una tecnica di smoothing
# La tecnica scelta è regression splines


# Funzioni di appoggio
source("Functions.R")
load("Territorio.RData")
library(fda)
library(KernSmooth)
library(geosphere)
library(fda)
library(RTriangle)
library(SDMTools)

##### CALCOLO DI TUTTI I BARICENTRI #####
xG<-NULL
yG<-NULL
for (i in 1:(dim(Triang)[1]))
{
    xG<-c(xG,(x[Triang[i,1]]+x[Triang[i,2]]+x[Triang[i,3]])/3)
    yG<-c(yG,(y[Triang[i,1]]+y[Triang[i,2]]+y[Triang[i,3]])/3)
}


##### DOWNLOAD DELLA FRONTIERA #####

# Scarico la struttura dei comuni
require(raster)
ProvinciaVZ =  subset(getData('GADM', country='ITA', level=3), NAME_2=="Venezia")
plot(ProvinciaVZ)

save(file="ProvinciaVZ.RData",ProvinciaVZ)
load(file="ProvinciaVZ.RData")

Nomi<-slot(ProvinciaVZ, 'data')[,c(1,10)]

# La frontiera è memorizzata come una unione di 131 poligoni (a causa delle isole della)
# laguna veneta)

##### CICLO DI ASSEGNAZIONE #####
tmp<- slot(ProvinciaVZ, 'polygons')
sub.tmp <- slot(tmp[[2]],'Polygons') 
order<-slot(sub.tmp[[17]],"coords")
labels<-numeric(length(xG))

for(i in 1:length(xG))
{
    xT<-xG[i]
    yT<-yG[i]
    
    done=FALSE
    for(j in 1:44)
    {
        #print(paste("j",j))
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

## ASSEGNAMENTI FINALI ##
Assegnamenti<-NULL
for(k in 1:(dim(Triang)[1]))
{
    if(labels[k]==0)
    {
        Assegnamenti<-rbind(Assegnamenti,c(k,"MANCANTE"))
    } else
    {
        Assegnamenti<-rbind(Assegnamenti,c(k,Nomi[labels[k],2]))
    }
}
# Ora hanno proprio i nomi dei comuni.
# Vado a modificare quelli mancanti

# STAMPO QUELLI NON IDENTIFICATI
for(k in 1:length(xG))
{
    if(labels[k]==0)
    {
        xlim1<-c(xG[k]-0.02,xG[k]+0.02)
        xlim2<-c(xG[k]-0.08,xG[k]+0.08)
        
        ylim1<-c(yG[k]-0.1,yG[k]+0.1)
        ylim2<-c(yG[k]-0.3,yG[k]+0.3)
        
        
        png(filename=paste("TR",k," 1.png",sep=""))
        plot(x,y,col="white",xlim=xlim1,ylim=ylim1,main=paste(k))
        for (ne in 1:dim(Triang)[1])
        {
            polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
        }
        points(xG[k],yG[k],pch=16,col="red")
        #Controllo se la triangolazione ha triangoli con solo punti di bordo
        dev.off()
        
        png(filename=paste("TR",k," 2.png",sep=""))
        plot(x,y,col="white",xlim=xlim2,ylim=ylim2,main=paste(k))
        for (ne in 1:dim(Triang)[1])
        {
            polygon(c(x[Triang[ne,1]],x[Triang[ne,2]],x[Triang[ne,3]]),c(y[Triang[ne,1]],y[Triang[ne,2]],y[Triang[ne,3]]))
        }
        #Controllo se la triangolazione ha triangoli con solo punti di bordo
        #Li coloro
        points(xG[k],yG[k],pch=16,col="red")
        dev.off()
    }
    
}

## 25 RESTANO FUORI. LI DEVO ASSEGNARE MAUALMENTE. CI SONO AD ESEMPIO PONTI.

# SONO NUMEROSISSIMI QUELLI CHE VANNO A VENEZIA (PONTI)
Assegnamenti[c(103,105,106,107,117,149,161,167,173),2]="Venezia"
Assegnamenti[c(14),2]="Cona"
Assegnamenti[c(19,21,34,36,54),2]="Chioggia"
Assegnamenti[c(56),2]="Campolongo Maggiore"
Assegnamenti[c(84),2]="Noale"
Assegnamenti[c(182),2]="Cavallino-Treporti"
Assegnamenti[c(183),2]="Cavallino-Treporti"
Assegnamenti[c(192),2]="Quarto d' Altino"
Assegnamenti[c(193),2]="Quarto d' Altino"
Assegnamenti[c(194),2]="Meolo"
Assegnamenti[c(198),2]="Jesolo"
Assegnamenti[c(239),2]="San Michele Al Tagliamento"
Assegnamenti[c(259),2]="Teglio Veneto"

ID<-as.numeric(Assegnamenti[Assegnamenti[,2]=="MANCANTE",1])
cbind(yG[ID],xG[ID])
