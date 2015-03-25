#MATPLOT DATI

load("Territorio.RData")

# Leggo i dati
CoordinateCovariate<-read.table("CoordinateCovariate.txt",header=T)
Risposta<-read.table("Risposta.txt",header=T)

# Punti di tempo
TimePoints<-1997:2011
InternalPoints<-!is.na(Codici)
nint<-sum(InternalPoints)

# Basi
TimeBasisObj<-Create.Bspline.Time.Basis(TimePoints,TimeOrder=4,DerivativeOrder=2,PlotIt=F)
SpaceBasisObj<-Create.FEM.Space.Basis(cbind(x,y),Triang,InternalPoints,1)

# Matrici di Dati
DataMatrix<-NULL
for(i in 1:length(Codici[1:nint]))
{
    DataTmp<-numeric(length(TimePoints))
    for(j in 1:length(TimePoints))
    {
        DataTmp[j]=Risposta$TotalePC[(Risposta$Codice==Codici[i]) & (Risposta$Anno==TimePoints[j])]
    }
    DataMatrix<-rbind(DataMatrix,DataTmp)
}

matplot(t(DataMatrix),type='l',col="black",lwd=1,xlab="Anni",ylab="Rifiuti",main="Produzione rifiuti",axes=F)
axis(side=2,at=c(0,500,1000,1500,2000),labels=c(0,500,1000,1500,2000))
axis(side=1,at=1:ncol(DataMatrix),labels=TimePoints)
