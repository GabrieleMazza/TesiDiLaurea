#ANALISI DELLA FRAZIONE DI RACCOLTA DIFFERENZIATA

#Si ha un aumento in percentuale della raccolta differenziata?

Risposta<-read.table(file="Risposta.txt",header=T)

#Ora creo una matrice con i dati
MatData<-NULL

for(i in 1:580)
{
    MatData<-cbind(MatData,Risposta$CoeffDiff[Risposta$Codice==i])
}

#Traccio il matplot di questi valori
png(filename="RapportoRaccoltaDifferenziata.png")
matplot(1997:2011,MatData,type='l',xlab="Anno",ylab="Diff/Tot",main="Rapporto di Raccolta Differenziata")
dev.off()

#Isolo un comune e provo a farci una regressione logistica
#Senza covariate
j<-4
Logistic<-NULL
for (i in 2:dim(MatData)[1])
{
    if(MatData[i,j]>=MatData[i-1,j])
    {
        Logistic<-c(Logistic,1)
    } else
    {
        Logistic<-c(Logistic,0)
    }
}

Years<-1998:2011
fit <- glm(Logistic ~ Years, family='binomial')
summary(fit)
