# Corso di Statistica Applicata 
# Docente: Prof. Piercesare Secchi
# Esercitatore: Ing. Simone Vantini


# 5 Giugno 2013
# An introduction to smoothing
# Laura Sangalli


#setwd("D:/Dati/politecnico/2009_scuolaSIS_FDA_Caserta")


# Parte del materiale che segue è basato sul libro 
# Ramsay, Hooker e Graves, "Functional Data Analysis with R and Matlab", Springer, 2009
# e parte dei codici utilizzati sono presi dagli scripts del pacchetto fda di R.


##########################################################################################
##########################################################################################

# PRIMA PARTE
# STIMA DEL DATO FUNZIONALE
# SMOOTHING

###noisycurvebis=truecurve$X0vera+rnorm(101,0,0.0025)
###
###plot(abscissa,Xobs0,xlab="t",ylab="observed data")
###plot(abscissa,noisycurvebis,col="blue")
###
###Dataset<-as.data.frame(cbind(Abscissa=abscissa,X0=noisycurvebis))
###write.table(Dataset,"noisycurvebis.txt", row.names=FALSE, col.names=TRUE, quote=FALSE)


noisycurve=read.table("noisycurvebis.txt",header=T)
Xobs0=noisycurve$X0
abscissa=noisycurve$Abscissa
NT=length(abscissa)


par (mfrow=c(1,1))
plot(abscissa,Xobs0,xlab="t",ylab="observed data")

truecurve=read.table("truecurve.txt",header=T)

points(abscissa,truecurve$X0vera,type="l")

#QUESTA E' LA CURVA VERA, ALLE OSSERVAZIONI E' STATO AGGIUNTO UN RUMORE
#LA VARIABILITA' E' PICCOLA RISPETTO ALLA CURVA

#SI PARLA DI DATI FUNZIONALI IN CASI COME QUESTO
#SE LA GRIGLIA FOSSE LASCA E L'ERRORE ALTO, SI FA ANALISI DI DATI LONGITUDINALI
#UN AMBITO DIVERSO

#STIME DELLE DERIVATE PRIME E SECONDE

rappincX1=(Xobs0[3:NT]-Xobs0[1:(NT-2)])/(abscissa[3:NT]-abscissa[1:(NT-2)])
rappincX2=((Xobs0[3:NT]-Xobs0[2:(NT-1)])/(abscissa[3:NT]-abscissa[2:(NT-1)])-(Xobs0[2:(NT-1)]-Xobs0[1:(NT-2)])/(abscissa[2:(NT-1)]-abscissa[1:(NT-2)]))*2/(abscissa[3:(NT)]-abscissa[1:(NT-2)])


par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")


######
# REGRESSION SPLINES
######

# ---> slides regression splines


library(fda)

nbasis=12
#ORDINE+1 DELLE SPLINES
m=3

basis = create.bspline.basis(c(0,1), nbasis, m)
basis
plot(basis)
basismat   = eval.basis(abscissa, basis)
basismat
lsfit(basismat, Xobs0, intercept=FALSE)$coef

#STIME LEAST SQUARES FIT
#SONO I COEFFICIENTI

par(mfrow=c(1,1))
Xsp0=basismat %*% lsfit(basismat, Xobs0, intercept=FALSE)$coef
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0 ,type="l",col="blue",lwd=2)
abline(v=basis$params)

#LA STIMA DELLA FUNZIONE E' QUELLA BLU
#STIMO POI LE DERIVATE DERIVANDO LE FUNZIONI DI BASE E TENENDO GLI STESSI COEFFICIENTI

# per ottenere derivata prima
basismat1   = eval.basis(abscissa, basis,Lfdobj=1)
Xsp1=basismat1 %*% lsfit(basismat, Xobs0, intercept=FALSE)$coef


# per ottenere derivata seconda
basismat2   = eval.basis(abscissa, basis,Lfdobj=2)
Xsp2=basismat2 %*% lsfit(basismat, Xobs0, intercept=FALSE)$coef


par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsp1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xsp2 ,type="l",col="blue",lwd=2)

#LA DERIVATA SECONDA HA UN ERRORE ALTISSIMO
#QUESTO L'ABBIAMO FATTO SU UNA UNITA' STATISTICA.. POI NE AVREMO n
#SE DEVO LAVORARE SU GRANDEZZE DERIVATE, APPENA INTRODUCO IL RUMORE QUESTO AUMENTA
#IL RAPPORTO INCREMENTALE

#QUINDI DEVO FARE SMOOTHING

############################
# codice alternativo
#USA TUTTE FUNZIONI DI fda

Xsp = smooth.basis(argvals=abscissa, y=Xobs0, basis)
Xsp0bis  = eval.fd(abscissa, Xsp$fd) #  the curve smoothing the data
Xsp1bis  = eval.fd(abscissa, Xsp$fd, Lfd=1)
Xsp2bis = eval.fd(abscissa, Xsp$fd, Lfd=2)
df  = Xsp$df   #  the degrees of freedom in the smoothing curve  ############## step =65 df=24 va bene
df


############################


par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0bis ,type="l",col="green",lwd=2)
points(abscissa,Xsp0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsp1bis ,type="l",col="green",lwd=2)
points(abscissa,Xsp1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xsp2bis ,type="l",col="green",lwd=2)
points(abscissa,Xsp2 ,type="l",col="blue",lwd=2)

#QUINDI USARE LE FUNZIONI DI FDA CORRISPONDE AI CALCOLI MANUALI

# altra possibilità per calcolare derivate: 
# derivata di spline è ancora spline con gli stessi nodi, con opportuno ordine e coefficienti direttamente calcolabili
# da quelli della spline originaria

#COSA SUCCEDE SE SI USANO TROPPE POCHE BASI? 7 INVECE CHE 12

nbasis=7

basisbis = create.bspline.basis(c(0,1), nbasis, m)
basismatbis   = eval.basis(abscissa, basisbis)

Xsp0bis=basismatbis %*% lsfit(basismatbis, Xobs0, intercept=FALSE)$coef

basismat1bis   = eval.basis(abscissa, basisbis,Lfdobj=1)
Xsp1bis=basismat1bis %*% lsfit(basismatbis, Xobs0, intercept=FALSE)$coef

basismat2bis   = eval.basis(abscissa, basisbis,Lfdobj=2)
Xsp2bis=basismat2bis %*% lsfit(basismatbis, Xobs0, intercept=FALSE)$coef


par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0bis ,type="l",col="green",lwd=2)
points(abscissa,Xsp0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsp1bis ,type="l",col="green",lwd=2)
points(abscissa,Xsp1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xsp2bis ,type="l",col="green",lwd=2)
points(abscissa,Xsp2 ,type="l",col="blue",lwd=2)

#IN VERDE VEDO LA STIMA CON SOLO 7 BASI
#SEMBRANO POCHI, SPAZIO TROPPO POVERO
#NON COLGO I PICCHI
#TAGLIA LE MONTAGNE E RIEMPIE LE VALLI


#SE CONSIDERO TROPPE FUNZIONI DI BASE:

# provare con
nbasis=30

basister = create.bspline.basis(c(0,1), nbasis, m)
basismatter   = eval.basis(abscissa, basister)

Xsp0ter=basismatter %*% lsfit(basismatter, Xobs0, intercept=FALSE)$coef

basismat1ter   = eval.basis(abscissa, basister,Lfdobj=1)
Xsp1ter=basismat1ter %*% lsfit(basismatter, Xobs0, intercept=FALSE)$coef

basismat2ter   = eval.basis(abscissa, basister,Lfdobj=2)
Xsp2ter=basismat2ter %*% lsfit(basismatter, Xobs0, intercept=FALSE)$coef



par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsp0ter ,type="l",col="red",lwd=2)
points(abscissa,Xsp0bis ,type="l",col="green",lwd=2)
points(abscissa,Xsp0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsp1ter ,type="l",col="red",lwd=2)
points(abscissa,Xsp1bis ,type="l",col="green",lwd=2)
points(abscissa,Xsp1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xsp2ter ,type="l",col="red",lwd=2)
points(abscissa,Xsp2bis ,type="l",col="green",lwd=2)
points(abscissa,Xsp2 ,type="l",col="blue",lwd=2)

#SOLITO BIAS-VARIANCE TRADEOFF
#STIAMO FACENDO OVERFITTING, SI VEDE BENE DALLA DERIVATA

#COME SI SCEGLIE IL NUMERO GIUSTO?

# generalized cross-validation

#GIA' CODIFICATO NEL PACCCHETTO

nbasis=6:30
gcv=numeric(length(nbasis))
for (i in 1:length(nbasis))
{basis = create.bspline.basis(c(0,1), nbasis[i], m)
gcv[i] = smooth.basis(abscissa, Xobs0, basis)$gcv
}
plot(nbasis,gcv)
nbasis[which.min(gcv)]

#VEDO IL GCV AL VARIARE DEL NUMERO DI BASI
#ESISTE UNA FORMA CHIUSA DI CROSS VALIDAZIONE
#SI HA UNA FORMA CHIUSA PER CALCOLARE UN INDICE DI CROSS VALIDAZIONE

#MA L'INDICE NON E' STATO PROGETTATO PER QUESTO AMBITO!
#NON FUNZIONA SEMPRE BENE

#PER SCEGLIERE IL NUMERO DI BASE E' MEGLIO GUARDARE LE STIME DELLE DERIVATE PRIME E SECONDE
#E VEDERE QUALE RACCONTA MEGLIO IN DATO

############################################################################
############################################################################

######
# SMOOTHING SPLINES
######

# ---> slides smoothing splines

#CIOE' CASO DI BASI CON PENALIZZAZIONE
#FAREMO UNA PENALIZZAZIONE DELLA DERIVATA TERZA

breaks=abscissa[((0:50)*2)+1]
#breaks=abscissa
basis = create.bspline.basis(breaks, norder=3)

#Qui si crea solo un oggetto, non la matrice che mi interessa
functionalPar = fdPar(fdobj=basis, Lfdobj=2, lambda=10e-8)  # functional parameter, da dare in pasto alla funzione smooth.basis
                                      # ha per argomenti: la base, l'ordine di derivata che si vuole penalizzare, 
                                      # il parametro di smoothing 

#DICO IN QUESTO COMANDO COME PENALIZZARE
#GLI DO GIA' LA BASE E IL LAMBDA
#SE VOGLIO POSSO PENALIZZARE ANCHE CON OPERATORI DIFFERENZIALI

#Da qui estrarrò la matrice di penalizzazione...
Xss=smooth.basis(abscissa, Xobs0, functionalPar)

Xss0 = eval.fd(abscissa, Xss$fd, Lfd=0)
Xss1 = eval.fd(abscissa, Xss$fd, Lfd=1)
Xss2 = eval.fd(abscissa, Xss$fd, Lfd=2)

df  = Xss$df   #  the degrees of freedom in the smoothing curve
df
gcv = Xss$gcv  #  the value of the gcv statistic
gcv

#OVVIAMENTE DOVRO' GUARDARLO FACENDO VARIARE LAMBDA

par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xss0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l",asp=1)
points(abscissa,Xss1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xss2 ,type="l",col="blue",lwd=2)

#ANCORA STO FACENDO OVERSMOOTHING?

#CON LAMBDA TROPPO GRANDE POTREI INDURRE TROPPA REGOLARIZZAZIONE

functionalParbis = fdPar(fdobj=basis, Lfdobj=3, lambda=10e-6)  # functional parameter, da dare in pasto alla funzione smooth.basis
                                      # ha per argomenti: la base, l'ordine di derivata che si vuole penalizzare, 
                                      # il parametro di smoothing 

Xssbis=smooth.basis(abscissa, Xobs0, functionalParbis)

Xss0bis = eval.fd(abscissa, Xssbis$fd, Lfd=0)
Xss1bis = eval.fd(abscissa, Xssbis$fd, Lfd=1)
Xss2bis = eval.fd(abscissa, Xssbis$fd, Lfd=2)

dfbis  = Xssbis$df   #  the degrees of freedom in the smoothing curve
dfbis
gcvbis = Xssbis$gcv  #  the value of the gcv statistic
gcvbis

#MOLTI MENO GRADI DI LIBERTA' E GCV AUMENTETO

#SE LAMBDA E' TROPPO PICCOLO FACCIO UNDERSMOOTHING E VADO AD INTERPOLARE ANCHE L'ERRORE

functionalParter = fdPar(fdobj=basis, Lfdobj=3, lambda=10e-11)  # functional parameter, da dare in pasto alla funzione smooth.basis
                                      # ha per argomenti: la base, l'ordine di derivata che si vuole penalizzare, 
                                      # il parametro di smoothing 

Xsster=smooth.basis(abscissa, Xobs0, functionalParter)

Xss0ter = eval.fd(abscissa, Xsster$fd, Lfd=0)
Xss1ter = eval.fd(abscissa, Xsster$fd, Lfd=1)
Xss2ter = eval.fd(abscissa, Xsster$fd, Lfd=2)

dfter  = Xsster$df   #  the degrees of freedom in the smoothing curve
dfter
gcvter = Xsster$gcv  #  the value of the gcv statistic
gcvter


par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xss0ter ,type="l",col="red",lwd=2)
points(abscissa,Xss0bis ,type="l",col="green",lwd=2)
points(abscissa,Xss0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xss1ter ,type="l",col="red",lwd=2)
points(abscissa,Xss1bis ,type="l",col="green",lwd=2)
points(abscissa,Xss1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xss2ter ,type="l",col="red",lwd=2)
points(abscissa,Xss2bis ,type="l",col="green",lwd=2)
points(abscissa,Xss2 ,type="l",col="blue",lwd=2)

#IN VERDE LAMBDA E' TROPPO ALTO
#SI VEDE BENE DALLA DERIVATA PRIMA

#LA ROSSA E' PICCOLO LAMBDA.
#NON STO INTRODUCENDO LA REGOLARIZZAZIONE E STO INTERPOLANDO ANCHE L'ERRORE

######
# REGRESSIONE POLINOMIALE LOCALE
######

#IN QUESTO CASO FITTIAMO UN POLINOMIO CON UN KERNEL (CE NE SONO TANTI, COME GAUSSIANO)
#INVECE CHE INTERPOLARE f CON CARATTERISTICHE DI REGOLARITA'
#NE FACCIO TANTE LOCALI (RETTA O PARABOLA)
#LA FACCIO SOLO SU UNA PICCOLA FINESTRINA DI DATI

#CONSIDERO TANTE FINESTRE IN CUI FACCIO QUESTA STIMA
#LO SI FA PESANDO I PUNTI CON UNA FUNZIONE DI PESO DETTA KERNEL
#PESO POI CON UN KERNEL

#POSSO TENERE LA LUNGHEZZA DELLA FINESTRA FISSA O POSSO APRIRLA O RESTRINGERLA IN BASE ALLA DENSITA' DEI PUNTI

#LA VARIABILITA' DELLE STIME DEI COEFFICIENTI DIPENDE DALLA VARIABILITA' DEI DATI
#A SECONDA DI QUANTI PUNTI CO SONO NELLA FINESTRA CAMBIA LA VARIABILITA' DEGLI STIMATORI

#SE LA VARIO MA MANTENTENDO LE NUMEROSITA' DEI COEFF NELLE FINESTRE, MANTERRO' PIU' O MENO
#COSTANTE LA VARIABILITA'

#L'EQUIVALENTE DI LAMBDA E' QUANTO E' GRANDE O PICCOLA LA FINESTRA DI STIMA

library(KernSmooth)

m=5           # ordine del polinomio
degree=m-1    # grado del polinomio

bw=0.05

Xsm0 <- locpoly(abscissa,Xobs0,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm0 <- locpoly(abscissa,Xobs0,degree=degree,bandwidth=bw,gridsize=150, range.x=range(abscissa))
Xsm0 <- Xsm0$y


par(mfrow=c(1,1))
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsm0 ,type="l",col="blue")


Xsm1 <- locpoly(abscissa,Xobs0,drv=1,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm1 <- Xsm1$y

Xsm2 <- locpoly(abscissa,Xobs0,drv=2,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm2 <- Xsm2$y

par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsm0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsm1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xsm2 ,type="l",col="blue",lwd=2)

#QUESTA E' LA STIMA OTTENUTA CON POLINOMI LOCALI. POSSO ANCHE ORA AVERE LE STIME DELLE DERIVATE


bw=0.15

Xsm0bis <- locpoly(abscissa,Xobs0,drv=0,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm0bis <- Xsm0bis$y

Xsm1bis <- locpoly(abscissa,Xobs0,drv=1,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm1bis <- Xsm1bis$y

Xsm2bis <- locpoly(abscissa,Xobs0,drv=2,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm2bis <- Xsm2bis$y


par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsm0bis ,type="l",col="green",lwd=2)
points(abscissa,Xsm0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsm1bis ,type="l",col="green",lwd=2)
points(abscissa,Xsm1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xsm2bis ,type="l",col="green",lwd=2)
points(abscissa,Xsm2 ,type="l",col="blue",lwd=2)

# oversmoothing

#SE SCELGO UN bw TROPPO GRANDE FACCIO TROPPO SMOOTHING
#SE LO SCELGO TROPPO PICCOLO HO UNDERSMOOTHING

# riprovare con 
bw=0.015 # undersmoothing
Xsm0ter <- locpoly(abscissa,Xobs0,drv=0,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm0ter <- Xsm0ter$y

Xsm1ter <- locpoly(abscissa,Xobs0,drv=1,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm1ter <- Xsm1ter$y

Xsm2ter <- locpoly(abscissa,Xobs0,drv=2,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm2ter <- Xsm2ter$y

par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(abscissa,Xobs0,xlab="t",ylab="observed data")
points(abscissa,Xsm0ter ,type="l",col="red",lwd=2)
points(abscissa,Xsm0bis ,type="l",col="green",lwd=2)
points(abscissa,Xsm0 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="t",ylab="first differences x",type="l")
points(abscissa,Xsm1ter ,type="l",col="red",lwd=2)
points(abscissa,Xsm1bis ,type="l",col="green",lwd=2)
points(abscissa,Xsm1 ,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX2,xlab="t",ylab="second differences x",type="l")
points(abscissa,Xsm2ter ,type="l",col="red",lwd=2)
points(abscissa,Xsm2bis ,type="l",col="green",lwd=2)
points(abscissa,Xsm2 ,type="l",col="blue",lwd=2)

# aiuta molto guardare stima derivata vs differenze finite, nella scelta del parametro di smoothing


######
# METODI ADATTIVI
######

# Free-knot regression splines
# ---> splides

#UNA ULTERIORE POSSIBILITA'
#UN NUMERO LIMITATO DI BASI MA SCELTE IN MODO ADATTIVO AI DATI

#CON LE SPLINES SONO DETTE DI NODI LIBERI, LASCERO' CHE I DATI SCELGANO DOVE PIAZZARLI

#OPPURE C'E' UNA BASE DATA DALLE COMPONENTI PRINCIPALI COME SCELTA ADATTATIVA

#DOVENDO FARLO PER n UNITA' STATISTICHE, POI AVRO' n FUNZIONI
#OGNUNA PUO' ESSERE OTTENUTA CON UNA TECNICA DIVERSA
#NON E' DETTO CHE USO SEMPRE LO STESO NUMERO DI BASI

#POI SI PUO' FARE DI TUTTO.. INTRODURRE SULLA PCA LA REGOLARIZZAZIONE...

#LA SPLINE CON NODI LIBERI E' FATTA CON REGOLARIZZAZIONE BASATA SUL NUMERO DI BASI CHE
#VOGLIO UTILIZZARE
#LI SCELGO IN MODO ADATTATIVO CON I DATI
#SI USANO ALGORITMI ITERATIVI PER MINIMIZZARE
#QUELLO DI ZHOU-SHEN AIUTA, INIZIA CON NODI EQUISPAZIATI E LI AGGIORNA IN CONTINUAZIONE

#VEDI ANCHE LE SLIDES


######
# SMOOTHING PER CURVE positive
######

#PER FARE SMOOTHING DI CURVE POSITIVE LE SI SCRIVE COSI:

# y_j = exp(w(t_j)) + e_j

#E w SI STIMA
#MA HO SOLO FORME NUMERICHE DI CALCOLO

# x(t) = exp(w(t)) 

# La funzione w(t) è unconstrained
# La funzione x(t) è monotona crescente

# w(t) viene modellizata via espansione in basi 
# coefficienti dell'espansione in basi vengono stimati con tecniche numeriche

#FUNZIONE CHE FA TUTTO
# smooth.pos


######
# SMOOTHING PER CURVE MONOTONE
######


# Esempio: Berkeley female growth data
 
#NEL PACCHETTO fda
growth

# Se si procede non considerando che le curve di crescita sono monotone crescenti

age     = growth$age
heightbasis12 = create.bspline.basis(c(1,18), 12, 6)
basismat   = eval.basis(growth$age, heightbasis12)
heightmat  = growth$hgtf
heightcoef = lsfit(basismat, heightmat, intercept=FALSE)$coef


height=basismat %*% lsfit(basismat, heightmat, intercept=FALSE)$coef

basismat1   = eval.basis(growth$age, heightbasis12,Lfdobj=1)
heightvelocity=basismat1 %*% lsfit(basismat, heightmat, intercept=FALSE)$coef

basismat2   = eval.basis(growth$age, heightbasis12,Lfdobj=2)
heightacceleration=basismat2 %*% lsfit(basismat, heightmat, intercept=FALSE)$coef


par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
matplot(age,height,type="l" )
matplot(age,heightvelocity,type="l" )
abline(h=0)
matplot(age,heightacceleration,type="l")


par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
matplot(age,height,type="l" )
matplot(age[-c(1,2,3,31)],heightvelocity[-c(1,2,3,31),],type="l" )
abline(h=0)
matplot(age[-c(1,2,3,31)],heightacceleration[-c(1,2,3,31),],type="l")

#C'E' UN PROBLEMA.. OCCORRE CAMBIARE MODELLO
# velocità negativa: insensato... significherebbe che ad un certo punto le bambine cominciano a decrescere...

# Modello per curve monotone

# y_j = b_0 + b_1 * x(t_j) + e_j

# x(t) = int_(t_0)^t exp(w(u)) du

#COSI' SI MODELLIZZA
#w E' CIO' CHE SI RACCONTA CON ESPANSIONE DI BASI
#L'INTEGRALE E' SEMBRE MONOTONO CRESCENTE

# La funzione w(t) è unconstrained
# La funzione x(t) è monotona crescente
# b_1>0 per funz monotone crescenti
# b_1<0 per funz monotone decrescenti
# b_0 è il valore della funz in t_0

# w(t) viene modellizata via espansione in basi 
# b_0, b_1 e w(t) vengono stimate con tecniche numeriche

#ANCORA PERDO LA CAPACITA' DI FARE CALCOLI ANALITICI

nage    = length(age)
ageRng  = range(age)
nfine   = 101
agefine = seq(ageRng[1], ageRng[2], length=nfine)

# consideriamo solo le prime 10 bambine

hgtf   = growth$hgtf[,1:10]
ncasef = dim(hgtf)[2]

#COSI' VELOCIZZIAMO I CALCOLI NUMERICI

# creiamo base Bspline di ordine 6 con nodi ad ogni età in cui abbiamo una misura 

norder = 6
nbasis = nage - 2 + norder 
wbasis = create.bspline.basis(ageRng, nbasis, norder, age)

# Costruiamo parametro funzionale con penalizzazione di derivata terza
Lfdobj    = 3          
lambda    = 10^(-0.5)  
cvecf     = matrix(0, nbasis, ncasef) # serve per valori iniziale tecniche numeriche
Wfd0      = fd(cvecf, wbasis)
growfdPar = fdPar(Wfd0, Lfdobj, lambda)

# smoothing monotono

growthMon = smooth.monotone(age, hgtf, growfdPar)

#QUESTA E' LA FUNZIONE CHE SI USA ED E' MOLTO LENTO

Wfd        = growthMon$Wfd
betaf      = growthMon$beta
hgtfhatfd  = growthMon$yhatfd

velocfdUN     = deriv.fd(hgtfhatfd, 1)
velocmeanfdUN = mean(velocfdUN)

accelfdUN     = deriv.fd(hgtfhatfd, 2)
accelmeanfdUN = mean(accelfdUN)


par (mfrow=c(2,2),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)
plot(hgtfhatfd, xlim=c(1,18), lty=1, lwd=2,
     cex=2, xlab="Age", ylab="Growth (cm)")
plot(velocfdUN, xlim=c(1,18),  lty=1, lwd=2,
     cex=2, xlab="Age", ylab="Velocity (cm/yr)")
plot(accelfdUN, xlim=c(1,18), ylim=c(-4,3), lty=1, lwd=2,
     cex=2, xlab="Age", ylab="Acceleration (cm/yr/yr)")

#ADESSO NON VA PIU' SOTTO ZERO
#QUESTO GRAZIED ALLA TIPOLOGIA DI MODELLO SCELTO




######
# ESTENSIONE AL CASO DI CURVE IN PIU' DIMENSIONI
######

#SMOOTHING DI SUPERFICI O DI CURVE SPAZIALI

noisycurve3D=read.table("noisycurve3D.txt",header=T)
Xobs0=noisycurve3D$X0
Yobs0=noisycurve3D$Y0
Zobs0=noisycurve3D$Z0
obs0=rbind(Xobs0,Yobs0,Zobs0)
abscissa=noisycurve3D$Abscissa
NT=length(abscissa)


# consideriamo il caso più semplice in cui si ha a disposizione il vettore dei tempi/posizioni delle osservazioni


truecurve3D=read.table("truecurve3D.txt",header=T)
Xtrue0=truecurve3D$X0
Ytrue0=truecurve3D$Y0
Ztrue0=truecurve3D$Z0
true0=rbind(Xtrue0,Ytrue0,Ztrue0)
abscissa=noisycurve3D$Abscissa




library(misc3d)
library(rgl)

open3d()
lines3d(t(true0[1,]),t(true0[2,]),t(true0[3,]),xlab="",ylab="",zlab="",size=3,axes=F)
points3d(t(obs0[1,]),t(obs0[2,]),t(obs0[3,]),xlab="",ylab="",zlab="",size=2,axes=F,pch=19,cex=2)
box3d()



rappincX1=(Xobs0[3:NT]-Xobs0[1:(NT-2)])/(abscissa[3:NT]-abscissa[1:(NT-2)])
rappincY1=(Yobs0[3:NT]-Yobs0[1:(NT-2)])/(abscissa[3:NT]-abscissa[1:(NT-2)])
rappincZ1=(Zobs0[3:NT]-Zobs0[1:(NT-2)])/(abscissa[3:NT]-abscissa[1:(NT-2)])

rappincX2=((Xobs0[3:NT]-Xobs0[2:(NT-1)])/(abscissa[3:NT]-abscissa[2:(NT-1)])-(Xobs0[2:(NT-1)]-Xobs0[1:(NT-2)])/(abscissa[2:(NT-1)]-abscissa[1:(NT-2)]))*2/(abscissa[3:(NT)]-abscissa[1:(NT-2)])
rappincY2=((Yobs0[3:NT]-Yobs0[2:(NT-1)])/(abscissa[3:NT]-abscissa[2:(NT-1)])-(Yobs0[2:(NT-1)]-Yobs0[1:(NT-2)])/(abscissa[2:(NT-1)]-abscissa[1:(NT-2)]))*2/(abscissa[3:(NT)]-abscissa[1:(NT-2)])
rappincZ2=((Zobs0[3:NT]-Zobs0[2:(NT-1)])/(abscissa[3:NT]-abscissa[2:(NT-1)])-(Zobs0[2:(NT-1)]-Zobs0[1:(NT-2)])/(abscissa[2:(NT-1)]-abscissa[1:(NT-2)]))*2/(abscissa[3:(NT)]-abscissa[1:(NT-2)])




par (mfrow=c(3,3),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)

plot(abscissa,obs0[1,],xlab=expression(tilde(s)),ylab="observed data x",cex=0.1,asp=1)
plot(abscissa,obs0[2,],xlab=expression(tilde(s)),ylab="observed data y",cex=0.1,asp=1)
plot(abscissa,obs0[3,],xlab=expression(tilde(s)),ylab="observed data z",cex=0.1,asp=1)
plot(abscissa[2:(NT-1)],rappincX1,xlab=expression(tilde(s)),ylab="first differences x",type="l",asp=1)
plot(abscissa[2:(NT-1)],rappincY1,xlab=expression(tilde(s)),ylab="first differences y",type="l",asp=1)
plot(abscissa[2:(NT-1)],rappincZ1,xlab=expression(tilde(s)),ylab="first differences z",type="l",asp=1)
plot(abscissa[2:(NT-1)],rappincX2,xlab=expression(tilde(s)),ylab="second differences x",type="l")
plot(abscissa[2:(NT-1)],rappincY2,xlab=expression(tilde(s)),ylab="second differences y",type="l")
plot(abscissa[2:(NT-1)],rappincZ2,xlab=expression(tilde(s)),ylab="second differences z",type="l")




bw=0.05

Xsm0 <- locpoly(abscissa,Xobs0,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm0 <- Xsm0$y

Xsm1 <- locpoly(abscissa,Xobs0,drv=1,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm1 <- Xsm1$y

Xsm2 <- locpoly(abscissa,Xobs0,drv=2,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Xsm2 <- Xsm2$y


Ysm0 <- locpoly(abscissa,Yobs0,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Ysm0 <- Ysm0$y

Ysm1 <- locpoly(abscissa,Yobs0,drv=1,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Ysm1 <- Ysm1$y

Ysm2 <- locpoly(abscissa,Yobs0,drv=2,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Ysm2 <- Ysm2$y


Zsm0 <- locpoly(abscissa,Zobs0,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Zsm0 <- Zsm0$y

Zsm1 <- locpoly(abscissa,Zobs0,drv=1,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Zsm1 <- Zsm1$y

Zsm2 <- locpoly(abscissa,Zobs0,drv=2,degree=degree,bandwidth=bw,gridsize=length(abscissa), range.x=range(abscissa))
Zsm2 <- Zsm2$y



par (mfrow=c(3,3),mar=c(6,5,2,1),mex=0.6, mgp=c(2.2,0.7,0),pty="m", font.main=1,font.lab=1, font.axis=1,cex.lab=1.3,cex.axis=1)

plot(abscissa,obs0[1,],xlab="s",ylab="x",cex=0.1,asp=1,xlim=c(0,1))
points(abscissa,Xsm0,type="l",col="blue",lwd=2)
plot(abscissa,obs0[2,],xlab="s",ylab="y",cex=0.1,asp=1,xlim=c(0,1))
points(abscissa,Ysm0,type="l",col="blue",lwd=2)
plot(abscissa,obs0[3,],xlab="s",ylab="z",cex=0.1,asp=1,xlim=c(0,1))
points(abscissa,Zsm0,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincX1,xlab="s",ylab="x'",type="l",ylim=c(-0.5,0.5),xlim=c(0,1))
points(abscissa,Xsm1,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincY1,xlab="s",ylab="y'",type="l",ylim=c(-0.5,0.5),xlim=c(0,1))
points(abscissa,Ysm1,type="l",col="blue",lwd=2)
plot(abscissa[2:(NT-1)],rappincZ1,xlab="s",ylab="z'",type="l",ylim=c(-0.5,0.5),xlim=c(0,1))
points(abscissa,Zsm1,type="l",col="blue",lwd=2)
#plot(abscissa[2:(NT-1)],rappincX2,xlab="s",ylab="x''",type="l")
#points(abscissa,Xsm2,type="l",col="blue",lwd=2)
#plot(abscissa[2:(NT-1)],rappincY2,xlab="s",ylab="y''",type="l")
#points(abscissa,Ysm2,type="l",col="blue",lwd=2)
#plot(abscissa[2:(NT-1)],rappincZ2,xlab="s",ylab="z''",type="l")
#points(abscissa,Zsm2,type="l",col="blue",lwd=2)


open3d()
lines3d(t(true0[1,]),t(true0[2,]),t(true0[3,]),xlab="",ylab="",zlab="",size=3,axes=F)
points3d(t(obs0[1,]),t(obs0[2,]),t(obs0[3,]),size=2,pch=19,cex=2)
lines3d(t(Xsm0),t(Ysm0),t(Zsm0),size=3,col="blue")
box3d()
