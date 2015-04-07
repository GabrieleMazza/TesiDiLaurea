#Ordina boxplot

load("Risultati.RData")
NEWLABEL<-NULL
NEWRMSE<-NULL
NEWBETA<-NULL

Order<-c("TPS","SOAP","ST-PDE")
LABEL[LABEL=="STR-PDE"]="ST-PDE"
LABEL<-factor(LABEL,Order)

png(filename="Confronto tra i metodi.png")
boxplot(RMSE ~ LABEL, main="Confronto tra metodi", xlab="Metodo",ylab="RMSE",ordered=F)
dev.off()

png(filename="Confronto tra Beta.png")
boxplot(BETA ~ LABEL, main="Confronto tra Beta", xlab="Metodo",ylab="Beta")
dev.off()
