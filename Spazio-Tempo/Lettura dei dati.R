#Importo gli oggetti che userò per l'analisi

##### LETTURA DEL DATASET #####

#Importo il dataset
Data<-read.table(file="DatasetTotale.txt",header=TRUE)

##### LETTURA DELLE COORDINATE DEL COMUNE #####

#Importo le coordinate di comune
Coord<-read.table(file="Coordinate.txt",header=TRUE)

##### LETTURA DELLE COVARIATE #####

#Importo le variabili delle covariate
Coord<-read.table(file="Covariate.txt",header=TRUE)