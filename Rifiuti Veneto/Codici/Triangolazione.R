#Esempio di triangolazione con i punti del comune

#Punti di confine
load("Boundaries.RData")
#Punti interni
load("Comuni.RData")

library(RTriangle)

#Prima devo creare l'oggetti pslg. Per poterlo fare mi serve
coord<-cbind(xcom,ycom)
dim(coord)

#Oggetto pslg
P<-pslg(coord)
names(P)
#Traccio i punti
plot(P,pch=16)
#L'oggetto contiene tutti e 581 i punti

#Creo la mesh
mesh<-triangulate(P,Y=TRUE,D=TRUE)
#Y è un'importante parametro. Se non è impostato a TRUE, sono aggiunti dei punti di
#Steiner al bordo e così ci sono molti punti non voluti per la triangolazione
#Invece D controlla che la triangolazione sia di Delaunay
plot(mesh)
#Senza avere dei bordi, ovviamente, è stato analizzato l'inviluppo convesso
names(mesh)

#Estrazione dei triangoli
mesh$T
#Con riferimento ai 581 punti iniziali
