load("ICResult.RData")
png(filename="Intervalli di confidenza.png")
plot(c(1,1,1),ICResult$ApproxCIBeta[1,],type='l',main="Intervalli di confidenza",xlab="Beta1",ylab="IC")
points(c(1,1,1),ICResult$ApproxCIBeta[1,],pch=8)
abline(h=0)
abline(h=Beta,col="blue")
dev.off()

