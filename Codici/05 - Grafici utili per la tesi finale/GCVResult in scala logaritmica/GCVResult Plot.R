
load("01.1 - GCVResult.RData")

GCVMat<-log(GCVResult$GCVMatrix)
png(filename="LogGCVMatrix.png")
image(LogS,LogT,GCVMat,xlab="LogLambdaS",ylab="LogLambdaT",main=paste("LogGCVMatrix (Best (",round(GCVResult$Best[1],2),",",round(GCVResult$Best[2],2),"))",sep=""))
dev.off()

GCVMatZoom<-GCVMat[7:10,4:8]
LogSZoom<-LogS[7:10]
LogTZoom<-LogT[4:8]
png(filename="LogGCVMatrix Zoom.png")
image(LogSZoom,LogTZoom,GCVMatZoom,xlab="LogLambdaS",ylab="LogLambdaT",main=paste("LogGCVMatrix (Best (",round(GCVResult$Best[1],2),",",round(GCVResult$Best[2],2),"))",sep=""))
dev.off()
