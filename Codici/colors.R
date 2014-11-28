val <- (1:950)/950
val2<-NULL
for(i in 1:950)
{
    val2<-c(val2,1)
}
val2[1:500]=val[451:950]
val2[501:950]=val[1:450]
val<-val2
from_col <- "red"
to_col <- "green"

val_col <- colorRampPalette(c(from_col,to_col))(length(val))

plot(c(0,1), c(0,1))
for(i in seq_along(val_col)){
    abline(h=val[i], col=val_col[i])
}
