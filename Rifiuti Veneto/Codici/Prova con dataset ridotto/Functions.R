# FUNZIONI DI APPOGGIO

CleanTriangulation = function (xPoints, yPoints, Triangles, Boundaries)
{
    library(SDMTools)
    #Devo processare tutti i triangoli, calcolarne il baricentro e vedere se è
    #un punto interno o no al contorno definito da boundaries 
    
    #Oggetto che conterrà i triangoli rimasti
    TNew<-NULL
    
    #Costruisco l'oggetto per la funzione
    PolyPoints<-NULL
    for(i in 1:dim(Boundaries)[1])
    {
        PolyPoints<-rbind(PolyPoints,c(xPoints[Boundaries[i,1]],yPoints[Boundaries[i,1]]))
    }
        
    for (t in 1:dim(Triangles)[1])
    {
        #Calcolo il baricentro
        xG=(xPoints[Triangles[t,1]]+xPoints[Triangles[t,2]]+xPoints[Triangles[t,3]])/3
        yG=(yPoints[Triangles[t,1]]+yPoints[Triangles[t,2]]+yPoints[Triangles[t,3]])/3
        
        #Il punto è interno al poligono definito da Boundaries? Come semiretta
        #scelgo quella che sale verso l'alto verticale
        
        out=pnt.in.poly(cbind(xG,yG),PolyPoints)
        # Il triangolo è interno se bool è risultato TRUE
        if (out$pip==1)
        {
            TNew<-rbind(TNew,Triangles[t,])
        }
    }
    return(TNew)
}