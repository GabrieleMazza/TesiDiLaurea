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

DuplicatedTriangulation = function (Triangles)
{
    newvalue=1
    Subsets=0
    for(t in 2:dim(Triangles)[1])
    {
        #Devo ancora processare questo triangolo?
        process=TRUE
        for(j in 1:(t-1))
        {
            if(process)
            {
                if(Triangles[t,1]==Triangles[j,1] && Triangles[t,2]==Triangles[j,2] && Triangles[t,3]==Triangles[j,3])
                {
                    #Devo assegnare il nuovo elemento in t
                    if(Subsets[j]==0)
                    {
                        Subsets[j]=newvalue
                        Subsets<-rbind(Subsets,newvalue)
                        newvalue=newvalue+1
                    } else
                    {
                        Subsets<-rbind(Subsets,Subsets[j])
                    }
                }
            }
            
        }
        if(process)
        {
            #Se non ha uguali a lui, è zero
            Subsets<-rbind(Subsets,0)
            process=FALSE
        }    
    }
    return(Subsets)
}
    