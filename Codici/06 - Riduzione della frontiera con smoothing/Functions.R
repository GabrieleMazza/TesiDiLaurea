# FUNZIONI DI APPOGGIO

SquareEuclideanDistance = function (p1,p2)
{
    dist=sqrt((p1[1]+p2[1])^2+(p1[1]+p2[1])^2)
    return(dist)
}

Intersections = function(x,y)
{
    # Studio se un poligono è semplice o complesso
    # Nel caso in cui sia complesso, restituisce gli indici delle coppie di 
    Intersect<-NULL
    if(length(x)!=length(y))
    {
        stop('Lengths of vectors is not the same')
    }
    if(length(x)<=3)
    {
        return(Intersect)
    }
    for(i in 1:(length(x)-2))
    {
        if(i==1)
        {
            #Non vado fino all'ultimo segmento..
            maximumj=length(x)-1
        } else
        {
            maximumj=length(x)
        }
        for(j in (i+2):maximumj)
        {
            #Inanzitutto cerco il punto di arrivo del lato che parte da j
            ind=j+1
            if(j==length(x))
            {
                ind=1
            }
            #### CASO DI SEGMENTI ENTRAMBI VERTICALI ####
            if((x[i]==x[i+1])&&(x[j]==x[ind]))
            {
                #Cioè se entrambi i segmenti sono verticali, si hanno due segmenti
                #Paralleli
                #Controllo che si intersecano
                if(x[i]==x[j])
                {
                    #Devono avere la stessa x
                    #E controllo se le y si intersecano
                    if(y[i]>y[i+1])
                    {
                        yimax=y[i]
                        yimin=y[i+1]
                    } else
                    {
                        yimax=y[i+1]
                        yimin=y[i]
                    }
                    if(y[j]>y[ind])
                    {
                        yjmax=y[j]
                        yjmin=y[ind]
                    } else
                    {
                        yjmax=y[ind]
                        yjmin=y[j]
                    }
                    #Si intersecano se, scelto il più alto, si ha
                    if (yimax>=yjmax)
                    {
                        if(yjmax>yimin)
                        {
                            Intersect<-rbind(Intersect,c(i,i+1,j,ind,"Segmenti sovrapposti"))
                        } else
                        {
                            if(yjmax==yimin)
                            {
                                Intersect<-rbind(Intersect,c(i,i+1,j,ind,"Intersezione"))
                            }
                        }
                    } else
                    {
                        if(yimax>yjmin)
                        {
                            Intersect<-rbind(Intersect,c(i,i+1,j,ind,"Segmenti sovrapposti"))
                        } else
                        {
                            if(yimax==yjmin)
                            {
                                Intersect<-rbind(Intersect,c(i,i+1,j,ind,"Intersezione"))
                            }
                        }
                    }
                }
            } else
            {
                #### CASO DI SEGMENTO CHE PARTE DA i VERTICALE ####
                if(x[i]==x[i+1])
                {
                    #Controllo che il segmento j lo intersechi
                    if(x[j]>x[ind])
                    {
                        xjmax=x[j]
                        xjmin=x[ind]
                    } else
                    {
                        xjmax=x[ind]
                        xjmin=x[j]
                    } 
                    if((x[i]>xjmin)&&(x[i]<xjmax))
                    {
                        #Ok, è compreso. Contollo che ci siano intersezioni
                        #Ricavo la retta che passa per il segmento j
                        mj=(y[j]-y[ind])/(x[j]-x[ind])
                        qj=y[j]-mj*x[j]
                        yint=mj*x[i]+qj
                        #Questa yint appartiene al segmento?
                        if(y[i]>y[i+1])
                        {
                            yimax=y[i]
                            yimin=y[i+1]
                        } else
                        {
                            yimax=y[i+1]
                            yimin=y[i]
                        }
                        if((yint>=yimin)&&(yint<=yimax))
                        {
                            #Se vale con il maggiore, intersezione
                            if((yint>yimin)&&(yint<yimax))
                            {
                                Intersect<-rbind(Intersect,c(i,i+1,j,ind,"Intersezione"))
                            } else
                            {
                                Intersect<-rbind(Intersect,c(i,i+1,j,ind,"Intersezione ad un estremo"))
                            }
                        }
                    }
                } else
                {
                    #### CASO DI SEGMENTO CHE PARTE DA j VERTICALE ####
                    if(x[j]==x[ind])
                    {
                        #Controllo che il segmento j lo intersechi
                        if(x[i]>x[i+1])
                        {
                            ximax=x[i]
                            ximin=x[i+1]
                        } else
                        {
                            ximax=x[i+1]
                            ximin=x[i]
                        } 
                        if((x[j]>ximin)&&(x[j]<ximax))
                        {
                            #Ok, è compreso. Contollo che ci siano intersezioni
                            #Ricavo la retta che passa per il segmento i
                            mi=(y[i]-y[i+1])/(x[i]-x[i+1])
                            qi=y[i]-mi*x[i]
                            yint=mi*x[j]+qi
                            #Questa yint appartiene al segmento?
                            if(y[j]>y[ind])
                            {
                                yjmax=y[j]
                                yjmin=y[ind]
                            } else
                            {
                                yjmax=y[ind]
                                yjmin=y[j]
                            }
                            if((yint>=yjmin)&&(yint<=yjmax))
                            {
                                #Se vale con il maggiore, intersezione
                                if((yint>yjmin)&&(yint<yjmax))
                                {
                                    Intersect<-rbind(Intersect,c(i,i+1,j,ind,"Intersezione"))
                                } else
                                {
                                    Intersect<-rbind(Intersect,c(i,i+1,j,ind,"Intersezione ad un estremo"))
                                }
                            }
                        }
                    } else
                    {
                        #Controllo se si intersecano i lati che partono dal punto i e dal punto j
                        #Ricavo la retta che passa per il lato che parte da i
                        mi=(y[i]-y[i+1])/(x[i]-x[i+1])
                        qi=y[i]-mi*x[i]
                        #Ricavo la retta del lato che parte da j
                        mj=(y[j]-y[ind])/(x[j]-x[ind])
                        qj=y[j]-mj*x[j]
                        #Cerco l'intersezione
                        if(mi!=mj)
                        {
                            #Delta
                            Delta=mi-mj
                            xint=-(qi-qj)/Delta
                            #yint=-(mi*qj-mj*qi)/Delta
                            #Basta la x, controllo che sia compresa in entrambi i segmenti
                            if(x[i]>x[i+1])
                            {
                                ximax=x[i]
                                ximin=x[i+1]
                            } else
                            {
                                ximax=x[i+1]
                                ximin=x[i]
                            }
                            if(x[j]>x[ind])
                            {
                                xjmax=x[j]
                                xjmin=x[ind]
                            } else
                            {
                                xjmax=x[ind]
                                xjmin=x[j]
                            }
                            if((xint>=ximin)&&(xint<=ximax)&&(xint>=xjmin)&&(xint<=xjmax))
                            {
                                #Se vale con il maggiore l'intersezione è interna
                                if ((xint>ximin)&&(xint<ximax)&&(xint>xjmin)&&(xint<xjmax))
                                {
                                    Intersect<-rbind(Intersect,c(i,i+1,j,ind,"Intersezione"))
                                } else
                                {
                                    Intersect<-rbind(Intersect,c(i,i+1,j,ind,"Intersezione ad un estremo"))
                                }
                            }
                        } else
                        {
                            #Allora sono paralleli. Ma si intersecano?
                            #Escludo il caso in cui siano paralleli ma distanti
                            if(qi==qj)
                            {
                                #E controllo se le y si intersecano
                                if(x[i]>x[i+1])
                                {
                                    ximax=x[i]
                                    ximin=x[i+1]
                                } else
                                {
                                    ximax=x[i+1]
                                    ximin=x[i]
                                }
                                if(x[j]>x[ind])
                                {
                                    xjmax=x[j]
                                    xjmin=x[ind]
                                } else
                                {
                                    xjmax=x[ind]
                                    xjmin=x[j]
                                }
                                #Si intersecano se, scelto il più alto, si ha
                                if (ximax>=xjmax)
                                {
                                    if(xjmax>ximin)
                                    {
                                        Intersect<-rbind(Intersect,c(i,i+1,j,ind,"Segmenti sovrapposti"))
                                    } else
                                    {
                                        if(xjmax==ximin)
                                        {
                                            Intersect<-rbind(Intersect,c(i,i+1,j,ind,"Intersezione"))
                                        }
                                    }
                                } else
                                {
                                    if(ximax>xjmin)
                                    {
                                        Intersect<-rbind(Intersect,c(i,i+1,j,ind,"Segmenti sovrapposti"))
                                    } else
                                    {
                                        if(ximax==xjmin)
                                        {
                                            Intersect<-rbind(Intersect,c(i,i+1,j,ind,"Intersezione"))
                                        }
                                    }
                                }
                            }
                        }
                    }
                }                
            } 
        }
    }
    return(Intersect)
}

BorderTriangles = function (Triangulation,Boundaries)
{
    firstindex=Boundaries[1,1]
    ReturnTriangles<-NULL
    for(i in 1:dim(Triangulation)[1])
    {
        count<-0
        for(j in 1:3)
        {
            if(Triangulation[i,j]>=firstindex)
            {
                count<-count+1
            }
        }
        if(count==3)
        {
            ReturnTriangles<-c(ReturnTriangles,i)
        }
    }
    return(ReturnTriangles)
}