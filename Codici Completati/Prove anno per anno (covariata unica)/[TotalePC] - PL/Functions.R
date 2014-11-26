# FUNZIONI DI APPOGGIO
source("2013_SSR_AllFunctions.R")

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

Duplicated = function (x,y)
{
    Output<-NULL
    Row<-NULL
    for(i in 1:(length(x)-1))
    {
        for(j in (i+1):length(x))
        {
            if((x[i]==x[j])&&(y[i]==y[j]))
            {
                Output<-rbind(Output,c(x[i],y[i]))
                Row<-rbind(Row,c(i,j))
            }
        }
    }
    Out<-list(Output,Row)
    names(Out)<-c("Points","Rows")
    return(Out)
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


smooth.FEM.fd.CovarCI = function (data,desmat,fdobj,lambda,CI_level)
{
    # SMOOTH.FEM.FD.Covar Compute a solution for a FELspline problem 
    #
    #     Arguments:
    # FELSPLOBJ a FELspline object, constructed by SOLVE_FELSPLINE.
    # LAMBDA    a scalar smoothing parameter
    # DATA      (optional) a n-by-2 new set of observations.  DATA(:,1)
    #           indexes the points (FELSPLOBJ.POINTS) at which the 
    #           values in DATA(:,2) were observed.
    # DESMAT    a n-by-q design matrix
    #
    #     Output:
    # FELSPLOBJ  ...  A FD object of the FEM type defined by the coefficient
    #                 vector resulting from smoothing
    # LAPLACEFD  ...  A FD object of the FEM type for the value of the 
    #                 Laplace operator if order == 2, or empty if order == 1
    #
    #
    #  Last modified 8 February 2011 by Laura Sangalli
    
    
    
    #  check arguments
    
    if (!is.fd(fdobj))
    {   stop('FDOBJ is not a FD object')
    }
    
    if (!is.numeric(lambda))
    {  
        stop('LAMBDA is not numeric')
    }  else if (length(lambda) != 1)
    {   
        stop('LAMBDA is not a scalar')
    }
    
    #  check data argument
    
    if (is.null(data))
    {   data=getdata(fdobj)
    } else 
    {if (length(dim(data))>0)
    {if (dim(data)[[2]]!=2)
    {if (dim(data)[[1]]!=2)
    {stop('DATA is not a n-by-2 array')}   else
    {data=t(data)}
    }
    } else {stop('DATA is not a n-by-2 array')}
    }
    
    
    
    #  Construct penalty matrix and 'b' vector for Ax=b.
    
    
    basisobj = fdobj$basis
    
    numnodes = dim(basisobj$params$nodes)[[1]]
    
    nodeStruct = basisobj$params
    
    #Numero di covariate
    q=dim(desmat)[[2]]
    #Numero di dati
    np = length(data[,1])
    
    #  ---------------------------------------------------------------
    # construct mass matrix K0 
    #  ---------------------------------------------------------------
    
    K0 = mass(nodeStruct)
    
    #  ---------------------------------------------------------------
    # construct stiffness matrix K0 and K1.
    #  ---------------------------------------------------------------
    
    K1 = stiff1(nodeStruct)
    
    #  ---------------------------------------------------------------
    # construct the penalty matrix P with ones on diagonal at data points
    #  ---------------------------------------------------------------
    # Projection matrix
    # ATTENZIONE: get the projection matrix without computing it as below here, via decomposition
    
    desmatprod = ( solve( t(desmat) %*% desmat ) ) %*% t(desmat)
    
    H= desmat %*% desmatprod
    
    L=matrix(0,numnodes,numnodes)
    L[data[,1],data[,1]]=diag(1,length(data[,1]))-H
    
    
    #  ---------------------------------------------------------------
    # construct vector b for system Ax=b
    #  ---------------------------------------------------------------
    
    b            = matrix(numeric(numnodes*2),ncol=1)
    b[data[,1],] = L[data[,1],data[,1]] %*% data[,2]
    
    
    #  ---------------------------------------------------------------
    # construct matrix A for system Ax=b.
    #  ---------------------------------------------------------------
    
    A  = rbind(cbind(L, -lambda*K1), cbind(K1, K0))
    
    
    # solve system
    bigsol = solve(A,b)
    u = bigsol[1:numnodes,]
    s = bigsol[(numnodes+1):(2*numnodes),]
    
    
    #Costruzione degli intervalli di confidenza
    #prima ricavo fhat
    fnhat    = bigsol[1:np]
    betahat  = desmatprod %*% (data[,2]-fnhat)
    zhat     = desmat %*% betahat + fnhat
    #Ottimo, ricavato fhat e quindi betahat e zhat
    
    #Ora devo modificare A
    A = solve( L + lambda * K1 %*% solve(K0) %*% K1)
    
    EDF = q + sum(diag(( A[data[,1],data[,1]] %*% L[data[,1],data[,1]] )))
    
    sigmahat2 = t(data[,2] - zhat) %*% (data[,2] - zhat) / ( np - EDF)
    
    varhatbetahat= sigmahat2 * ( solve( t(desmat) %*% desmat )  +  desmatprod %*% ( A[data[,1],data[,1]] %*% L[data[,1],data[,1]] %*% A[data[,1],data[,1]] %*% t(desmatprod)))
    
    approxCIbeta=matrix(0,nrow=length(CI_level),ncol=3)
    
    for(i in 1:length(CI_level))
    {
        CI_l = betahat - qnorm(1-(1-CI_level[i])/2) * sqrt(diag(varhatbetahat))
        CI_u = betahat + qnorm(1-(1-CI_level[i])/2) * sqrt(diag(varhatbetahat))
        approxCIbeta[i,] = cbind(CI_l, betahat, CI_u)
    }
    colnames(approxCIbeta)=c("low","betahat","up")
    rownames(approxCIbeta)=CI_level
    
    # Make FELspline object
    
    felsplobj  = fd(u, basisobj)
    laplacefd = fd(s, basisobj)
    
    reslist=list(felsplobj=felsplobj,laplacefd=laplacefd,approxCIbeta=approxCIbeta)
    
    
    return(reslist)
}

plot.FEM.2D = function(fdobj, zlimits ,X=NULL, Y=NULL)  
{
    # PLOT  Plots a FEM object FDOBJ over a rectangular grid defined by 
    # vectors X and Y;
    #
    
    #  Last modified 4 February 2011 by Laura Sangalli.
    
    
    if (!is.fd(fdobj))
    {
        stop('FDOBJ is not an FD object')
    }
    
    coefmat = fdobj$coef
    
    basisobj = fdobj$basis
    
    params = basisobj$params
    
    p = params$p
    t = params$t
    t = t[,1:3]
    
    
    
    if (is.null(X))
    {
        xmin = min(p[,1])
        xmax = max(p[,1])
        nx   = 201
        X    = matrix(seq(xmin, xmax, len=nx),ncol=1)
    } else
    {
        xmin = min(X)
        xmax = max(X)
        nx   = length(X)
    }
    
    if (is.null(Y))
    {
        ymin = min(p[,2])
        ymax = max(p[,2])
        ny   = 201
        Y    = matrix(seq(ymin, ymax, len=ny),ncol=1)    
    } else
    {
        ymin = min(Y)
        ymax = max(Y)
        ny   = length(Y)
    }
    
    
    
    
    Xmat = X %*% matrix(1,nrow=1,ncol=ny)
    Ymat = matrix(1,nrow=nx,ncol=1) %*% t(Y)
    Xvec = NULL
    for (numc in 1:nx)
    {Xvec=c(Xvec,Xmat[,numc])}
    Yvec = NULL
    for (numc in 1:ny)
    {Yvec=c(Yvec,Ymat[,numc])}
    
    
    
    evalmat = eval.FEM.fd(Xvec, Yvec, fdobj)
    
    nsurf = dim(coefmat)[[2]]
    for (isurf in 1:nsurf)
    {
        evalmati = matrix(evalmat[,isurf],nrow=nx, ncol=ny, byrow=F)
        image(X,Y,evalmati,col=heat.colors(100), xlab="", ylab="", asp=1,zlim=zlimits)
        contour(X,Y,evalmati,add=T)
        if (nsurf > 1)
        {pause}
    }
    
    
}