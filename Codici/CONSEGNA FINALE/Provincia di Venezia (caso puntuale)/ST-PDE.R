
##### LOAD PACKAGES #####
library(fda)




##### TEMPORARY FUNCTIONS #####

basisfd = function (type, rangeval, nbasis, params, 
                    dropind = vector("list", 0), quadvals = vector("list", 0), 
                    values = vector("list", 0), basisvalues = vector("list", 0)) 
{
    #  Modified on 23 June 2010 by Laura Sangalli
    #  from R code in fda package.
    
    if (nargs() == 0) {
        type <- "bspline"
        rangeval <- c(0, 1)
        nbasis <- 2
        params <- vector("list", 0)
        dropind <- vector("list", 0)
        quadvals <- vector("list", 0)
        values <- vector("list", 0)
        basisvalues <- vector("list", 0)
        basisobj <- list(type = type, rangeval = rangeval, nbasis = nbasis, 
                         params = params, dropind = dropind, quadvals = quadvals, 
                         values = values, basisvalues = basisvalues)
        oldClass(basisobj) <- "basisfd"
        return(basisobj)
    }
    if (class(type) == "basisfd") {
        basisobj <- type
        return(basisobj)
    }
    if (type == "bspline" || type == "Bspline" || type == "spline" || 
            type == "Bsp" || type == "bsp") {
        type = "bspline"
    }
    else if (type == "con" || type == "const" || type == "constant") {
        type = "const"
    }
    else if (type == "exp" || type == "expon" || type == "exponential") {
        type = "expon"
    }
    else if (type == "Fourier" || type == "fourier" || type == 
                 "Fou" || type == "fou") {
        type = "fourier"
    }
    else if (type == "mon" || type == "monom" || type == "monomial") {
        type = "monom"
    }
    else if (type == "polyg" || type == "polygon" || type == 
                 "polygonal") {
        type = "polyg"
    }
    else if (type == "poly" || type == "pol" || type == "polynom" || 
                 type == "polynomial") {
        type = "polynom"
    }
    else if (type == "pow" || type == "power") {
        type = "power"
    }
    else if (type == "FEM") {
        type = "FEM"
    }
    else {
        type = "unknown"
    }
    if (type == "unknown") {
        stop("'type' unrecognizable.")
    }
    if (missing(quadvals)) 
        quadvals <- vector("list", 0)
    else if (!(length(quadvals) == 0 || is.null(quadvals))) {
        nquad <- dim(quadvals)[1]
        ncol <- dim(quadvals)[2]
        if ((nquad == 2) && (ncol > 2)) {
            quadvals <- t(quadvals)
            nquad <- dim(quadvals)[1]
            ncol <- dim(quadvals)[2]
        }
        if (nquad < 2) 
            stop("Less than two quadrature points are supplied.")
        if (ncol != 2) 
            stop("'quadvals' does not have two columns.")
    }
    if (!(length(values) == 0 || missing(values) || is.null(values))) {
        n <- dim(values)[1]
        k <- dim(values)[2]
        if (n != nquad) 
            stop(paste("Number of rows in 'values' not equal to number of", 
                       "quadrature points."))
        if (k != nbasis) 
            stop(paste("Number of columns in 'values' not equal to number of", 
                       "basis functions."))
    }
    else values <- vector("list", 0)
    if (!(length(basisvalues) == 0 || missing(basisvalues) || 
              !is.null(basisvalues))) {
        if (!is.list(basisvalues)) 
            stop("BASISVALUES is not a list object.")
        sizevec <- dim(basisvalues)
        if (length(sizevec) != 2) 
            stop("BASISVALUES is not 2-dimensional.")
        for (i in 1:sizevec[1]) {
            if (length(basisvalues[[i, 1]]) != dim(basisvalues[[i, 
                                                                2]])[1]) 
                stop(paste("Number of argument values not equal number", 
                           "of values."))
        }
    }
    else basisvalues <- vector("list", 0)
    if (missing(dropind)) 
        dropind <- vector("list", 0)
    if (type == "fourier") {
        paramvec <- rangeval[2] - rangeval[1]
        period <- params[1]
        if (period <= 0) 
            stop("Period must be positive for (a Fourier basis")
        params <- period
        if ((2 * floor(nbasis/2)) == nbasis) 
            nbasis <- nbasis + 1
    }
    else if (type == "bspline") {
        if (!missing(params)) {
            nparams <- length(params)
            if (nparams > 0) {
                if (params[1] <= rangeval[1]) 
                    stop("Smallest value in BREAKS not within RANGEVAL")
                if (params[nparams] >= rangeval[2]) 
                    stop("Largest value in BREAKS not within RANGEVAL")
            }
        }
    }
    else if (type == "expon") {
        if (length(params) != nbasis) 
            stop("No. of parameters not equal to no. of basis fns for (exponential basisobj$")
    }
    else if (type == "polyg") {
        if (length(params) != nbasis) 
            stop("No. of parameters not equal to no. of basis fns for (polygonal basisobj$")
    }
    else if (type == "power") {
        if (length(params) != nbasis) 
            stop("No. of parameters not equal to no. of basis fns for (power basisobj$")
    }
    else if (type == "const") {
        params <- 0
    }
    else if (type == "monom") {
        if (length(params) != nbasis) 
            stop("No. of parameters not equal to no. of basis fns for (monomial basisobj$")
    }
    else if (type == "polynom") {
        if (length(params) > 1) 
            stop("More than one parameter for (a polynomial basisobj$")
    }
    else if (type == "FEM") {
        if (!is.list(params)) 
            stop("params is not a list")
    }
    else stop("Unrecognizable basis")
    obj.call <- match.call()
    basisobj <- list(call = obj.call, type = type, rangeval = rangeval, 
                     nbasis = nbasis, params = params, dropind = dropind, 
                     quadvals = quadvals, values = values, basisvalues = basisvalues)
    oldClass(basisobj) <- "basisfd"
    basisobj
}



makenodes=function(p,t,order=2)
{
    #MAKENODES produces:
    #  a matrix NODES containing coordinates for all of the nodes to be used, 
    #  a matrix NODEINDEX defining which nodes correspond to each element.  
    #  If NORDER is 2, the midpoint of each edge is computed and added to 
    #  POINTS to obtain matrix NODES. 
    #  The row index of that midpoint is then added to the rows of TRIANGLES 
    #  containing that edge to define NODEINDEX.
    #  If NORDER is 1, nodes corresponds to vertices and NODEINDEX is 
    #  identical to TRIANGLES.
    #
    #  nvert:  number of vertices
    #  nele:   number of triangles or elements
    #
    #Input: POINTS is an nvert by 2 matrix containing the x and y
    #   coordinates of each of the nvert points in the right-hand rule mesh.
    #   POINTS = P' where P is the points matrix for pde
    #   The call can use P directly (see below).
    #   TRIANGLES is T(1:3,:)' where T is the triangle index matrix from pde.
    #   Vertices must be numbered in the counterclockwise direction.
    #   NORDER is the order of elements, and may be either 1 or 2 (default)
    #
    #Output: 
    #     NODES:  a numnodes*2 matrix whose i'th row contains
    #       the coordinates of the i'th nodal variable.
    #       Nodes for the second order element consist of vertices and 
    #       midpoints of edges, that is, 6 per triangle.
    #       The first NVER rows of NODES is POINTS, and the remainder are
    #       the edge midpoints.
    #       Nodes for the first order element consist of only vertices.
    #
    #     NODEINDEX:  for NORDER == 2, an nele*6 matrix whose i'th row
    #       contains the row numbers (in NODES) of the
    #       nodal variables defining the i'th finite 
    #       element.  If the i'th row of FMESH is [V1 V2 V3]
    #       then the i'th row of nodeindex is
    #       [V1 V(12) V2 V(23) V3 V(31)], where Vi is the
    #       row number of the i'th point and V(ij) is the 
    #       row number of the midpoint of the edge defined
    #       by the i'th and j'th points.
    #       If NORDER == 1, NODEINDEX is TRIANGLES.
    #
    #  Last modified 4 February 2011 by Laura Sangalli.
    
    #
    #  The first rows of nodes are the vertices
    
    
    if  (dim(p)[[2]]>2)
    {
        nodes=t(p)
        t=t(t[1:3,])
    } else
    {
        nodes = p
    }
    
    
    nele = dim(t)[[1]]
    nver = dim(p)[[1]]
    
    
    
    
    Jvec   = matrix(0,nele,1)      #  vector of jacobian values
    metric = array(0,c(nele,2,2))  #  3-d array of metric matrices
    
    if (order == 2)
    {  
        
        rec  = matrix (0, nrow=nver, ncol=nver)     
        ############# ATTENZIONE! for the moment I am not using sparse matrix notation
        ind  = rbind( c(1,2), c(2,3), c(3,1) )
        nodeindex = matrix (0, nrow=nele, ncol=6) 
        nodeindex[,c(1,3,5)] = t[,c(1,2,3)]
        
        
        for (i in 1:nele)
        {for (j in 1:3)
        {if (rec[t[i,ind[j,1]],t[i,ind[j,2]]]==0)
        {nodes = rbind(nodes, as.vector(t(c(.5,.5))) %*% as.matrix(nodes[t[i,ind[j,]],]))
         rec[t[i,ind[j,1]],t[i,ind[j,2]]] = dim(nodes)[[1]]
         rec[t[i,ind[j,2]],t[i,ind[j,1]]] = dim(nodes)[[1]]
         nodeindex[i,2*j] = dim(nodes)[[1]]
        }  
        else
        {
            nodeindex[i,2*j] = rec[t[i,ind[j,1]],t[i,ind[j,2]]]
        }   
        }
        
        #  deviations of vertices 2 and 3 from vertex 1
        
        diff1x = nodes[nodeindex[i,3],1] - nodes[nodeindex[i,1],1]
        diff1y = nodes[nodeindex[i,3],2] - nodes[nodeindex[i,1],2]
        diff2x = nodes[nodeindex[i,5],1] - nodes[nodeindex[i,1],1]
        diff2y = nodes[nodeindex[i,5],2] - nodes[nodeindex[i,1],2]
        
        #  Jacobian or area of triangle
        
        Jvec[i] = diff1x*diff2y - diff2x*diff1y
        
        #  Compute controvariant transformation matrix
        
        Ael = matrix(c(diff2y, -diff1y, -diff2x,  diff1x),nrow=2,ncol=2,byrow=T)/Jvec[i]
        
        #  Compute metric matrix
        
        metric[i,,] = t(Ael)%*%Ael
        
        } 
    } else if (order == 1)
    {nodeindex = t[,1:3]
     
     for (i in 1:nele)
     { #  deviations of vertices 2 and 3 from vertex 1
         
         diff1x = nodes[nodeindex[i,2],1]-nodes[nodeindex[i,1],1]
         diff1y = nodes[nodeindex[i,2],2]-nodes[nodeindex[i,1],2]
         diff2y = nodes[nodeindex[i,3],2]-nodes[nodeindex[i,1],2]
         diff2x = nodes[nodeindex[i,3],1]-nodes[nodeindex[i,1],1]
         
         #  Jacobian or area of triangle
         
         Jvec[i] = diff1x*diff2y - diff2x*diff1y
         
         #  Compute contravariant transformation matrix
         
         Ael = matrix(c(diff2y, -diff1y, -diff2x,  diff1x),nrow=2,ncol=2,byrow=T)/Jvec[i]
         
         #  Compute metric matrix
         
         metric[i,,] = t(Ael)%*%Ael
     }
    }  else
    {stop("ORDER not 1 or 2")}
    
    
    nodeStruct <- list(order=order, nodes=nodes, nodeindex = nodeindex, J=Jvec, metric=metric)
    return(nodeStruct)
}


mass=function(nodeStruct)
{
    #MASS produces the mass matrix containing integrals of products of
    #  nodal functions.  
    #
    #Input: NODESTRUCT is a struct object produced by function makenodes.
    #    It contains:
    #        ORDER     ... The order of the element (1 or 2)
    #        NODES     ... Coordinates of node points
    #        NODEINDEX ... indices of node points for each element
    #        JVEC      ... Jacobian of the affine transformation of each
    #                      element to the master element
    #
    #Output: K0: the NNOD by NNOD matrix of sums of products of nodal basis
    #        functions.
    #        For each element i, the integral of the product 
    #        of the j'th and k'th shape functions over the i'th element is
    #        computed.  Then that value is the 
    #        (NODEINDEX(i,j),NODEINDEX(i,k))'th entry of the i'th elemental 
    #        mass matrix.
    
    #  Last modified 4 February 2011 by Laura Sangalli.
    
    #  retrieve arrays from nodeStruct
    
    order     = nodeStruct$order
    nodes     = nodeStruct$nodes
    nodeindex = nodeStruct$nodeindex
    Jvec      = nodeStruct$J
    
    nele  = dim(nodeindex)[[1]]
    nnod  = dim(nodes)[[1]]
    
    
    if (order ==2)
    {   
        #  the integrals of products of basis functions for master element:
        
        K0M = matrix( c( 6,  0, -1, -4, -1,  0,
                         0, 32,  0, 16, -4, 16,
                         -1,  0,  6,  0, -1, -4,
                         -4, 16,  0, 32,  0, 16,
                         -1, -4, -1,  0,  6,  0,
                         0, 16, -4, 16,  0, 32), ncol=6, nrow=6, byrow=T)/360
        
        #  assemble the mass matrix
        
        
        K0 = matrix(0,nrow=nnod,ncol=nnod)
        for (el in 1:nele)
        {  ind = nodeindex[el,]
           K0[ind,ind] = K0[ind,ind] + K0M * Jvec[el]
        }
        
    } 
    else if (order == 1)
    {  
        #  the integrals of products of basis functions for master element:
        
        K0M = matrix( c( 2,  1,  1,
                         1,  2,  1,
                         1 , 1,  2), ncol=3, nrow=3, byrow=T) / 24
        
        #  assemble the mass matrix
        
        K0 = matrix(0,nrow=nnod,ncol=nnod)
        for (el in 1:nele)
        {  ind = nodeindex[el,]
           K0[ind,ind] = K0[ind,ind] + K0M * Jvec[el]
        }
    } 
    else
    {stop("ORDER not 1 or 2")}
    
    K0
}





stiff1= function(nodeStruct)
{
    #STIFF1 produces the nnod*nnod stiffness matrix K1
    #defined (K1)jk = int(dpsik/da*dpsij/da + dpsik/db*dpsij/db).
    #
    #Input: NODESTRUCT is a struct object produced by function makenodes.
    #    It contains:
    #        ORDER     ... The order of the element (1 or 2)
    #        NODES     ... Coordinates of node points
    #        NODEINDEX ... indices of node points for each element
    #        JVEC      ... Jacobian of the affine transformation of each
    #                      element to the master element
    #        METRIC    ... The crossproduct of the inverse of the linear
    #                      part of the transformation
    #
    #Output: K1 is an nnod*nnod matrix out which is
    #        the sum of the nele element stiffness matrices
    #        and the penalty stiffness matrix.
    #        These i'th element matrix has (ij)'th element defined
    #        as follows:  
    #        Let psita and psitb be the partial derivatives of the
    #        t'th shape function with respect to a and b (1<=t<=6).
    #        Then the integral of the sum of products
    #        (psija*psika+psijb+psikb) over the i'th element is
    #        computed.  Then that value is assigned to the
    #        (nodeindex(i,j),nodeindex(i,k))'th entry of the i'th elemental 
    #        stiffness matrix and the other elements are given the value zero.
    #
    #
    #  Last modified 4 February 2011 by Laura Sangalli.
    
    #  retrieve arrays from nodeStruct
    
    order     = nodeStruct$order
    nodes     = nodeStruct$nodes
    nodeindex = nodeStruct$nodeindex
    Jvec      = nodeStruct$J
    metric    = nodeStruct$metric
    
    nele  = dim(nodeindex)[[1]]
    nnod  = dim(nodes)[[1]]
    
    
    
    K1   = matrix(0,nrow=nnod,ncol=nnod)      ############### ATTENZIONE: sparse matrix representation not used for the moment
    ############### 
    
    #  assemble the stiffness matrix
    
    if (order == 2)
    {   
        #  values of K1 for master elements
        
        KXX = matrix( c( 3, -4,  1,  0,  0,  0, 
                         -4,  8, -4,  0,  0,  0, 
                         1, -4,  3,  0,  0,  0, 
                         0,  0,  0,  8,  0, -8, 
                         0,  0,  0,  0,  0,  0, 
                         0,  0,  0, -8,  0,  8), ncol=6, nrow=6, byrow=T)/6
        
        KXY = matrix( c( 3,  0,  0,  0,  1, -4,
                         -4,  4,  0, -4,  0,  4,
                         1, -4,  0,  4, -1,  0,
                         0, -4,  0,  4,  4, -4,
                         0,  0,  0,  0,  0,  0,
                         0,  4,  0, -4, -4,  4), ncol=6, nrow=6, byrow=T)/6
        
        KYY = matrix( c( 3,  0,  0,  0,  1, -4,
                         0,  8,  0, -8,  0,  0,
                         0,  0,  0,  0,  0,  0,
                         0, -8,  0,  8,  0,  0,
                         1,  0,  0,  0,  3, -4,
                         -4,  0,  0,  0, -4,  8), ncol=6, nrow=6, byrow=T)/6
        
        for (el in 1:nele)
        {
            ind    = nodeindex[el,]
            K1M = (metric[el,1,1]*KXX    + metric[el,1,2]*KXY +
                       metric[el,2,1]*t(KXY) + metric[el,2,2]*KYY)
            K1[ind,ind] = K1[ind,ind] + K1M*Jvec[el]
        }
    }    
    else if (order == 1)
    {   
        KXX = matrix( c(  1, -1,  0,
                          -1,  1,  0,
                          0,  0,  0), ncol=3, nrow=3, byrow=T) /2
        
        KXY = matrix( c(  1,  0, -1,
                          -1,  0,  1,
                          0,  0,  0), ncol=3, nrow=3, byrow=T) /2
        
        KYY = matrix( c(  1,  0, -1,
                          0,  0,  0,
                          -1,  0,  1), ncol=3, nrow=3, byrow=T) /2
        
        #  assemble the stiffness matrix
        
        for (el in 1:nele)
        {
            ind      = nodeindex[el,]
            K1M = (metric[el,1,1]*KXX    + metric[el,1,2]*KXY + 
                       metric[el,2,1]*t(KXY) + metric[el,2,2]*KYY)
            K1[ind,ind] = K1[ind,ind] + K1M*Jvec[el]
        }
    }  
    else
    {stop("ORDER not 1 or 2")}
    
    
    
    
    K1
}





eval.FEM.fd = function(X,Y,fdobj)
{
    # EVAL_FEM_FD evaluates the FEM fd object at points (X,Y)
    #
    #        arguments:
    # FELSPLOBJ a FELspline object
    # X         an array of x-coordinates.
    # Y         an array of y-coordinates.
    # 
    #        output:
    # EVALMAT   an array of the same size as X and Y containing the value of 
    #           FELSPLOBJ at (X,Y).
    
    
    #  Last modified 4 February 2011 by Laura Sangalli.
    
    
    #  check X
    
    if (!is.numeric(X))
    {
        stop('X is not a numerical array')
    } else
    {
        X  = matrix(X,ncol=1)     # treat X as a column vector
    }
    
    #  check Y
    
    if (!is.numeric(Y))
    {
        stop('Y is not a numerical array')
    } else if (length(Y)!=length(X))
    {
        stop('Y is not the same length as X')
    } else
    {
        Y  = matrix(Y,ncol=1)    # treat Y as a column vector
    }
    
    
    N = length(X)
    
    #  matrix of local coordinates 
    
    Pgpts = cbind((numeric(N)+1), X, Y)
    
    #  get basis
    
    basisobj = fdobj$basis
    
    #  get nodes and mesh
    
    params   = basisobj$params
    order    = params$order
    nodes    = params$nodes
    nodeindex = params$nodeindex
    p        = params$p
    t        = params$t
    t        = t[,1:3]
    
    Jvec      = params$J
    
    #  get the coefficient matrix
    
    coefmat = fdobj$coef
    
    nsurf = dim(coefmat)[[2]]
    
    # 1st, 2nd, and 3rd vertices of triangles
    
    if  (order == 2)
    {
        pg1 = nodes[nodeindex[,1],]
        pg2 = nodes[nodeindex[,3],]
        pg3 = nodes[nodeindex[,5],]
    }
    else if (order == 1)
    { 
        pg1 = nodes[nodeindex[,1],]
        pg2 = nodes[nodeindex[,2],]
        pg3 = nodes[nodeindex[,3],]
    }
    else
    {stop('ORDER is neither 1 nor 2.')}
    
    
    # denominator of change-of-coordinates change matrix
    
    ########modJac = Jvec
    ########
    ########modJacMat = modJac%*%t(ones3)
    
    
    
    # 1st, 2nd, and 3rd rows of change-of-coordinates matrix
    
    M1 = cbind(((pg2[,1]*(pg3[,2]/Jvec))-(pg3[,1]*(pg2[,2]/Jvec))),
               ((pg2[,2]-pg3[,2])/Jvec),                        
               ((pg3[,1]-pg2[,1])/Jvec))
    
    M2 = cbind(((pg3[,1]*(pg1[,2]/Jvec))-(pg1[,1]*(pg3[,2]/Jvec))),
               ((pg3[,2]-pg1[,2])/Jvec),    
               ((pg1[,1]-pg3[,1])/Jvec))
    
    M3 = cbind(((pg1[,1]*(pg2[,2]/Jvec))-(pg2[,1]*(pg1[,2]/Jvec))),
               ((pg1[,2]-pg2[,2])/Jvec),   
               ((pg2[,1]-pg1[,1])/Jvec))
    
    # identify element containing point in vector (Xvec(i),Yvec(i))
    # if no element contains a point, ind(i) is NaN
    
    
    tricoef = tricoefCal(p, t)
    
    ind = numeric(N)
    for (i in 1:N)
    {
        ind[i] = insideIndex(X[i,], Y[i,], p, t, tricoef)
    }
    
    
    #  interpolate values
    
    ones3 = matrix(c(1,1,1),ncol=1)
    evalmat = matrix(NA, nrow=N, ncol=nsurf)
    
    
    for (isurf in 1:nsurf)
    {for (i in 1:N)
    {indi = ind[i]
     if (!is.na(indi))
     {  
         #  change to barycentric coordinates
         baryc1 = (M1[indi,]*Pgpts[i,])%*%ones3
         baryc2 = (M2[indi,]*Pgpts[i,])%*%ones3
         baryc3 = (M3[indi,]*Pgpts[i,])%*%ones3
         if (order == 2)
         { 
             c1 = coefmat[nodeindex[indi,1],isurf]
             c2 = coefmat[nodeindex[indi,3],isurf]
             c3 = coefmat[nodeindex[indi,5],isurf]
             c4 = coefmat[nodeindex[indi,2],isurf]
             c5 = coefmat[nodeindex[indi,4],isurf]
             c6 = coefmat[nodeindex[indi,6],isurf]
             fval = c1*(2*baryc1^2 - baryc1) + 
                 c2*(2*baryc2^2 - baryc2) + 
                 c3*(2*baryc3^2 - baryc3) + 
                 c4*(4*baryc1* baryc2)    + 
                 c5*(4*baryc2* baryc3)    + 
                 c6*(4*baryc3* baryc1)
             evalmat[i,isurf] = fval
         } 
         else
         {
             c1 = coefmat[nodeindex[indi,1],isurf]
             c2 = coefmat[nodeindex[indi,2],isurf]
             c3 = coefmat[nodeindex[indi,3],isurf]
             fval = c1*baryc1 + c2*baryc2 + c3*baryc3
             evalmat[i,isurf] = fval
         }
     }
    }
    }
    
    
    
    evalmat
}







tricoefCal = function (p, t)
{
    #  TRICOEFCAL compute the coefficient matrix TRICOEF
    #  required to test of a point is indside a triangle
    
    #  Last modified 24 June 2010 by Laura Sangalli
    
    
    ntri   = dim(t)[[1]]
    
    #  compute coefficients for computing barycentric coordinates if
    #  needed
    
    tricoef = matrix(0, nrow=ntri, ncol=4)
    tricoef[,1] = p[t[,1],1]-p[t[,3],1]
    tricoef[,2] = p[t[,2],1]-p[t[,3],1]
    tricoef[,3] = p[t[,1],2]-p[t[,3],2]
    tricoef[,4] = p[t[,2],2]-p[t[,3],2]
    detT = matrix((tricoef[,1]*tricoef[,4] - tricoef[,2]*tricoef[,3]),ncol=1)
    tricoef = tricoef/(detT %*% matrix(1,nrow=1,ncol=4))
    
    tricoef
}






insideIndex = function (X, Y, p, t, tricoef)
{
    #  insideIndex returns the index of the triangle containing the point
    # (X,Y) if such a triangle exists, and NaN otherwise.
    #  TRICOEF may have already been calculated for efficiency,
    #  but if the function is called with four arguments, it is calculated.
    
    
    #  Last modified 24 June 2010 by Laura Sangalli
    
    eps=2.2204e-016
    small = 10000*eps
    
    ntri   = dim(t)[[1]]
    indtri   = matrix(1:ntri,ncol=1)
    
    #  compute coefficients for computing barycentric coordinates if needed
    
    if (is.null(tricoef))
    {
        tricoef = matrix(0, nrow=ntri, ncol=4)
        tricoef[,1] = p[t[,1],1]-p[t[,3],1]
        tricoef[,2] = p[t[,2],1]-p[t[,3],1]
        tricoef[,3] = p[t[,1],2]-p[t[,3],2]
        tricoef[,4] = p[t[,2],2]-p[t[,3],2]
        detT = matrix((tricoef[,1]*tricoef[,4] - tricoef[,2]*tricoef[,3]),ncol=1)
        tricoef = tricoef/(detT %*% matrix(1,nrow=1,ncol=4))
    }
    
    #  compute barycentric coordinates
    
    r3 = X - p[t[,3],1]
    s3 = Y - p[t[,3],2]
    lam1 = ( tricoef[,4]*r3 - tricoef[,2]*s3)
    lam2 = (-tricoef[,3]*r3 + tricoef[,1]*s3)
    lam3 = 1 - lam1 - lam2
    
    #  test these coordinates for a triple that are all between 0 and 1
    
    int  = (-small <= lam1 & lam1 <= 1+small) & 
        (-small <= lam2 & lam2 <= 1+small) & 
        (-small <= lam3 & lam3 <= 1+small)
    
    #  return the index of this triple, or NaN if it doesn't exist
    
    indi = indtri[int]
    if (length(indi)<1)
    {
        ind = NA
    } else
    {ind = min(indi)
    }
    
    
    ind
}



##### FUNCTION FOR SPATIO-TEMPORAL ANALYSIS #####

create.FEM.basis = function(p, e, t, order, dl=NULL)
{
    # 
    #  CREATE.FEM.BASIS sets up a finite element basis for the analysis
    #  of spatial data.  It requires a triangular mesh as input.
    #  The triangular mesh is defined by a Delaunay triangulation of the domain 
    #  of the basis functions, along with, optionally, a decomposition of the 
    #  domain into subdomains (optional).
    #  The finite elements used for functional data analysis are first or second
    #  order Lagrangian elements.  These are triangles covering a region,
    #  and the basis system is piecewise polinomials (linear or quadratic). There is a basis
    #  function associated with each node in the system.
    #  When ORDER = 1 the basis system is piecewise linear and the nodes are the vertices of the triangles.
    #  When ORDER = 2 the basis system is piecewise quadratic and the nodes are points that are etiher 
    #  the vertices of the triangles or midpoints of edges of triangles.
    #
    #  Arguments:
    #  P  ...  The NBASIS by 2 matrix of vertices of triangles containing
    #          the X- and Y-coordinates of the vertices.
    #          P may also be provided as a 2 by NBASIS matrix.
    #  E  ...  The number of edges by 7 matrix defining the segments of the 
    #          boundary of the region which are also edges of the triangles
    #          adjacent to the boundary.  
    #          The values in matrix E as follows:
    #          Columns 1 and 2:  the indices of the vertices in matrix P of
    #          the starting and ending points of the edges
    #          Columns 3 and 4:  the starting and ending parameter values.  
    #          These are initialized to 0 and 1.
    #          Column  5:        the edge segment number 
    #          Columns 6 and 7:  the left and right-hand side subdomain numbers
    #          E may also be provided as a 7 by number of edges matrix.
    #  T  ...  The 4 by no. of triangles matrix specifying triangles and 
    #          their properties:
    #          Columns 1 to 3:   the indices in P of the vertices of each 
    #          triangle in counter-clockwise order.
    #          Column  4:        the subdomain number of the triangle
    #          T may also be provided with 3 columns, in which case the
    #          final column is set to ones.  It may also be provided as
    #          either a 3 or 4 by number of triangles matrix.
    #  ORDER.  Order of elements, which may be either 1 or 2 (2 is default)
    #  DL ...  Decomposed geometry object.  This defines the characteristics of the subdomains
    #          and the boundary segments.  It is optional.
    #  
    #  Returns:
    #  An object of the basis class with parameters stored in member params,
    #  which in this case is a struct object with members p, e and t.
    
    #  Last modified 4 February 2011 by Laura Sangalli.
    
    
    #  check t for having 3 rows or columns, and add 1's to
    #  make the dimension 4
    
    t<-as.matrix(t)
    dimnames(t)=NULL
    if  (dim(t)[[2]] == 3)
    {    #  add a column of 1's
        t = cbind(t, numeric(dim(t)[[1]])+1 )
    }
    
    if  (dim(t)[[1]] == 3)
    {   #  add a row of 1's
        t = rbind(t, numeric(dim(t)[[1]])+1 )
    }
    
    #  check dimensions of P, E and T and transpose if necessary
    
    if  (is.null(e))
    {
        if ((dim(p)[[1]] != 2 && dim(p)[[2]] != 2) || (dim(t)[[1]] != 4 && dim(t)[[2]] != 4))
        { stop('Dimensions of at least one of P and T are not correct.')
        }
    } else
    {
        if ((dim(p)[[1]] != 2 && dim(p)[[2]] != 2) || (dim(e)[[1]] != 7 && dim(e)[[2]] != 7) || (dim(t)[[1]] != 4 && dim(t)[[2]] != 4))
        { stop('Dimensions of at least one of P, E and T are not correct.')
        }
    }
    
    
    if (dim(p)[[2]] != 2 && dim(p)[[1]] == 2)
    {   p = t(p)
    }
    
    if  (!is.null(e))
    { if (dim(e)[[2]] != 7 && dim(e)[[1]] == 7)
    {   e = t(e)
    }
    }
    
    if (dim(t)[[2]] != 4 && dim(t)[[1]] == 4)
    {   t = t(t)
    }
    
    type <- "FEM"
    
    #  Argument rangeval is not needed for an FEM basis since domain
    #  boundary is defined by the triangular mesh.
    rangeval = NULL  
    
    #  The number of basis functions corresponds to the number of vertices
    #  for order = 1, and to vertices plus edge midpoints for order = 2
    
    #  set up the nodal points and corresponding mesh:  
    #    this is p' and t(1:3,:)' in the order 1 case, but
    #    includes mid-points in the order 2 case
    
    
    nodeStruct= makenodes(p,t[,1:3],order)
    
    #  The params argument is a struct object 
    petstr=NULL
    petstr$p        = p
    petstr$e        = e
    petstr$t        = t
    petstr$order     = order
    petstr$nodes     = nodeStruct$nodes
    petstr$nodeindex = nodeStruct$nodeindex
    petstr$J         = nodeStruct$J
    petstr$metric    = nodeStruct$metric
    
    if (!is.null(dl))
    {
        petstr$dl       = dl
    } 
    
    params = petstr
    
    nbasis = dim(petstr$nodes)[[1]]
    
    basisobj = basisfd(type, rangeval, nbasis, params)
    
    basisobj
}


ST.Smooth = function(Data,DesMat=NA,SpacePoints,SpaceBasisObj,TimePoints,TimeBasisObj,LogLambdaS,LogLambdaT)
{
    # ST.SMOOTH Computes GCV and solution for ST-PDE model
    # This function can be used for model with covariates or without covariates
    # If covariates are not used in model, DesMat must be a NA
    # If space/time smoothing parameter is not single, GCV will be computed
    # and best value will be found
    #
    # Arguments:
    # DATA              nm * 1 vector with data for each spatio-temporal observation.
    #                   if i is Space index and j Time Index 
    #                   DATA MUST BE ORDERED FIRST VARYING TIME INDEX, THAN SPACE INDEX
    #                   z11, z12, z13, ..., z1m, z21, z22, ..., z2m, ..., znm
    # DESMAT            nm * p matrix, each row is the vector of p covariates
    #                   associated to the corresponding observation in Data
    #                   For the rows, use the same order of Data.
    #                   IF YOU WANT TO USE MODEL WITHOUT COVARIATES, NA IN REQUIRED!
    #                   NA is default value
    # SPACEPOINTS       n*2 matrix with spatial points
    # SPACEBASISOBJ     space basis object created with create.FEM.basis
    # TIMEPOINTS        m-dimensional vector with time points
    # TIMEBASISOBJ      time basis object created with create.bspline.basis (fda)
    # LOGLAMBDAS        log10 of space smoothing parameter
    #                   If is given more than a single value, GCV will be computed and
    #                   best value will be found
    # LOGLAMBDAT        log10 of time smoothing parameter
    #                   If is given more than a single value, GCV will be computed and
    #                   best value will be found
    
    #
    # Output:
    # A SolutionObj class object with this elements:
    # CHAT              cHat vector for function estimate
    #                   elements are ordered like Data (but index does not refers to
    #                   given observation, but to basis functions)
    # BETAHAT           estimation of Beta
    # BESTLOGLAMBDA     Best values for (LogLambdaS,LogLambdaT).
    #                   If GCV isn't computed, these are given values
    #                   If GCV is computed, these are best values
    # EDF               (if GCV isn't computed)
    #                   Equivalent degrees of freedom (trace of smoothing matrix)
    # SIGMA2HAT         (if GCV isn't computed)
    #                   Variance estimated
    # EDFMATRIX         (if GCV is computed)
    #                   Equivalent degrees of freedom for each GCV case
    # SIGMA2HATMATRIX   (if GCV is computed)
    #                   Variance estimated for each GCV case
    # GCVMATRIX         (if GCV is computed)
    #                   GCV index value for each GCV case
    #
    # Last modified on 15 April 2015 by Gabriele Mazza
    
    # SPACE OBJECTS
    if (class(SpaceBasisObj)!="basisfd")
    {  
        stop('Space basis object not valid')
    }
    if (SpaceBasisObj$type!="FEM")
    {  
        stop('Space basis are not FE')
    }
    if (!is.numeric(LogLambdaS))
    {  
        stop('LOGLAMBDAS is not numeric')
    }
    if (dim(SpacePoints)[2]!=2)
    {  
        stop('Space Points are not bidimensional')
    }
    
    # TIME OBJECTS
    if(class(TimeBasisObj)!="basisfd")
    {  
        stop('Time basis object not valid')
    }
    if(TimeBasisObj$type!="bspline")
    {  
        stop('Time basis are not bspline')
    }
    if (!is.numeric(LogLambdaT))
    {  
        stop('LOGLAMBDAT is not numeric')
    }
    
    # DATA & DESIGN MATRIX
    Data<-as.matrix(Data)
    # Check dimensions
    if ((dim(Data)[2]) != 1)
    {   
        stop('Data must be a vector')
    }
    if ((dim(Data)[1]) != ((dim(SpacePoints)[1])*(length(as.vector(TimePoints)))))
    {   
        stop('Incorrect number of data')
    }
    if(!is.na(DesMat)[1])
    {
        DesMat<-as.matrix(DesMat)
        # Check dimensions
        if ((dim(Data)[1]) != (dim(DesMat)[1]))
        {   
            stop('DesMat dimensions not correct')
        }
    }
    
    
    
    ### CREATE B MATRIX ###
    
    # I create Psi and Phi and use kronecker product
    Phi<-NULL
    for(i in 1:SpaceBasisObj$nbasis)
    {
        # basis function i
        coeff<-numeric(SpaceBasisObj$nbasis)
        coeff[i]<-1
        fdObj<-fd(coeff,SpaceBasisObj)
        eval<-eval.FEM.fd(SpacePoints[,1],SpacePoints[,2],fdObj)
        Phi<-cbind(Phi,eval)
    }
    Psi=eval.basis(TimePoints,TimeBasisObj)
    dimnames(Psi)[[2]]<-NULL
    B=kronecker(Phi,Psi)
    #Remove useless matrixes
    rm(Phi,Psi)
    
    ### CREATE OBJECTS FOR PTIME & PSPACE ###
    # P matrix has LambdaS and LambdaT
    # So I cant' compute it now (in GCV P matrix change)
    # I build some matrices that will be used in all the algorithm
    K0 = mass(SpaceBasisObj$params)
    K1 = stiff1(SpaceBasisObj$params)
    Pspace = K1%*%solve(K0)%*%K1
    Ptime=eval.penalty(TimeBasisObj, 2)
    KronSpace=kronecker(Pspace,diag(dim(Ptime)[1]))
    KronTime=kronecker(diag(dim(Pspace)[1]),Ptime)
    #Remove useless matrixes
    rm(K0,K1,Pspace,Ptime)
    
    
    ### MODEL WITHOUT COVARIATES ###
    
    if(is.na(DesMat)[1])
    {
        # Compute GCV?
        if(length(LogLambdaS)*length(LogLambdaT)!=1)
        {
            ### GCV ###
            print("GCV phase...")
            n=dim(Data)[1]

            # Give names to matrices... Using LogLambda values
            # So output will be better
            NameS<-NULL
            NameT<-NULL
            for(i in 1:length(LogLambdaS))
            {
                NameS<-c(NameS,paste("S ",round(LogLambdaS[i],2),sep=""))
            }
            for(i in 1:length(LogLambdaT))
            {
                NameT<-c(NameT,paste("T ",round(LogLambdaT[i],2),sep=""))
            }
            List<-list(NameS,NameT)
            # Generalized Cross Validation Matrix
            GCVMatrix<-matrix(nrow=length(LogLambdaS),ncol=length(LogLambdaT),dimnames=List)
            # Effective Degrees of Freedom Matrix
            EDFMatrix<-matrix(data=0,nrow=length(LogLambdaS),ncol=length(LogLambdaT),dimnames=List)
            # Sigma^2 Hat Matrix
            Sigma2HatMatrix<-matrix(data=0,nrow=length(LogLambdaS),ncol=length(LogLambdaT),dimnames=List)
            
            BestS<-NULL
            BestT<-NULL
            BestGCV<-Inf
            
            # Progress bar (so in Rstudio i can see the progress)
            ntot=length(LogLambdaS)*length(LogLambdaT)+1
            k=1
            pb <- txtProgressBar(min = 1, max = ntot, style = 3)
            setTxtProgressBar(pb, k)
            
            for(i in 1:length(LogLambdaS))
            {
                LambdaS=10^LogLambdaS[i]
                for(j in 1:length(LogLambdaT))
                {
                    LambdaT<-10^LogLambdaT[j]
                    # Compute Hat Matrix
                    P=LambdaS*KronSpace+LambdaT*KronTime
                    H=B%*%solve(t(B)%*%B+P)%*%t(B)
                    
                    # Effective Degrees of Freedom
                    EDFMatrix[i,j]=sum(diag(H))
                    DataHat<-H%*%Data
                    tmp=t((Data-DataHat))%*%(Data-DataHat)
                    # Generalized Cross Validation
                    GCVMatrix[i,j]=(n/((n-EDFMatrix[i,j])^2))*tmp
                    # Sigma2Hat
                    Sigma2HatMatrix[i,j]=(1/(n-EDFMatrix[i,j]))*tmp
                    # Is the Best?
                    if(GCVMatrix[i,j]<BestGCV)
                    {
                        BestS=LogLambdaS[i]
                        BestT=LogLambdaT[j]
                        BestGCV=GCVMatrix[i,j]
                    }
                    k=k+1
                    setTxtProgressBar(pb, k)
                }
            }
            #Remove useless matrixes
            rm(P,H)
            # GCV is optimum on grid extremes?
            if((BestS-LogLambdaS[1])*(BestS-LogLambdaS[length(LogLambdaS)])*(BestT-LogLambdaT[1])*(BestT-LogLambdaT[length(LogLambdaT)])==0)
            {
                print("WARNING: GCV optimum on grid extremes")
            }
            print("End of GCV phase.")
        } else
        {
            # IF GCV isn't computed, best values are given values
            BestS=LogLambdaS
            BestT=LogLambdaT
        }
        
        ### SOLUTION ### 
        print(paste("Analysis will be computed with Lambda=(10^",BestS,";10^", BestT,")",sep=""))
        LambdaS=10^BestS
        LambdaT=10^BestT
        print("Solving system...")
        P=LambdaS*KronSpace+LambdaT*KronTime
        Temp=solve(t(B)%*%B+P)%*%t(B)
        cHat=Temp%*%Data
        print("Done.")
        #Remove useless matrixes
        rm(P)
        if(length(LogLambdaS)*length(LogLambdaT)!=1)
        {
            # If GCV is computed, I have matrices
            SolutionObj<-list(cHat=cHat,BestLogLambda=c(BestS,BestT),GCVMatrix=GCVMatrix,EDFMatrix=EDFMatrix,Sigma2HatMatrix=Sigma2HatMatrix)
        } else
        {
            # If GCV is not computed, I have to compute variance and EDFG
            H=B%*%Temp
            DataHat=H%*%Data
            EDF=sum(diag(H))
            Sigma2Hat=(1/((dim(Data)[1])-EDF))*t((Data-DataHat))%*%(Data-DataHat)
            SolutionObj<-list(cHat=cHat,BestLogLambda=c(BestS,BestT),EDF=EDF,Sigma2Hat=Sigma2Hat)
            #Remove useless matrixes
            rm(H,Temp)
        }        
    } else
    {
        
                
        ### MODEL WITH COVARIATES ###
                
        ### GCV ###
        if(length(LogLambdaS)*length(LogLambdaT)!=1)
        {
            print("GCV phase...")
            n=dim(Data)[1]
            # Give names to matrices... Using LogLambda values
            # So output will be better
            NameS<-NULL
            NameT<-NULL
            for(i in 1:length(LogLambdaS))
            {
                NameS<-c(NameS,paste("S ",round(LogLambdaS[i],2),sep=""))
            }
            for(i in 1:length(LogLambdaT))
            {
                NameT<-c(NameT,paste("T ",round(LogLambdaT[i],2),sep=""))
            }
            List<-list(NameS,NameT)
            # Generalized Cross Validation Matrix
            GCVMatrix<-matrix(nrow=length(LogLambdaS),ncol=length(LogLambdaT),dimnames=List)
            # Effective Degrees of Freedom Matrix
            EDFMatrix<-matrix(data=0,nrow=length(LogLambdaS),ncol=length(LogLambdaT),dimnames=List)
            # Sigma^2 Hat Matrix
            Sigma2HatMatrix<-matrix(data=0,nrow=length(LogLambdaS),ncol=length(LogLambdaT),dimnames=List)
            
            BestS<-NULL
            BestT<-NULL
            BestGCV<-Inf
            
            #Temporary matrix (to have improve efficiency)
            TempDM<-DesMat%*%solve(t(DesMat)%*%DesMat)%*%t(DesMat)
            Q=(diag(dim(DesMat)[1])-TempDM)
            Temp1<-t(B)%*%Q%*%B
            Temp2=t(B)%*%Q
            
            # progress bar... So in rstudio I can see the progress
            ntot=length(LogLambdaS)*length(LogLambdaT)+1
            k=1
            pb <- txtProgressBar(min = 1, max = ntot, style = 3)
            setTxtProgressBar(pb, k)
            for(i in 1:length(LogLambdaS))
            {
                LambdaS=10^LogLambdaS[i]
                for(j in 1:length(LogLambdaT))
                {
                    LambdaT=10^LogLambdaT[j]
                    P=LambdaS*KronSpace+LambdaT*KronTime
                    HSm=B%*%solve(P+Temp1)%*%Temp2
                    HCovar=TempDM%*%(diag(dim(HSm)[1])-HSm)
                    # Effective Degrees of Freedom
                    EDFMatrix[i,j]<-sum(diag(HSm))+sum(diag(HCovar))
                    
                    DataHat<-(HSm+HCovar)%*%Data
                    tmp=t((Data-DataHat))%*%(Data-DataHat)
                    # Generalized Cross Validation
                    GCVMatrix[i,j]=(n/((n-EDFMatrix[i,j])^2))*tmp
                    # Sigma2Hat
                    Sigma2HatMatrix[i,j]=(1/(n-EDFMatrix[i,j]))*tmp
                    # Is the Best?
                    if(GCVMatrix[i,j]<BestGCV)
                    {
                        BestS=LogLambdaS[i]
                        BestT=LogLambdaT[j]
                        BestGCV=GCVMatrix[i,j]
                    }
                    k=k+1
                    setTxtProgressBar(pb, k)
                }
            }
            #Remove useless matrixes
            rm(TempDM,Q,Temp1,Temp2,P,HSm,HCovar)
            # GCV is optimum on grid extremes?            
            if((BestS-LogLambdaS[1])*(BestS-LogLambdaS[length(LogLambdaS)])*(BestT-LogLambdaT[1])*(BestT-LogLambdaT[length(LogLambdaT)])==0)
            {
                print("WARNING: GCV optimum on grid extremes")
            }
            print("End of GCV phase.")
        } else
        {
            # if GCV is not computed, given values are best values
            BestS=LogLambdaS
            BestT=LogLambdaT
        }
        print(paste("Analysis will be computed with Lambda=(10^",BestS,";10^", BestT,")",sep=""))
        LambdaS=10^BestS
        LambdaT=10^BestT
        print("Solving system...")
        P=LambdaS*KronSpace+LambdaT*KronTime
        TempDM<-DesMat%*%solve(t(DesMat)%*%DesMat)%*%t(DesMat)
        Q=(diag(dim(TempDM)[1])-TempDM)
        cHat=solve(P+t(B)%*%Q%*%B)%*%t(B)%*%Q%*%Data
        BetaHat=solve(t(DesMat)%*%DesMat)%*%t(DesMat)%*%(Data-B%*%cHat)
        print("Done.")
        
        if(length(LogLambdaS)*length(LogLambdaT)!=1)
        {
            # If GCV is computed, I have matrices
            SolutionObj<-list(cHat=cHat,BetaHat=BetaHat,BestLogLambda=c(BestS,BestT),GCVMatrix=GCVMatrix,EDFMatrix=EDFMatrix,Sigma2HatMatrix=Sigma2HatMatrix)
        } else
        {
            # If GCV is not computed, I have to compute variance and EDFG
            HSm=B%*%solve(P+t(B)%*%Q%*%B)%*%t(B)%*%Q
            HCovar=TempDM%*%(diag(dim(HSm)[1])-HSm)
            DataHat=(HSm+HCovar)%*%Data
            EDF=sum(diag(HSm+HCovar))
            Sigma2Hat=(1/((dim(Data)[1])-EDF))*t((Data-DataHat))%*%(Data-DataHat)
            SolutionObj<-list(cHat=cHat,BetaHat=BetaHat,BestLogLambda=c(BestS,BestT),EDF,Sigma2Hat=Sigma2Hat)
            #Remove useless matrixes
            rm(HSm,HCovar)
        }       
    }
    class(SolutionObj)="SolutionObj"
    return(SolutionObj)
}



ST.Eval = function(X,Y,Time,SpaceBasisObj,TimeBasisObj,SolutionObj)
{
    # ST.EVAL Evaluates the estimated function in some given spatio-temporal points
    # In both cases (with or without covariates) this function is used to
    # evaluate functional part (only f(p,t))
    #
    # Arguments:
    # X                 Vector of x coordinates (same length than Y and TIME)              
    # Y                 Vector of y coordinates (same length than X and TIME) 
    # TIME              Vector of time istants (same length than X and Y) 
    # SPACEBASISOBJ     space basis object created with create.FEM.basis
    # TIMEBASISOBJ      time basis object created with create.bspline.basis (fda)
    # SOLUTIONOBJ       Result of ST.Smooth (estimated function)
    
    #
    # Output:
    # RESULT            Values of f corresponding to gives space points and time istants
    # Last modified on 15 April 2015 by Gabriele Mazza
    
    # Check arguments
    if (class(SolutionObj)!="SolutionObj")
    {  
        stop('Not valid SolutionObj in FixedPointPlot')
    }
    if (class(SpaceBasisObj)!="basisfd")
    {  
        stop('Space basis object not valid')
    }
    if (SpaceBasisObj$type!="FEM")
    {  
        stop('Space basis are not FE')
    }
    if(class(TimeBasisObj)!="basisfd")
    {  
        stop('Time basis object not valid')
    }
    if(TimeBasisObj$type!="bspline")
    {  
        stop('Time basis are not bspline')
    }
    # Check X
    if (!is.numeric(X))
    {
        stop('X is not a numerical array')
    } else
    {
        X  = matrix(X,ncol=1)     # treat X as a column vector
    }
    
    N = length(X)
    
    # Check Y
    if (!is.numeric(Y))
    {
        stop('Y is not a numerical array')
    } else if (length(Y)!=N)
    {
        stop('Y is not the same length as X')
    } else
    {
        Y  = matrix(Y,ncol=1)    # treat Y as a column vector
    }
    # Check Time
    if (!is.numeric(Time))
    {
        stop('TIME is not a numerical array')
    } else if (length(Time)!=N)
    {
        stop('TIME is not the same length as X and Y')
    }
    # TIME must be a vector, not a matrix
    
    # I need Phi matrix in used time points
    Psi=eval.basis(Time,TimeBasisObj)
    
    Result<-numeric(N)
    
    CMat<-matrix(data=SolutionObj$cHat,nrow=SpaceBasisObj$nbasis,ncol=TimeBasisObj$nbasis,byrow=T)
    for (j in 1:(dim(CMat)[2]))
    {
        fdObj=fd(CMat[,j],SpaceBasisObj)
        temp<-eval.FEM.fd(X,Y,fdObj)
        Result=Result+temp*Psi[,j]
    }
    
    return(Result)
}




ST.IC = function(Data,DesMat,SpaceBasisObj,TimeBasisObj,LogLambdaS,LogLambdaT,Alpha=0.05,Correct=TRUE)
{
    # ST.IC Computes confidence intervals for BetaHat components
    # Confidence intervals are approssimated.
    # If Beta han more than one component, confidence can be corrected with Bonferroni
    #
    # Arguments:
    # DATA              nm * 1 vector with data for each spatio-temporal observation.
    #                   if i is Space index and j Time Index 
    #                   DATA MUST BE ORDERED FIRST VARYING TIME INDEX, THAN SPACE INDEX
    #                   z11, z12, z13, ..., z1m, z21, z22, ..., z2m, ..., znm
    # DESMAT            nm * p matrix, each row is the vector of p covariates
    #                   associated to the corresponding observation in Data
    #                   For the rows, use the same order of Data.
    #                   IF YOU WANT TO USE MODEL WITHOUT COVARIATES, NA IN REQUIRED!
    #                   NA is default value
    # SPACEBASISOBJ     space basis object created with create.FEM.basis
    # LOGLAMBDAS        log10 of space smoothing parameter
    #                   If is given more than a single value, GCV will be computed and
    #                   best value will be found
    # LOGLAMBDAT        log10 of time smoothing parameter
    #                   If is given more than a single value, GCV will be computed and
    #                   best value will be found
    # ALPHA             Confidence (default 0.05)
    # CORRECT           Une Bonferroni correction? (default TRUE)
    #
    # Output:
    # An objact with the following elements:
    # APPROXCIBETA      cHat vector for function estimate
    #                   elements are ordered like Data (but index does not refers to
    #                   given observation, but to basis functions)
    # CORRECTEDALPHA    Alpha after correction (real alpha used for intervals).
    #                   If CORRECT parameter is false, CORRECTEDALPHA will be
    #                   equal to ALPHA
    # SIGMA2HAT         (if GCV isn't computed)
    #
    # Last modified on 15 April 2015 by Gabriele Mazza

    # Check arguments
    # Space
    if (class(SpaceBasisObj)!="basisfd")
    {  
        stop('Space basis object not valid')
    }
    if (SpaceBasisObj$type!="FEM")
    {  
        stop('Space basis are not FE')
    }
    if(class(TimeBasisObj)!="basisfd")
    {  
        stop('Time basis object not valid')
    }
    if(TimeBasisObj$type!="bspline")
    {  
        stop('Time basis are not bspline')
    }
    if (!is.numeric(LogLambdaS))
    {  
        stop('LOGLAMBDAS is not numeric')
    }
    if (length(LogLambdaS) != 1)
    {   
        stop('LOGLAMBDAS is not a scalar')
    }
    if (!is.numeric(LogLambdaT))
    {  
        stop('LOGLAMBDAT is not numeric')
    }
    if (length(LogLambdaT) != 1)
    {   
        stop('LOGLAMBDAT is not a scalar')
    }
    if (Alpha<0 | Alpha >1)
    {   
        stop('ALPHA not valid.')
    }
    
    
    # DATA & DESIGN MATRIX
    Data<-as.matrix(Data)
    # Check dimensions
    if ((dim(Data)[2]) != 1)
    {   
        stop('Data must be a vector')
    }
    if ((dim(Data)[1]) != ((dim(SpacePoints)[1])*(length(as.vector(TimePoints)))))
    {   
        stop('Incorrect number of data')
    }
    DesMat<-as.matrix(DesMat)
    # Check dimensions
    if ((dim(Data)[1]) != (dim(DesMat)[1]))
    {   
        stop('DesMat dimensions not correct')
    }
        
    ### CREATE B MATRIX ###
    
    # Now i build an fdobj for each basis function
    Phi<-NULL
    for(i in 1:SpaceBasisObj$nbasis)
    {
        # basis function i
        coeff<-numeric(SpaceBasisObj$nbasis)
        coeff[i]<-1
        fdObj<-fd(coeff,SpaceBasisObj)
        eval<-eval.FEM.fd(SpacePoints[,1],SpacePoints[,2],fdObj)
        Phi<-cbind(Phi,eval)
    }
    Psi=eval.basis(TimePoints,TimeBasisObj)
    dimnames(Psi)[[2]]<-NULL
    B=kronecker(Phi,Psi)
    #Remove useless matrixes
    rm(Phi,Psi)
    
    ### CREATE OBJECTS FOR PTIME & PSPACE ###
    
    K0 = mass(SpaceBasisObj$params)
    K1 = stiff1(SpaceBasisObj$params)
    Pspace = K1%*%solve(K0)%*%K1
    Ptime=eval.penalty(TimeBasisObj, 2)
    P=(10^LogLambdaS)*kronecker(Pspace,diag(dim(Ptime)[1]))+(10^LogLambdaT)*kronecker(diag(dim(Pspace)[1]),Ptime)
    #Remove useless matrixes
    rm(K0,K1,Pspace,Ptime)
    
    n=dim(Data)[1]
    
 
    Temp<-solve(t(DesMat)%*%DesMat)%*%t(DesMat)
    Q= (diag(dim(DesMat)[1])-DesMat%*%Temp)
    A= solve(t(B)%*%Q%*%B+P)%*%t(B)
    
    PSm=A%*%Q
    PCovar=Temp%*%(diag(dim(B)[1])-B%*%PSm)
    
    # First of all, I need BetaHat
    BetaHat=PCovar%*%Data
        
    # Then, I need Sigma2Hat
    HatMatrix=B%*%PSm+DesMat%*%PCovar
    DataHat=HatMatrix%*%Data
    
    Sigma2Hat=(1/(n-sum(diag(HatMatrix))))*t(Data-DataHat)%*%(Data-DataHat)
    
    # Then, I need Variance Matrix
    VarMatrix=
        Sigma2Hat*(
            solve(t(DesMat)%*%DesMat)+
                Temp%*%B%*%A%*%Q%*%t(A)%*%t(B)%*%t(Temp)
        )
    
    if(Correct)
    {
        CorrectedAlpha=Alpha/(dim(BetaHat)[1])
        
        CI_l = BetaHat - qnorm(1-CorrectedAlpha/2) * sqrt(diag(VarMatrix))
        CI_u = BetaHat + qnorm(1-CorrectedAlpha/2) * sqrt(diag(VarMatrix))
        ApproxCIBeta = cbind(CI_l, BetaHat, CI_u)
    } else
    {
        CorrectedAlpha=NULL
        
        CI_l = BetaHat - qnorm(1-Alpha/2) * sqrt(diag(VarMatrix))
        CI_u = BetaHat + qnorm(1-Alpha/2) * sqrt(diag(VarMatrix))
        ApproxCIBeta = cbind(CI_l, BetaHat, CI_u)
    }
    
    # Now I create names for rows and colums... 
    colnames(ApproxCIBeta)=c("Low","BetaHat","Up")
    rows<-NULL
    for(i in 1:(dim(BetaHat)[1]))
    {
        rows<-c(rows,paste("Beta",i,sep=""))
    }
    rownames(ApproxCIBeta)=rows
    
    # Now I can build the result..
    Result=list(ApproxCIBeta=ApproxCIBeta,Alpha=Alpha,CorrectedAlpha=CorrectedAlpha,Sigma2Hat=Sigma2Hat)
    return(Result)
    }



##### RESULT ANALYSIS #####

FixedPointPlot=function(x,y,SpaceBasisObj,TimeBasisObj,SolutionObj,lwd=1,NameLocation=NA,ylim=NA,N=100)
{
    # Plot of f(p,t) for a given p
    # Function will have only time variation
    #
    # Arguments:
    # X                 x coordinate of selected point
    # Y                 y coordinate of selected point
    # SPACEBASISOBJ     space basis object created with create.FEM.basis
    # TIMEBASISOBJ      time basis object created with create.bspline.basis (fda)
    # SOLUTIONOBJ       Result of ST.Smooth (estimated function)
    # LWD               lwd parameter of plot function
    # NAMELOCATION      Name of selected location (for main parameter of plot)
    # YLIM              ylim parameter of plot function
    # N                 number of time istants in range for the plot
    #                   (default 100)
    # Last modified on 15 April 2015 by Gabriele Mazza 
        
    # Check arguments
    if (class(SolutionObj)!="SolutionObj")
    {  
        stop('Not valid SolutionObj in FixedPointPlot')
    }
    if (class(SpaceBasisObj)!="basisfd")
    {  
        stop('Space basis object not valid')
    }
    if (SpaceBasisObj$type!="FEM")
    {  
        stop('Space basis are not FE')
    }
    if(class(TimeBasisObj)!="basisfd")
    {  
        stop('Time basis object not valid')
    }
    if(TimeBasisObj$type!="bspline")
    {  
        stop('Time basis are not bspline')
    }
    
    # For the title
    if(is.na(NameLocation))
    {
        title=paste("(",x,",",y,")",sep="")
    } else
    {
        title=NameLocation
    }
    # Create vectors for ST.Eval function
    xvec<-rep(x,N)
    yvec<-rep(y,N)
    Time<-seq(TimeBasisObj$rangeval[1],TimeBasisObj$rangeval[2],length.out=N)

    # Plot
    eval<-ST.Eval(xvec,yvec,Time,SpaceBasisObj,TimeBasisObj,SolutionObj)
    if(is.na(ylim)[1])
    {
        plot(Time,eval,type='l',xlab="Time",ylab=" ",main=paste("Time evolution in ",title,sep=""),lwd=2)
    } else
    {
        plot(Time,eval,type='l',xlab="Time",ylab=" ",main=paste("Time evolution in ",title,sep=""),lwd=lwd,ylim=ylim)
    }
}


FixedTimePlot=function(t,SpaceBasisObj,TimeBasisObj,SolutionObj,Nx=100,Ny=100)
{
    # Plot of f(p,t) for a given t
    # Function will be a bidimensional plot
    #
    # Arguments:
    # T                 selected time
    # SPACEBASISOBJ     space basis object created with create.FEM.basis
    # TIMEBASISOBJ      time basis object created with create.bspline.basis (fda)
    # SOLUTIONOBJ       Result of ST.Smooth (estimated function)
    # NX                number of x for spatial grid of plot (default 100)
    # NY                number of y for spatial grid of plot (default 100)
    #
    # Last modified on 15 April 2015 by Gabriele Mazza
    
    # Check arguments
    if (class(SolutionObj)!="SolutionObj")
    {  
        stop('Not valid SolutionObj in FixedPointPlot')
    }
    if (class(SpaceBasisObj)!="basisfd")
    {  
        stop('Space basis object not valid')
    }
    if (SpaceBasisObj$type!="FEM")
    {  
        stop('Space basis are not FE')
    }
    if(class(TimeBasisObj)!="basisfd")
    {  
        stop('Time basis object not valid')
    }
    if(TimeBasisObj$type!="bspline")
    {  
        stop('Time basis are not bspline')
    }
    if((t>TimeBasisObj$rangeval[2])|(t<TimeBasisObj$rangeval[1]))
    {  
        stop('Time is not in range values')
    }
    
    xmin = min(SpaceBasisObj$params$Points[,1])
    xmax = max(SpaceBasisObj$params$Points[,1])
    X    = matrix(seq(xmin, xmax, len=Nx),ncol=1)
    
    ymin = min(SpaceBasisObj$params$Points[,2])
    ymax = max(SpaceBasisObj$params$Points[,2])
    Y    = matrix(seq(ymin, ymax, len=Ny),ncol=1)    
    
    Xmat = X %*% matrix(1,nrow=1,ncol=Ny)
    Ymat = matrix(1,nrow=Nx,ncol=1) %*% t(Y)
    Xvec = NULL
    for (numc in 1:Ny)
    {
        Xvec=c(Xvec,Xmat[,numc])
    }
    Yvec = NULL
    for (numc in 1:Ny)
    {
        Yvec=c(Yvec,Ymat[,numc])
    }
    
    evalmat = ST.Eval(Xvec, Yvec, rep(t,length(Xvec)),SpaceBasisObj,TimeBasisObj,SolutionObj)
    evalmat = matrix(evalmat, nrow=Nx, ncol=Ny, byrow=F)
    
    image(X,Y,evalmat,col=heat.colors(100), xlab="", ylab="", asp=1,main=paste("Function at t=",t,sep=""))
    contour(X,Y,evalmat,add=T)
}
