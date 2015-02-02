# Creo le funzioni per l'analisi spazio-temporale


bwsubs=function(Q,b)
{
    N=dim(Q)[1]
    Sol=numeric(N)
    Sol[N]=b[N]/Q[N,N]
    for(i in (N-1):1)
    {
        temp=b[i]
        for(j in N:(i+1))
        {
            temp=temp-Sol[j]*Q[i,j]
        }
        Sol[i]=temp/Q[i,i]
    }
    return(as.matrix(Sol))
}

Create.Bspline.Time.Basis = function(TimePoints, TimeOrder, NumBasis=NA, PlotIt=FALSE)
{
    # CREATE.BsPLINE.TIME.BASIS 
    # computes Bspline time basis for spatio-temporal regression problem
    # 
    #
    # ARGUMENTS:
    #   TIMEPOINTS  times for which data are present
    #   TIMEORDER   one plus spline order
    #   PLOTIT      logical value, if TRUE, basis are plotted      
    #
    # RETURNS:
    # a list of three objects
    #   BASISOBJ    object containing the time basis (like fda implemented in
    #               fda package)
    #   TIMEPOINTS  times for which data are present
    #   TIMEORDER   one plus spline order
    
    
    # If NumBasis hasn't been set yet, default value is the number of points
    if(is.na(NumBasis))
    {
        NumBasis=length(TimePoints)
    }
    
    # Create the basis
    BasisObj = create.bspline.basis(c(min(TimePoints),max(TimePoints)), NumBasis, TimeOrder)
    
    TimeBasisObj<-list(BasisObj=BasisObj,TimePoints=TimePoints,TimeOrder=TimeOrder)   

    class(TimeBasisObj)<-"TimeBasisObj"
    
    # Plot of the basis?
    if(PlotIt==TRUE)
        plot(TimeBasisObj$BasisObj,main="Time Spline Basis",xlab="Time",ylab="Basis")
    
    return(TimeBasisObj)
}



Create.FEM.Space.Basis = function(Points, Triang, Internal, SpaceOrder, E=NULL, Dl=NULL)
{
    # CREATE.FEM.SPACE.BASIS
    # sets up a finite element basis for the analysis
    # of spatial data.  It requires a triangular mesh as input.
    # The triangular mesh is defined by a Delaunay triangulation of the domain 
    # of the basis functions, along with, optionally, a decomposition of the 
    # domain into subdomains (optional).
    # The finite elements used for functional data analysis are first or second
    # order Lagrangian elements. These are triangles covering a region,
    # and the basis system is piecewise polinomials (linear or quadratic). There is a basis
    # function associated with each node in the system.
    # When ORDER = 1 the basis system is piecewise linear and the nodes are the vertices of the triangles.
    # When ORDER = 2 the basis system is piecewise quadratic and the nodes are points that are etiher 
    # the vertices of the triangles or midpoints of edges of triangles.
    #
    # ARGUMENTS:
    #   Points      the NBASIS by 2 matrix of vertices of triangles
    #               containing the X- and Y-coordinates of the vertices.
    #               P may also be provided as a 2 by NBASIS matrix.
    #   Triang      The 4 by number of triangles matrix specifying triangles
    #               and their properties:
    #               Columns 1 to 3:   the indices in P of the vertices of each 
    #               triangle in counter-clockwise order.
    #               Column  4:        the subdomain number of the triangle
    #               T may also be provided with 3 columns, in which case the
    #               final column is set to ones.  It may also be provided as
    #               either a 3 or 4 by number of triangles matrix.
    #   INTERNAL    For every point, is TRUE if the point is internal (so data
    #               is available on the point), FALSE if it is on the boundary
    #   SPACEORDER  Order of elements, which may be either 1 or 2 (2 is default)
    #   E           The number of edges by 7 matrix defining the segments of the 
    #               boundary of the region which are also edges of the triangles
    #               adjacent to the boundary.  
    #               The values in matrix E as follows:
    #               Columns 1 and 2:  the indices of the vertices in matrix P of
    #               the starting and ending points of the edges
    #               Columns 3 and 4:  the starting and ending parameter values.  
    #               These are initialized to 0 and 1.
    #               Column  5:        the edge segment number 
    #               Columns 6 and 7:  the left and right-hand side subdomain numbers
    #               E may also be provided as a 7 by number of edges matrix.
    #   DL          Decomposed geometry object.  This defines the characteristics
    #               of the subdomains and the boundary segments.  It is optional.
    #  
    #  RETURNS:
    
    
    #  Check of inputs
    if(dim(Triang)[[2]] == 3)
    {
        # add a column of 1's
        Triang = cbind(Triang, numeric(dim(Triang)[[1]])+1)
    }
    
    if(dim(Triang)[[1]] == 3)
    {
        # add a row of 1's
        Triang = rbind(Triang, numeric(dim(Triang)[[1]])+1)
    }
    
    #  check dimensions of P, E and TRIANG and transpose if necessary
    
    if(is.null(E))
    {
        if((dim(Points)[[1]] != 2 && dim(Points)[[2]] != 2) || (dim(Triang)[[1]] != 4 && dim(Triang)[[2]] != 4))
        {
            stop('Dimensions of at least one of POINTS and TRIANG are not correct.')
        }
    } else
    {
        if((dim(Points)[[1]] != 2 && dim(Points)[[2]] != 2) || (dim(E)[[1]] != 7 && dim(E)[[2]] != 7) || (dim(Triang)[[1]] != 4 && dim(Triang)[[2]] != 4))
        {
            stop('Dimensions of at least one of POINTS, E and TRIANG are not correct.')
        }
    }
    
    
    if(dim(Points)[[2]] != 2 && dim(Points)[[1]] == 2)
    {
        Points = t(Points)
    }
    
    if(!is.null(E))
    {
        if(dim(E)[[2]] != 7 && dim(E)[[1]] == 7)
        {
            E = t(E)
        }
    }
    
    if (dim(Triang)[[2]] != 4 && dim(Triang)[[1]] == 4)
    {
        Triang = t(Triang)
    }
    
    type <- "FEM"
    
    # Argument rangeval is not needed for an FEM basis since domain
    # boundary is defined by the triangular mesh.
    rangeval = NULL  
    
    # The number of basis functions corresponds to the number of vertices
    # for order = 1, and to vertices plus edge midpoints for order = 2
    
    # set up the nodal points and corresponding mesh:  
    # this is p' and t(1:3,:)' in the order 1 case, but
    # includes mid-points in the order 2 case
    
    
    NodeStruct= makenodes(Points,Triang[,1:3], SpaceOrder)
    
    #  The params argument is a struct object
    params=NULL
    params$Points    = Points
    params$Internal  = Internal
    params$E         = E
    params$Triang    = Triang
    params$SpaceOrder= SpaceOrder
    params$Nodes     = NodeStruct$nodes
    params$NodeIndex = NodeStruct$nodeindex
    params$J         = NodeStruct$J
    params$Metric    = NodeStruct$metric
    
    if (!is.null(Dl))
    {
        params$Dl = Dl
    } 
    
    nbasis = dim(params$Nodes)[[1]]
    
    SpaceBasisObj = basisfd(type, rangeval, nbasis, params)
    
    return(SpaceBasisObj)
}




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



smooth.ST.fd = function(Data,SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT)
{
    # SMOOTH.ST.FD Compute a solution for a Spatial-Time Spline problem 
    #
    # Arguments:
    # FELSPLOBJ     a FELspline object.
    # LAMBDAS       a scalar smoothing parameter for space term.
    # LAMBDAT       a scalar smoothing parameter for time term.
    # DATA          a (n+1)-by-(m+1) set of noisy observations of the surface values.  
    #               DATA(:,1) indexes the points at which the 
    #               values in DATA(:,2:end) were observed.
    #               DATA(1,:) indexes the times at which the values of
    #               DATA(2:end,:) are observed
    #
    # Output:
    # FELSPLOBJ  ...  A FD object of the FEM type defined by the coefficient
    #                 vector resulting from smoothing
    # LAPLACEFD  ...  A FD object of the FEM type for the value of the 
    #                 Laplace operator if order == 2, or empty if order == 1
    #
    #
    # Last modified on 8 February 2011 by Laura Sangalli
    
    #  check arguments
    
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
    if (!is.numeric(LambdaS))
    {  
        stop('LAMBDAS is not numeric')
    }
    if (length(LambdaS) != 1)
    {   
        stop('LAMBDAS is not a scalar')
    }
    # Time
    if(class(TimeBasisObj)!="TimeBasisObj")
    {  
        stop('Time basis object not valid')
    }
    if(TimeBasisObj$BasisObj$type!="bspline")
    {  
        stop('Time basis are not bspline')
    }
    if (!is.numeric(LambdaT))
    {  
        stop('LAMBDAT is not numeric')
    }
    if (length(LambdaT) != 1)
    {   
        stop('LAMBDAT is not a scalar')
    }
    
    # Create Pi matrix
    t1<-Sys.time()
    Pi=MakePi(SpaceBasisObj,TimeBasisObj)
    t2<-Sys.time()
    DeltaT<-t2-t1
    print(paste("Pi Matrix: ",DeltaT," sec",sep=""))
    
    # Create S matrix
    t1<-Sys.time()
    S=MakeS(SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT,1)
    t2<-Sys.time()
    DeltaT<-t2-t1
    print(paste("S Matrix: ",DeltaT," sec",sep=""))
    
    # From data matrix, build z vector
    z<-NULL
    for(j in 1:dim(Data)[2])
    {
        z<-c(z,Data[,j])
    }
    z<-as.matrix(z)
    #Sol=solve(t(Pi)%*%Pi+S)%*%t(Pi)%*%z
    
    t1<-Sys.time()
    QR=qr((t(Pi)%*%Pi+S))
    Q=qr.Q(QR)
    R=qr.R(QR)
    t2<-Sys.time()
    DeltaT<-t2-t1
    print(paste("QR factorization: ",DeltaT," sec",sep=""))
    
    t1<-Sys.time()
    b=t(Q)%*%t(Pi)%*%z
    Sol<-bwsubs(R,b)
    t2<-Sys.time()
    DeltaT<-t2-t1
    print(paste("Solution: ",DeltaT," sec",sep=""))
    
    # Now that I've found the solution, i create an object with coefficients
    Timenbasis=TimeBasisObj$BasisObj$nbasis
    Spacenbasis=SpaceBasisObj$nbasis
    
    # Check for the sistem solution
    if(dim(Sol)[1]!=Timenbasis*Spacenbasis)
    {
        stop('Something wrong with solution dimensions')
    }
    
    # C is a matrix with rows containing time coefficients for a fixed space point,
    # columns containig space coefficients for a fixed time istant
    C<-NULL
    for(j in 1:Timenbasis)
    {
        C<-cbind(C,Sol[((j-1)*Spacenbasis+1):(j*Spacenbasis),1])
    }
    
    return(C)
    
    # Make FELspline object
    
#     felsplobj  = fd(u, basisobj)
#     laplacefd = fd(s, basisobj)  
#     
#     reslist=list(felsplobj=felsplobj,laplacefd=laplacefd)   
#     return(reslist)

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

order     = nodeStruct$SpaceOrder
nodes     = nodeStruct$Nodes
nodeindex = nodeStruct$NodeIndex
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

order     = nodeStruct$SpaceOrder
nodes     = nodeStruct$Nodes
nodeindex = nodeStruct$NodeIndex
Jvec      = nodeStruct$J
metric    = nodeStruct$Metric

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





eval.FEM.fd = function(X,Y,fdObj)
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

basisobj = fdObj$basis

#  get nodes and mesh

params      = basisobj$params
order       = params$SpaceOrder
nodes       = params$Nodes
nodeindex   = params$NodeIndex
p           = params$Points
t           = params$Triang
t           = t[,1:3]

Jvec      = params$J

#  get the coefficient matrix

coefmat = fdObj$coef

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


MakePi = function(SpaceBasisObj,TimeBasisObj)
{

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
    # Time
    if(class(TimeBasisObj)!="TimeBasisObj")
    {  
        stop('Time basis object not valid')
    }
    if(TimeBasisObj$BasisObj$type!="bspline")
    {  
        stop('Time basis are not bspline')
    }
    
    
    ### PHI MATRIX - SPACE ###
    
    # Now i build an fdobj for each basis function
    Phi<-NULL
    
    #Trovo i punti
    X<-SpaceBasisObj$params$Points[,1]
    Y<-SpaceBasisObj$params$Points[,2]
    
    #I want only points with data available
    X<-X[SpaceBasisObj$params$Internal]
    Y<-Y[SpaceBasisObj$params$Internal]
    
    #I build Phi as indentity matrix
    #     Phi<-matrix(data=0,nrow=length(X),ncol=length(SpaceBasisObj$params$Points[,1]))
    #     for(i in 1:length(X))
    #     {
    #         Phi[i,i]=1
    #     }
    #     
    
    for(i in 1:SpaceBasisObj$nbasis)
    {
        # basis function i
        coeff<-numeric(SpaceBasisObj$nbasis)
        coeff[i]<-1
        fdObj<-fd(coeff,SpaceBasisObj)
        
        eval<-eval.FEM.fd(X,Y,fdObj)
        
        Phi<-cbind(Phi,eval)
    }
    
    
    ### PSI MATRIX - TIME ###
    
    Psi=eval.basis(TimeBasisObj$TimePoints,TimeBasisObj$BasisObj)
    
    # I don't want basis names...
    dimnames(Psi)[[2]]<-NULL
    
    
    ### PI MATRIX - SPACE+TIME ###

    Pi=kronecker(Phi,Psi)
    return(Pi)
}




MakeS = function(SpaceBasisObj,TimeBasisObj,LambdaS,LambdaT,DerivativeOrder=1)
{
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
    # Time
    if(class(TimeBasisObj)!="TimeBasisObj")
    {  
        stop('Time basis object not valid')
    }
    if(TimeBasisObj$BasisObj$type!="bspline")
    {  
        stop('Time basis are not bspline')
    }

    
    ### SSPACE ###
    
    NodeStruct = SpaceBasisObj$params
    
    K0 = mass(NodeStruct)
    
    K1 = stiff1(NodeStruct)
    
    Sspace = K1%*%solve(K0)%*%K1
    
    
    ### STIME ###
    
    Stime=eval.penalty(TimeBasisObj$BasisObj, DerivativeOrder)
    
    
    ### S=SSPACE+STIME ###
    
    # S is created using kronaker product
    S=LambdaS*kronecker(diag(dim(Stime)[1]),Sspace) + LambdaT*kronecker(Stime,diag(dim(Sspace)[1]))
    
    return(S)
}



eval.ST.fd = function(X,Y,Time,C,SpaceBasisObj,TimeBasisObj)
{
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
    # Time
    if(class(TimeBasisObj)!="TimeBasisObj")
    {  
        stop('Time basis object not valid')
    }
    if(TimeBasisObj$BasisObj$type!="bspline")
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
    # Check Y
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
    # Check Time
    if (!is.numeric(Time))
    {
        stop('TIME is not a numerical array')
    } else if (length(Time)!=length(X))
    {
        stop('TIME is not the same length as X and Y')
    }
    # TIME must be a vector, not a matrix
    
    # I Want to use the old function eval.FEM.fd
    N = length(X)
    
    # I need Phi matrix in used time points
    Psi=eval.basis(Time,TimeBasisObj$BasisObj)
    if ((dim(Psi)[2])!=(dim(C)[2]))
    {  
        stop('Something wrong with dimensions')
    }
    
    Result<-numeric(N)
    
    for (j in 1:(dim(C)[2]))
    {
        fdObj=fd(C[,j],SpaceBasisObj)
        temp<-eval.FEM.fd(X,Y,fdObj)
        Result=Result+temp*Psi[,j]
    }
    
    return(Result)
}