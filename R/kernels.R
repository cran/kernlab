# Define the kernel objects,
# functions with an additional slot for the kernel parameter list.

setClass("kernel",representation("function",kpar="list"))


rbfdot<- function(sigma=1)
  {

    rval <- function(x,y=NULL)
    {
       if(!is.vector(x)) stop("x must be a vector")
      if(!is.vector(y)&&!is.null(y)) stop("y must a vector")
      if (is.vector(x) && is.null(y)){
        return(1)
      }
      if (is.vector(x) && is.vector(y)){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")
        return(exp(sigma*(2*crossprod(x,y) - crossprod(x) - crossprod(y))))  
        # sigma/2 or sigma ??
      }
    }
     return(new("rbfkernel",.Data=rval,kpar=list(sigma=sigma)))
  }
setClass("rbfkernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

laplacedot<- function(sigma=1)
  {
    rval <- function(x,y=NULL)
    {
       if(!is.vector(x)) stop("x must be a vector")
      if(!is.vector(y)&&!is.null(y)) stop("y must a vector")
      if (is.vector(x) && is.null(y)){
        return(1)
      }
      if (is.vector(x) && is.vector(y)){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")
        return(exp(-sigma*sqrt(-(round(2*crossprod(x,y) - crossprod(x) - crossprod(y),9)))))
      }
    }
     return(new("laplacekernel",.Data=rval,kpar=list(sigma=sigma)))
  }

setClass("laplacekernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

besseldot<- function(sigma = 1, order = 1, degree = 1)
  {
    rval <- function(x,y=NULL)
    {
      if(!is.vector(x)) stop("x must be a vector")
      if(!is.vector(y)&&!is.null(y)) stop("y must a vector")
      if (is.vector(x) && is.null(y)){
        return(1)
      }
      if (is.vector(x) && is.vector(y)){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")
        lim <- 1/(gamma(order+1)*2^(order))
        bkt <- sigma*sqrt(-(2*crossprod(x,y) - crossprod(x) - crossprod(y)))
        if(bkt < 10e-5)
          res <- lim
        else
          res <- besselJ(bkt,order)*(bkt^(-order))
        return((res/lim)^degree)
      }
    }
     return(new("besselkernel",.Data=rval,kpar=list(sigma=sigma ,order = order ,degree = degree)))
  }

setClass("besselkernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

anovadot<- function(sigma = 1, degree = 1)
  {
    rval <- function(x,y=NULL)
    {
      if(!is.vector(x)) stop("x must be a vector")
      if(!is.vector(y)&&!is.null(y)) stop("y must a vector")
      if (is.vector(x) && is.null(y)){
        return(1)
      }
      if (is.vector(x) && is.vector(y)){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")

          res <- sum(exp(- sigma * (x - y)^2))
        return((res)^degree)
      }
    }
     return(new("anovakernel",.Data=rval,kpar=list(sigma=sigma ,degree = degree)))
  }

setClass("anovakernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))


splinedot<- function()
  {
    rval <- function(x,y=NULL)
    {
      if(!is.vector(x)) stop("x must be a vector")
      if(!is.vector(y)&&!is.null(y)) stop("y must a vector")
      if (is.vector(x) && is.null(y)){
        return(1)
      }
      if (is.vector(x) && is.vector(y)){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")
        minv <- pmin(x,y)
        res <- 1 + x*y*(1+minv) - ((x+y)/2)*minv^2 + (minv^3)/3
          fres <- prod(res)
        return(fres)
      }
    }
     return(new("splinekernel",.Data=rval,kpar=list()))
  }

setClass("anovakernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))



tanhdot <- function(scale = 1, offset = 1)
{
  rval<- function(x, y = NULL)
    {
      if(!is.vector(x)) stop("x must be a  vector")
      if(!is.vector(y)&&!is.null(y)) stop("y must be a vector")
      if (is.vector(x) && is.null(y)){
        tanh(scale*crossprod(x)+offset)
      }
      if (is.vector(x) && is.vector(y)){
        if (!length(x)==length(y))
          stop("number of dimension must be the same on both data points")
        tanh(scale*crossprod(x,y)+offset)
      }
    }
  return(new("tanhkernel",.Data=rval,kpar=list(scale=scale,offset=offset)))
}
setClass("tanhkernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

setClass("polykernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

polydot <- function(degree = 1, scale = 1, offset = 1)
{
  rval<- function(x, y = NULL)
    {
      if(!is.vector(x)) stop("x must be a vector")
      if(!is.vector(y)&&!is.null(y)) stop("y must be a vector")
      if (is.vector(x) && is.null(y)){
        (scale*crossprod(x)+offset)^degree
        }

      if (is.vector(x) && is.vector(y)){
        if (!length(x)==length(y))
            stop("number of dimension must be the same on both data points")
        (scale*crossprod(x,y)+offset)^degree
          }

      }
  return(new("polykernel",.Data=rval,kpar=list(degree=degree,scale=scale,offset=offset)))
}

setClass("vanillakernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

vanilladot <- function( )
{
  rval<- function(x, y = NULL)
    {
      if(!is.vector(x)) stop("x must be a vector")
      if(!is.vector(y)&&!is.null(y)) stop("y must be a vector")
      if (is.vector(x) && is.null(y)){
        crossprod(x)
        }

      if (is.vector(x) && is.vector(y)){
        if (!length(x)==length(y))
            stop("number of dimension must be the same on both data points")
        crossprod(x,y)
          }

      }
  return(new("vanillakernel",.Data=rval,kpar=list()))
}

#create accesor function as in "S4 Classses in 15 pages more or less", well..

if (!isGeneric("kpar")){
  if (is.function("kpar"))
    fun <- kpar
  else fun <- function(object) standardGeneric("kpar")
  setGeneric("kpar",fun)
}

setMethod("kpar","kernel", function(object) object@kpar)




# Functions that return usefull kernel calculations (kernel matrix etc.)


kernelMatrix <- function(kernel, x, y=NULL)
{
  if(is.vector(x))
    x <- as.matrix(x)
  if(is.vector(y))
    y <- as.matrix(y)
  if(!is.matrix(x)) stop("x must be a matrix")
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix")
  n <- nrow(x)
  res1 <- matrix(rep(0,n*n), ncol = n)
  if(is.null(y)){
    for(i in 1:n) {
      for(j in 1:n) {
        res1[i,j] <- kernel(x[i,],x[j,])
      }
    }
  }
  if (is.matrix(y)){
    m<-dim(y)[1]
    res1 <- matrix(0,dim(x)[1],dim(y)[1])
    for(i in 1:n) {
      for(j in 1:m) {
        res1[i,j] <- kernel(x[i,],y[j,])
      }
    }
  }
  return(res1)
}

setGeneric("kernelMatrix",function(kernel, x, y = NULL) standardGeneric("kernelMatrix"))



kernelMatrix.rbfkernel <- function(kernel, x, y = NULL)
{
  if(is.vector(x))
    x <- as.matrix(x)
  if(is.vector(y))
    y <- as.matrix(y)
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix or a vector")
  sigma = kpar(kernel)$sigma
  n <- dim(x)[1]
  dota <- rowSums(x*x)/2
  if (is.matrix(x) && is.null(y)){
    res <- crossprod(t(x))
    for (i in 1:n)
      res[i,]<- exp(2*sigma*(res[i,] - dota - rep(dota[i],n)))
    return(res)
  }
  if (is.matrix(x) && is.matrix(y)){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    m <- dim(y)[1]
    dotb <- rowSums(y*y)/2
    res <- x%*%t(y)
    for( i in 1:m)
      res[,i]<- exp(2*sigma*(res[,i] - dota - rep(dotb[i],n)))
    return(res)
  }
}
setMethod("kernelMatrix",signature(kernel="rbfkernel",x="matrix"),kernelMatrix.rbfkernel)

kernelMatrix.laplacekernel <- function(kernel, x, y = NULL)
{
  if(is.vector(x))
    x <- as.matrix(x)
  if(is.vector(y))
    y <- as.matrix(y)
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix or a vector")
  sigma = kpar(kernel)$sigma
  n <- dim(x)[1]
  dota <- rowSums(x*x)/2
  if (is.matrix(x) && is.null(y)){
    res <- crossprod(t(x))
    for (i in 1:n)
      res[i,]<- exp(-sigma*sqrt(round(-2*(res[i,] - dota - rep(dota[i],n)),9)))
    return(res)
  }
  if (is.matrix(x) && is.matrix(y)){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    m <- dim(y)[1]
    dotb <- rowSums(y*y)/2
    res <- x%*%t(y)
    for( i in 1:m)
      res[,i]<- exp(-sigma*sqrt(round(-2*(res[,i] - dota - rep(dotb[i],n)),9)))
    return(res)
  }
}
setMethod("kernelMatrix",signature(kernel="laplacekernel",x="matrix"),kernelMatrix.laplacekernel)

kernelMatrix.besselkernel <- function(kernel, x, y = NULL)
{
  if(is.vector(x))
    x <- as.matrix(x)
  if(is.vector(y))
    y <- as.matrix(y)
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix or a vector")
  sigma = kpar(kernel)$sigma
  nu = kpar(kernel)$order
  ni = kpar(kernel)$degree
  n <- dim(x)[1]
  lim <- 1/(gamma(nu+1)*2^(nu))
  dota <- rowSums(x*x)/2
  if (is.matrix(x) && is.null(y)){
    res <- crossprod(t(x))
    for (i in 1:n){
      xx <- sigma*sqrt(round(-2*(res[i,] - dota - rep(dota[i],n)),9))
      res[i,] <- besselJ(xx,nu)*(xx^(-nu))
      res[i,which(xx<10e-5)] <- lim
    }
    return((res/lim)^ni)
  }
  if (is.matrix(x) && is.matrix(y)){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    m <- dim(y)[1]
    dotb <- rowSums(y*y)/2
    res <- x%*%t(y)
    for( i in 1:m){
      xx <- sigma*sqrt(round(-2*(res[,i] - dota - rep(dotb[i],n)),9))
      res[,i] <- besselJ(xx,nu)*(xx^(-nu))
      res[which(xx<10e-5),i] <- lim
    }     
    return((res/lim)^ni)
  }
}
setMethod("kernelMatrix",signature(kernel="besselkernel",x="matrix"),kernelMatrix.besselkernel)


kernelMatrix.anovakernel <- function(kernel, x, y = NULL)
{
  if(is.vector(x))
    x <- as.matrix(x)
  if(is.vector(y))
    y <- as.matrix(y)
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix or a vector")
  sigma = kpar(kernel)$sigma
  degree = kpar(kernel)$degree
  n <- dim(x)[1]
  if (is.matrix(x) && is.null(y)){
    a <- matrix(0,  dim(x)[2], n)
    res <- matrix(0, n ,n)
    for (i in 1:n)
      {
        a[rep(TRUE,dim(x)[2]), rep(TRUE,n)] <- x[i,]
        res[i,]<- colSums(exp( - sigma*(a - t(x))^2))^degree
      }
    return(res)
  }
  if (is.matrix(x) && is.matrix(y)){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    
    m <- dim(y)[1]
    b <- matrix(0, dim(x)[2],m)
    res <- matrix(0, dim(x)[1],m)
    for( i in 1:n)
      {
        b[rep(TRUE,dim(x)[2]), rep(TRUE,m)] <- x[i,]
        res[i,]<- colSums(exp( - sigma*(b - t(y))^2))^degree
      }
    return(res)
  }
}
setMethod("kernelMatrix",signature(kernel="anovakernel",x="matrix"),kernelMatrix.anovakernel)


kernelMatrix.polykernel <- function(kernel, x, y = NULL)
{
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix")
  scale = kpar(kernel)$scale
  offset = kpar(kernel)$offset
  degree = kpar(kernel)$degree
  if (is.matrix(x) && is.null(y))
    {
      res <- (scale*crossprod(t(x))+offset)^degree
      return(res)
    }
      if (is.matrix(x) && is.matrix(y)){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    res <- (scale*crossprod(t(x),t(y)) + offset)^degree
    return(res)
  }
}
setMethod("kernelMatrix",signature(kernel="polykernel",x="matrix"),kernelMatrix.polykernel)

kernelMatrix.vanilla <- function(kernel, x, y = NULL)
{
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix")
  if (is.matrix(x) && is.null(y)){
    res <- crossprod(t(x))
    return(res)
  }
  if (is.matrix(x) && is.matrix(y)){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    res <- crossprod(t(x),t(y))
    return(res)
  }
}
setMethod("kernelMatrix",signature(kernel="vanillakernel",x="matrix"),kernelMatrix.vanilla)

kernelMatrix.tanhkernel <- function(kernel, x, y = NULL)
{
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix")
  if (is.matrix(x) && is.null(y)){
    scale = kpar(kernel)$scale
    offset = kpar(kernel)$offset
    res <- tanh(scale*crossprod(t(x)) + offset)
    return(res)
  }
  if (is.matrix(x) && is.matrix(y)){
    if (!(dim(x)[2]==dim(y)[2]))
      stop("matrixes must have the same number of columns")
    res <- tanh(scale*crossprod(t(x),t(y)) + offset)
    return(res)
  }
}
setMethod("kernelMatrix",signature(kernel="tanhkernel",x="matrix"),kernelMatrix.tanhkernel)

# Function computing <x,x'> * z (<x,x'> %*% z)


kernelMult <- function(kernel, x, y=NULL, z, blocksize = 128)
{
#  if(is.function(kernel)) ker <- deparse(substitute(kernel))
#      kernel <- do.call(kernel, kpar)

  if(!is.matrix(x)) stop("x must be a matrix")
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix")
  if(!is.matrix(z)&&!is.vector(z)) stop("z must ba a matrix or a vector")
  n <- nrow(x)

  if(is.null(y))
    {
      if(is.vector(z))
        {if(is.null(y)&&!length(z)==n) stop("vector z length must be equal to x rows")
         z<-matrix(z,n,1)}
      if(is.null(y)&&!dim(z)[1]==n)
        stop("z must have the length equal to x colums")
      res1 <- matrix(rep(0,n*n), ncol = n)     
      
      for(i in 1:n)
        {
          for(j in 1:n)
            {
              res1[i,j] <- kernel(x[i,],x[j,])
      }
    }
  }
  if (is.matrix(y))
    {
      
      m<-dim(y)[1]
      if(is.vector(z))
       z <- as.matrix(z)

      if(!dim(z)[1] == m) stop("z has wrong dimension")
      res1 <- matrix(rep.int(0,m*n),ncol=m)
      for(i in 1:n)
        {
          for(j in 1:m)
            {
              res1[i,j] <- kernel(x[i,],y[j,])
            }
        }
    }
  return(res1%*%z)
}

setGeneric("kernelMult",function(kernel, x, y = NULL, z, blocksize = 256) standardGeneric("kernelMult"))

kernelMult.rbfkernel <- function(kernel, x, y=NULL, z, blocksize = 256)
{
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix")
  if(!is.matrix(z)&&!is.vector(z)) stop("z must be a matrix or a vector")
  sigma <- kpar(kernel)$sigma
  n <- dim(x)[1]
  m <- dim(x)[2]
  nblocks <- floor(n/blocksize)
  lowerl <- 1
  upperl <- 0
  dota <- as.matrix(rowSums(x^2))

  if (is.null(y))
    {
      if(is.vector(z))
        {
          if(!length(z) == n) stop("vector z length must be equal to x rows")
          z <- matrix(z,n,1)
        }
      if(!dim(z)[1]==n)
        stop("z rows must equal x rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      if(nblocks > 0)
        {
          dotab <- rep(1,blocksize)%*%t(dota)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            res[lowerl:upperl,] <- exp(sigma*(2*x[lowerl:upperl,]%*%t(x) - dotab - dota[lowerl:upperl]%*%t(rep.int(1,n))))%*%z
            lowerl <- upperl + 1
          
          }
      }
      if(lowerl <= n)
        res[lowerl:n,] <- exp(sigma*(2*x[lowerl:n,]%*%t(x) - rep.int(1,n+1-lowerl)%*%t(dota) - dota[lowerl:n]%*%t(rep.int(1,n))))%*%z

    }
  if(is.matrix(y))
    {
      n2 <- dim(y)[1]
      if(is.vector(z))
        {
          if(!length(z) == n2) stop("vector z length must be equal to y rows")
          z <- matrix(z,n2,1)
        }
      if(!dim(z)[1]==n2)
        stop("z length must equal y rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      dotb <- as.matrix(rowSums(y*y))
      
       if(nblocks > 0)
         {
           dotbb <-  rep(1,blocksize)%*%t(dotb)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            res[lowerl:upperl,] <- exp(sigma*(2*x[lowerl:upperl,]%*%t(y) - dotbb - dota[lowerl:upperl]%*%t(rep.int(1,n2))))%*%z
            lowerl <- upperl + 1
          }
      }
      if(lowerl <= n)
        res[lowerl:n,] <- exp(sigma*(2*x[lowerl:n,]%*%t(y) - rep.int(1,n+1-lowerl)%*%t(dotb) - dota[lowerl:n]%*%t(rep.int(1,n2))))%*%z
    }
  return(res)
}  
setMethod("kernelMult",signature(kernel="rbfkernel", x="matrix"),kernelMult.rbfkernel)


kernelMult.laplacekernel <- function(kernel, x, y=NULL, z, blocksize = 256)
{
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix")
  if(!is.matrix(z)&&!is.vector(z)) stop("z must be a matrix or a vector")
  sigma <- kpar(kernel)$sigma
  n <- dim(x)[1]
  m <- dim(x)[2]
  nblocks <- floor(n/blocksize)
  lowerl <- 1
  upperl <- 0
  dota <- as.matrix(rowSums(x^2))

  if (is.null(y))
    {
      if(is.vector(z))
        {
          if(!length(z) == n) stop("vector z length must be equal to x rows")
          z <- matrix(z,n,1)
        }
      if(!dim(z)[1]==n)
        stop("z rows must equal x rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      if(nblocks > 0)
        {
          dotab <- rep(1,blocksize)%*%t(dota)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            res[lowerl:upperl,] <- exp(-sigma*sqrt(-round(2*x[lowerl:upperl,]%*%t(x) - dotab - dota[lowerl:upperl]%*%t(rep.int(1,n)),9)))%*%z
            lowerl <- upperl + 1
          
          }
      }
      if(lowerl <= n)
        res[lowerl:n,] <- exp(-sigma*sqrt(-round(2*x[lowerl:n,]%*%t(x) - rep.int(1,n+1-lowerl)%*%t(dota) - dota[lowerl:n]%*%t(rep.int(1,n)),9)))%*%z

    }
  if(is.matrix(y))
    {
      n2 <- dim(y)[1]
      if(is.vector(z))
        {
          if(!length(z) == n2) stop("vector z length must be equal to y rows")
          z <- matrix(z,n2,1)
        }
      if(!dim(z)[1]==n2)
        stop("z length must equal y rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      dotb <- as.matrix(rowSums(y*y))
      
       if(nblocks > 0)
         {
           dotbb <-  rep(1,blocksize)%*%t(dotb)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            res[lowerl:upperl,] <- exp(-sigma*sqrt(-round(2*x[lowerl:upperl,]%*%t(y) - dotbb - dota[lowerl:upperl]%*%t(rep.int(1,n2)),9)))%*%z
            lowerl <- upperl + 1
          }
      }
      if(lowerl <= n)
        res[lowerl:n,] <- exp(-sigma*sqrt(-round(2*x[lowerl:n,]%*%t(y) - rep.int(1,n+1-lowerl)%*%t(dotb) - dota[lowerl:n]%*%t(rep.int(1,n2)),9)))%*%z
    }
  return(res)
}  
setMethod("kernelMult",signature(kernel="laplacekernel", x="matrix"),kernelMult.laplacekernel)



kernelMult.besselkernel <- function(kernel, x, y=NULL, z, blocksize = 256)
{
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix")
  if(!is.matrix(z)&&!is.vector(z)) stop("z must be a matrix or a vector")
  sigma <- kpar(kernel)$sigma
  nu <- kpar(kernel)$order
  ni <- kpar(kernel)$degree
  n <- dim(x)[1]
  m <- dim(x)[2]
  nblocks <- floor(n/blocksize)
  lowerl <- 1
  upperl <- 0
  lim <- 1/(gamma(nu+1)*2^(nu))
  dota <- as.matrix(rowSums(x^2))

  if (is.null(y))
    {
      if(is.vector(z))
        {
          if(!length(z) == n) stop("vector z length must be equal to x rows")
          z <- matrix(z,n,1)
        }
      if(!dim(z)[1]==n)
        stop("z rows must equal x rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      if(nblocks > 0)
        {
          dotab <- rep(1,blocksize)%*%t(dota)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            
            xx <- sigma*sqrt(-round(2*x[lowerl:upperl,]%*%t(x) - dotab - dota[lowerl:upperl]%*%t(rep.int(1,n)),9))
            res1 <- besselJ(xx,nu)*(xx^(-nu))
            res1[which(xx<10e-5)] <- lim
           
            res[lowerl:upperl,] <- ((res1/lim)^ni)%*%z
            lowerl <- upperl + 1
          }
      }
      if(lowerl <= n)
        {
          xx <- sigma*sqrt(-round(2*x[lowerl:n,]%*%t(x) - rep.int(1,n+1-lowerl)%*%t(dota) - dota[lowerl:n]%*%t(rep.int(1,n)),9))
          res1 <- besselJ(xx,nu)*(xx^(-nu))
          res1[which(xx<10e-5)] <- lim
          res[lowerl:n,] <- ((res1/lim)^ni)%*%z
      }
    }
  if(is.matrix(y))
    {
      n2 <- dim(y)[1]
      if(is.vector(z))
        {
          if(!length(z) == n2) stop("vector z length must be equal to y rows")
          z <- matrix(z,n2,1)
        }
      if(!dim(z)[1]==n2)
        stop("z length must equal y rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      dotb <- as.matrix(rowSums(y*y))
      
       if(nblocks > 0)
         {
           dotbb <-  rep(1,blocksize)%*%t(dotb)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            xx <- sigma*sqrt(-round(2*x[lowerl:upperl,]%*%t(y) - dotbb - dota[lowerl:upperl]%*%t(rep.int(1,n2)),9))
            res1 <- besselJ(xx,nu)*(xx^(-nu))
            res1[which(xx < 10e-5)] <- lim

            res[lowerl:upperl,] <- ((res1/lim)^ni)%*%z
            lowerl <- upperl + 1
          }
      }
      if(lowerl <= n)
      {
        xx <- sigma*sqrt(-round(2*x[lowerl:n,]%*%t(y) - rep.int(1,n+1-lowerl)%*%t(dotb) - dota[lowerl:n]%*%t(rep.int(1,n2)),9))
        res1 <- besselJ(xx,nu)*(xx^(-nu))
        res1[which(xx < 10e-5)] <- lim
        res[lowerl:n,] <- ((res1/lim)^ni)%*%z 
      }
    }
  return(res)
}  
setMethod("kernelMult",signature(kernel="besselkernel", x="matrix"),kernelMult.besselkernel)

kernelMult.anovakernel <- function(kernel, x, y=NULL, z, blocksize = 256)
{
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix")
  if(!is.matrix(z)&&!is.vector(z)) stop("z must be a matrix or a vector")
  sigma <- kpar(kernel)$sigma
  degree <- kpar(kernel)$degree
  n <- dim(x)[1]
  m <- dim(x)[2]
  nblocks <- floor(n/blocksize)
  lowerl <- 1
  upperl <- 0
 

  if (is.null(y))
    {
      if(is.vector(z))
        {
          if(!length(z) == n) stop("vector z length must be equal to x rows")
          z <- matrix(z,n,1)
        }
      if(!dim(z)[1]==n)
        stop("z rows must equal x rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      
      if(nblocks > 0)
        {
          a <- matrix(0,m,blocksize)
          re <- matrix(0, n, blocksize)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            for(j in 1:n)
              {
                a[rep(TRUE,m),rep(TRUE,blocksize)] <- x[j,]
                re[j,] <-  colSums(exp( - sigma*(a - t(x[lowerl:upperl,]))^2))^degree
              }
            res[lowerl:upperl,] <- t(re)%*%z
            lowerl <- upperl + 1
          
          }
        }
      if(lowerl <= n){
        a <- matrix(0,m,n-lowerl+1)
        re <- matrix(0,n,n-lowerl+1)
        for(j in 1:n)
          {
            a[rep(TRUE,m),rep(TRUE,n-lowerl+1)] <- x[j,]
            re[j,] <-  colSums(exp( - sigma*(a - t(x[lowerl:n,,drop=FALSE]))^2))^degree
          }
        res[lowerl:n,] <- t(re)%*%z
      }
    }
  if(is.matrix(y))
    {
      n2 <- dim(y)[1]
      nblocks <- floor(n2/blocksize)
      if(is.vector(z))
        {
          if(!length(z) == n2) stop("vector z length must be equal to y rows")
          z <- matrix(z,n2,1)
        }
      if(!dim(z)[1]==n2)
        stop("z length must equal y rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])

      if(nblocks > 0)
        {
          b <- matrix(0, m, blocksize)
          re <- matrix(0, n, blocksize)
          for(i in 1:nblocks)
            {
              upperl = upperl + blocksize
              for(j in 1:n)
                {
                  b[rep(TRUE,dim(x)[2]), rep(TRUE,blocksize)] <- x[j,]
                  re[j,]<- colSums(exp( - sigma*(b - t(y[lowerl:upperl,]))^2)^degree)
                }
              res[,1] <- res[,1] + re %*%z[lowerl:upperl,]
              lowerl <- upperl + 1
            }
        }
      if(lowerl <= n)
        {
          b <- matrix(0, dim(x)[2], n2-lowerl+1)
          re <- matrix(0, n, n2-lowerl+1)
          for( i in 1:n)
            {
              b[rep(TRUE,dim(x)[2]),rep(TRUE,n2-lowerl+1)] <- x[i,]
              re[i,]<- colSums(exp( - sigma*(b - t(y[lowerl:n2,,drop=FALSE]))^2)^degree)
            }
       
          res[,1] <- res[,1] + re%*%z[lowerl:n2]
        }
    }
  return(res)
}  
setMethod("kernelMult",signature(kernel="anovakernel", x="matrix"),kernelMult.anovakernel)



kernelMult.polykernel <- function(kernel, x, y=NULL, z, blocksize = 256)
{
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix")
  if(!is.matrix(z)&&!is.vector(z)) stop("z must be a matrix or a vector")
  degree <- kpar(kernel)$degree
  scale <- kpar(kernel)$scale
  offset <- kpar(kernel)$offset
  n <- dim(x)[1]
  m <- dim(x)[2]
  nblocks <- floor(n/blocksize)
  lowerl <- 1
  upperl <- 0
  if (is.null(y))
    {
      if(is.vector(z))
        {
          if(!length(z) == n) stop("vector z length must be equal to x rows")
          z <- matrix(z,n,1)
        }
      if(!dim(z)[1]==n)
        stop("z rows must equal x rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      if(nblocks > 0)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            res[lowerl:upperl,] <- ((scale*x[lowerl:upperl,]%*%t(x) + offset)^degree) %*% z
            lowerl <- upperl + 1
          }
      if(lowerl <= n)
        res[lowerl:n,] <- ((scale*x[lowerl:n,]%*%t(x) +offset)^degree)%*%z
    }
  if(is.matrix(y))
    {
      n2 <- dim(y)[1]
      if(is.vector(z))
        {
          if(!length(z) == n2) stop("vector z length must be equal to y rows")
          z <- matrix(z,n2,1)
        }
      if(!dim(z)[1]==n2)
        stop("z length must equal y rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
 
       if(nblocks > 0)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            res[lowerl:upperl,] <- ((scale*x[lowerl:upperl,]%*%t(y) + offset)^degree)%*%z
            lowerl <- upperl + 1
          }
      if(lowerl <= n)
        res[lowerl:n,] <- ((scale*x[lowerl:n,]%*%t(y) + offset)^degree)%*%z
    }
  return(res)
} 
setMethod("kernelMult",signature(kernel="polykernel", x="matrix"),kernelMult.polykernel)


kernelMult.tanhkernel <- function(kernel, x, y=NULL, z, blocksize = 256)
{
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix")
  if(!is.matrix(z)&&!is.vector(z)) stop("z must be a matrix or a vector")
  scale <- kpar(kernel)$scale
  offset <- kpar(kernel)$offset
  n <- dim(x)[1]
  m <- dim(x)[2]
  nblocks <- floor(n/blocksize)
  lowerl <- 1
  upperl <- 0
 

  if (is.null(y))
    {
      if(is.vector(z))
        {
          if(!length(z) == n) stop("vector z length must be equal to x rows")
          z <- matrix(z,n,1)
        }
      if(!dim(z)[1]==n)
        stop("z rows must equal x rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      if(nblocks > 0)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            res[lowerl:upperl,] <- tanh(scale*x[lowerl:upperl,]%*%t(x) + offset) %*% z
            lowerl <- upperl + 1
          }
      if(lowerl <= n)
        res[lowerl:n,] <- tanh(scale*x[lowerl:n,]%*%t(x) +offset)%*%z
    }
  if(is.matrix(y))
    {
      n2 <- dim(y)[1]
      if(is.vector(z))
        {
          if(!length(z) == n2) stop("vector z length must be equal to y rows")
          z <- matrix(z,n2,1)
        }
      if(!dim(z)[1]==n2)
        stop("z length must equal y rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
 
       if(nblocks > 0)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            res[lowerl:upperl,] <- tanh(scale*x[lowerl:upperl,]%*%t(y) + offset)%*%z
            lowerl <- upperl + 1
          }
      if(lowerl <= n)
        res[lowerl:n,] <- tanh(scale*x[lowerl:n,]%*%t(y) + offset)%*%z
    }
  return(res)
} 
setMethod("kernelMult",signature(kernel="tanhkernel", x="matrix"),kernelMult.tanhkernel)


kernelMult.vanillakernel <- function(kernel, x, y=NULL, z, blocksize = 256)
{
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix")
  if(!is.matrix(z)&&!is.vector(z)) stop("z must be a matrix or a vector")
  n <- dim(x)[1]
  m <- dim(x)[2]
  nblocks <- floor(n/blocksize)
  lowerl <- 1
  upperl <- 0
 

  if (is.null(y))
    {
      if(is.vector(z))
        {
          if(!length(z) == n) stop("vector z length must be equal to x rows")
          z <- matrix(z,n,1)
        }
      if(!dim(z)[1]==n)
        stop("z rows must equal x rows")
      res <- t(crossprod(crossprod(x,z),t(x)))
    }
  if(is.matrix(y))
    {
      n2 <- dim(y)[1]
      if(is.vector(z))
        {
          if(!length(z) == n2) stop("vector z length must be equal to y rows")
          z <- matrix(z,n2,1)
        }
      if(!dim(z)[1]==n2)
        stop("z length must equal y rows")
      res <- t(crossprod(crossprod(y,z),t(x)))
    }
  return(res)
} 
setMethod("kernelMult",signature(kernel="vanillakernel", x="matrix"),kernelMult.vanillakernel)

# kernelPol returns the scalar product of x y componentwise with polarities
# of z and k 

kernelPol <- function(kernel, x, y=NULL, z, k=NULL)
{
  if(!is.matrix(x)) stop("x must be a matrix")
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix")
  if(!is.matrix(z)&&!is.vector(z)) stop("z must ba a matrix or a vector")
  n <- nrow(x)
   if(is.vector(z))
      {
        if(!length(z)==n) stop("vector z length must be equal to x rows")
        z<-matrix(z,n,1)
      }
  if(!dim(z)[1]==n)
    stop("z must have the length equal to x colums")
  res1 <- matrix(rep(0,n*n), ncol = n)
  if (is.null(y))
    {
      for(i in 1:n)
        {
          for(j in 1:n)
            {
              res1[i,j] <- kernel(x[i,],x[j,])*z[j]*z[i]
      }
    }
  }
  if (is.matrix(x) && is.matrix(y)){
    m <- dim(y)[1]
    if(is.null(k)) stop("k not specified!")
    if(is.vector(k))
      {
        if(!length(k)==m) stop("vector k length must be equal to x rows")
        k<-as.matrix(k,n,1)
      }
    if(!dim(x)[2]==dim(y)[2])
      stop("matrixes must have the same number of columns")
    if(!dim(z)[2]==dim(k)[2])
      stop("z and k vectors must have the same number of columns")
    if(!dim(x)[1]==dim(z)[1])
      stop("z and x must have the same number of rows")
    if(!dim(y)[1]==dim(k)[1])
      stop("y and k must have the same number of rows")
    res1 <- matrix(0,dim(x)[1],dim(y)[1])
    for(i in 1:n)
      {
        for(j in 1:m)
          {
            res1[i,j] <- kernel(x[i,],y[j,])*z[i]*k[j]
          }
      }
  }
  return(res1)
}

setGeneric("kernelPol", function(kernel, x, y=NULL, z, k = NULL) standardGeneric("kernelPol"))


kernelPol.rbfkernel <- function(kernel, x, y=NULL, z, k=NULL)
{
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix or NULL")
  if(!is.matrix(z)&&!is.vector(z)) stop("z must be a matrix or a vector")
  if(!is.matrix(k)&&!is.vector(k)&&!is.null(k)) stop("k must be a matrix or a vector")
  sigma <- kpar(kernel)$sigma
  n <- dim(x)[1]
  dota <- rowSums(x*x)/2
  if(is.vector(z))
    {
      if(!length(z)==n) stop("vector z length must be equal to x rows")
      z<-matrix(z,n,1)
    }
  if(!dim(z)[1]==n)
    stop("z must have the length equal to x colums")
  if (is.null(y))
    {
      if(is.matrix(z)&&!dim(z)[1]==n)
       stop("z must have size equal to x colums")
      res <- crossprod(t(x))
      for (i in 1:n)
        res[i,] <- z[i,]*(exp(2*sigma*(res[i,] - dota - rep(dota[i],n)))*z)
      return(res)
    }
  if (is.matrix(y))
    {
      if(is.null(k)) stop("k not specified!")
      m <- dim(y)[1]
      if(!dim(k)[1]==m)
        stop("k must have equal rows to y")
      if(is.vector(k))
        {
          if(!length(k)==m) stop("vector k length must be equal to x rows")
          k<-matrix(k,n,1)
        }
      if(!dim(x)[2]==dim(y)[2])
        stop("matrixes must have the same number of columns")
      dotb <- rowSums(y*y)/2
      res <- x%*%t(y)
      for( i in 1:m)#2*sigma or sigma
        res[,i]<- k[i,]*(exp(2*sigma*(res[,i] - dota - rep(dotb[i],n)))*z)
      return(res)
    }
}
setMethod("kernelPol",signature(kernel="rbfkernel", x="matrix"),kernelPol.rbfkernel)

kernelPol.laplacekernel <- function(kernel, x, y=NULL, z, k=NULL)
{
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix or NULL")
  if(!is.matrix(z)&&!is.vector(z)) stop("z must be a matrix or a vector")
  if(!is.matrix(k)&&!is.vector(k)&&!is.null(k)) stop("k must be a matrix or a vector")
  sigma <- kpar(kernel)$sigma
  n <- dim(x)[1]
  dota <- rowSums(x*x)/2
  if(is.vector(z))
    {
      if(!length(z)==n) stop("vector z length must be equal to x rows")
      z<-matrix(z,n,1)
    }
  if(!dim(z)[1]==n)
    stop("z must have the length equal to x colums")
  if (is.null(y))
    {
      if(is.matrix(z)&&!dim(z)[1]==n)
       stop("z must have size equal to x colums")
      res <- crossprod(t(x))
      for (i in 1:n)
        res[i,] <- z[i,]*(exp(-sigma*sqrt(-round(2*(res[i,] - dota - rep(dota[i],n)),9)))*z)
      return(res)
    }
  if (is.matrix(y))
    {
      if(is.null(k)) stop("k not specified!")
      m <- dim(y)[1]
      if(!dim(k)[1]==m)
        stop("k must have equal rows to y")
      if(is.vector(k))
        {
          if(!length(k)==m) stop("vector k length must be equal to x rows")
          k<-matrix(k,n,1)
        }
      if(!dim(x)[2]==dim(y)[2])
        stop("matrixes must have the same number of columns")
      dotb <- rowSums(y*y)/2
      res <- x%*%t(y)
      for( i in 1:m)#2*sigma or sigma
        res[,i]<- k[i,]*(exp(-sigma*sqrt(-round(2*(res[,i] - dota - rep(dotb[i],n)),9)))*z)
      return(res)
    }
}
setMethod("kernelPol",signature(kernel="laplacekernel", x="matrix"),kernelPol.laplacekernel)


kernelPol.besselkernel <- function(kernel, x, y=NULL, z, k=NULL)
{
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix or NULL")
  if(!is.matrix(z)&&!is.vector(z)) stop("z must be a matrix or a vector")
  if(!is.matrix(k)&&!is.vector(k)&&!is.null(k)) stop("k must be a matrix or a vector")
  sigma <- kpar(kernel)$sigma
  nu <- kpar(kernel)$order
  ni <- kpar(kernel)$degree
  n <- dim(x)[1]
  lim <- 1/(gamma(nu + 1)*2^nu)
  dota <- rowSums(x*x)/2
  if(is.vector(z))
    {
      if(!length(z)==n) stop("vector z length must be equal to x rows")
      z<-matrix(z,n,1)
    }
  if(!dim(z)[1]==n)
    stop("z must have the length equal to x colums")
  if (is.null(y))
    {
      if(is.matrix(z)&&!dim(z)[1]==n)
       stop("z must have size equal to x colums")
      res <- crossprod(t(x))
      for (i in 1:n)
        {
        xx <- sigma*sqrt(-round(2*(res[i,] - dota - rep(dota[i],n)),9))
        res[i,] <- besselJ(xx,nu)*(xx^(-nu))
        res[i,which(xx < 10e-5)] <- lim
        res[i,] <- z[i,]*(((res[i,]/lim)^ni)*z)
      }
      return(res)
    }
  if (is.matrix(y))
    {
      if(is.null(k)) stop("k not specified!")
      m <- dim(y)[1]
      if(!dim(k)[1]==m)
        stop("k must have equal rows to y")
      if(is.vector(k))
        {
          if(!length(k)==m) stop("vector k length must be equal to x rows")
          k<-matrix(k,n,1)
        }
      if(!dim(x)[2]==dim(y)[2])
        stop("matrixes must have the same number of columns")
      dotb <- rowSums(y*y)/2
      res <- x%*%t(y)
      for( i in 1:m){#2*sigma or sigma
        xx <- sigma*sqrt(-round(2*(res[,i] - dota - rep(dotb[i],n)),9))
        res[,i] <- besselJ(xx,nu)*(xx^(-nu))
        res[which(xx<10e-5),i] <- lim
        res[,i]<- k[i,]*(((res[,i]/lim)^ni)*z)
      }
      return(res)
    }
}
setMethod("kernelPol",signature(kernel="besselkernel", x="matrix"),kernelPol.besselkernel)


kernelPol.anovakernel <- function(kernel, x, y=NULL, z, k=NULL)
{
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix or NULL")
  if(!is.matrix(z)&&!is.vector(z)) stop("z must be a matrix or a vector")
  if(!is.matrix(k)&&!is.vector(k)&&!is.null(k)) stop("k must be a matrix or a vector")
  sigma <- kpar(kernel)$sigma
  degree <- kpar(kernel)$degree
  n <- dim(x)[1]
   if(is.vector(z))
    {
      if(!length(z)==n) stop("vector z length must be equal to x rows")
      z<-matrix(z,n,1)
    }
  if(!dim(z)[1]==n)
    stop("z must have the length equal to x colums")
  if (is.null(y))
    {
      if(is.matrix(z)&&!dim(z)[1]==n)
       stop("z must have size equal to x colums")
      a <- matrix(0, dim(x)[2], n)
      res <- matrix(0,n,n)
      for (i in 1:n)
        {
          a[rep(TRUE,dim(x)[2]), rep(TRUE,n)] <- x[i,]
          res[i,]<- z[i,]*((colSums(exp( - sigma*(a - t(x))^2))^degree)*z)
        }
      return(res)
    }
  if (is.matrix(y))
    {
      if(is.null(k)) stop("k not specified!")
      m <- dim(y)[1]
      if(!dim(k)[1]==m)
        stop("k must have equal rows to y")
      if(is.vector(k))
        {
          if(!length(k)==m) stop("vector k length must be equal to x rows")
          k<-matrix(k,n,1)
        }
      if(!dim(x)[2]==dim(y)[2])
        stop("matrixes must have the same number of columns")

      b <- matrix(0, dim(x)[2],m)
      res <- matrix(0, dim(x)[1],m)
      for( i in 1:n)
        {
          b[rep(TRUE,dim(x)[2]), rep(TRUE,m)] <- x[i,]
          res[i,] <- z[i,]*((colSums(exp( - sigma*(b - t(y))^2))^degree)*k)
        }
      return(res)
    }
}
setMethod("kernelPol",signature(kernel="anovakernel", x="matrix"),kernelPol.anovakernel)


kernelPol.polykernel <- function(kernel, x, y=NULL, z, k=NULL)
{
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix or NULL")
  if(!is.matrix(z)&&!is.vector(z)) stop("z must be a matrix or a vector")
  if(!is.matrix(k)&&!is.vector(k)&&!is.null(k)) stop("k must be a matrix or a vector")
  degree <- kpar(kernel)$degree
  scale <- kpar(kernl)$scale
  offset <- kpar(kernel)$offset
  n <- dim(x)[1]
  if(is.vector(z))
    {
      if(!length(z)==n) stop("vector z length must be equal to x rows")
      z<-matrix(z,n,1)
    }
  if(!dim(z)[1]==n)
    stop("z must have the length equal to x colums")
  if (is.null(y))
    {
      if(is.matrix(z)&&!dim(z)[1]==n)
       stop("z must have size equal to x colums")
      for (i in 1:n)
        res <- z*(((scale*crossprod(t(x))+offset)^degree)*z)
      return(res)
    }
  if (is.matrix(y))
    {
      if(is.null(k)) stop("k not specified!")
      m <- dim(y)[1]
      if(!dim(k)[1]==m)
        stop("k must have equal rows to y")
      if(is.vector(k))
        {
          if(!length(k)==m) stop("vector k length must be equal to x rows")
          k<-matrix(k,n,1)
        }
      if(!dim(x)[2]==dim(y)[2])
        stop("matrixes must have the same number of columns")
      for( i in 1:m)#2*sigma or sigma
        res<- k*(((scale*x%*%t(y) + offset)^degree)*z)
      return(res)
    }
}
setMethod("kernelPol",signature(kernel="polykernel", x="matrix"),kernelPol.polykernel)


kernelPol.tanhkernel <- function(kernel, x, y=NULL, z, k=NULL)
{
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix or NULL")
  if(!is.matrix(z)&&!is.vector(z)) stop("z must be a matrix or a vector")
  if(!is.matrix(k)&&!is.vector(k)&&!is.null(k)) stop("k must be a matrix or a vector")
  scale <- kpar(kernel)$scale
  offset <- kpar(kernel)$offset
  n <- dim(x)[1]
  if(is.vector(z))
    {
      if(!length(z)==n) stop("vector z length must be equal to x rows")
      z<-matrix(z,n,1)
    }
  if(!dim(z)[1]==n)
    stop("z must have the length equal to x colums")
  if (is.null(y))
    {
      if(is.matrix(z)&&!dim(z)[1]==n)
       stop("z must have size equal to x colums")
      for (i in 1:n)
        res <- z*(tanh(scale*crossprod(t(x))+offset)*z)
      return(res)
    }
  if (is.matrix(y))
    {
      if(is.null(k)) stop("k not specified!")
      m <- dim(y)[1]
      if(!dim(k)[1]==m)
        stop("k must have equal rows to y")
      if(is.vector(k))
        {
          if(!length(k)==m) stop("vector k length must be equal to x rows")
          k<-matrix(k,n,1)
        }
      if(!dim(x)[2]==dim(y)[2])
        stop("matrixes must have the same number of columns")
      for( i in 1:m)#2*sigma or sigma
        res<- k*(tanh(scale*x%*%t(y) + offset)*z)
      return(res)
    }
}
setMethod("kernelPol",signature(kernel="tanhkernel", x="matrix"),kernelPol.tanhkernel)


kernelPol.vanillakernel <- function(kernel, x, y=NULL, z, k=NULL)
{
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix or NULL")
  if(!is.matrix(z)&&!is.vector(z)) stop("z must be a matrix or a vector")
  if(!is.matrix(k)&&!is.vector(k)&&!is.null(k)) stop("k must be a matrix or a vector")
  sigma <- kpar(kernel)$sigma
  n <- dim(x)[1]
  if(is.vector(z))
    {
      if(!length(z)==n) stop("vector z length must be equal to x rows")
      z<-matrix(z,n,1)
    }
  if(!dim(z)[1]==n)
    stop("z must have the length equal to x colums")
  if (is.null(y))
    {
      if(is.matrix(z)&&!dim(z)[1]==n)
       stop("z must have size equal to x colums")
      for (i in 1:n)
        res <- z*(crossprod(t(x))*z)
      return(res)
    }
  if (is.matrix(y))
    {
      if(is.null(k)) stop("k not specified!")
      m <- dim(y)[1]
      if(!dim(k)[1]==m)
        stop("k must have equal rows to y")
      if(is.vector(k))
        {
          if(!length(k)==m) stop("vector k length must be equal to x rows")
          k<-matrix(k,n,1)
        }
      if(!dim(x)[2]==dim(y)[2])
        stop("matrixes must have the same number of columns")
      for( i in 1:m)
        res<- k*(x%*%t(y)*z)
      return(res)
    }
}
setMethod("kernelPol",signature(kernel="vanillakernel", x="matrix"),kernelPol.vanillakernel)

setMethod("show","kernel",
          function(object)
          {
            switch(class(object),
                   "rbfkernel" = cat(paste("Gaussian Radial Basis kernel function.", "\n","Hyperparameter :" ,"sigma = ", kpar(object)$sigma,"\n")),
		   "laplacekernel" = cat(paste("Laplace kernel function.", "\n","Hyperparameter :" ,"sigma = ", kpar(object)$sigma,"\n")),
                   "besselkernel" = cat(paste("Bessel kernel function.", "\n","Hyperparameter :" ,"sigma = ", kpar(object)$sigma,"order = ",kpar(object)$order, "degree = ", kpar(object)$degree,"\n")),
                    "anovakernel" = cat(paste("Anova RBF kernel function.", "\n","Hyperparameter :" ,"sigma = ", kpar(object)$sigma, "degree = ", kpar(object)$degree,"\n")),
                   "tanhkernel" = cat(paste("Hyperbolic Tangent kernel function.", "\n","Hyperparameters :","scale = ", kpar(object)$scale," offset = ", kpar(object)$offset,"\n")),
                   "polykernel" = cat(paste("Polynomial kernel function.", "\n","Hyperparameters :","degree = ",kpar(object)$degree," scale = ", kpar(object)$scale," offset = ", kpar(object)$offset,"\n")),
                   "vanillakernel" = cat(paste("Linear (vanilla) kernel function.", "\n"))
                     )
                 })
