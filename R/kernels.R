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
setClass("tanhkernel",prototype=structure(.Data=function(){},kpar=list()),contains=c("kernel"))

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
  return(new("tanhkernel",.Data=rval,kpar=list(degree=degree,scale=scale,offset=offset)))
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
      res <- t(crossprod((x%*%z),t(a)))
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
      res <- t(crossprod((y%*%z),t(x)))
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
                   "tanhkernel" = cat(paste("Hyperbolic Tangent kernel function.", "\n","Hyperparameters :","scale = ", kpar(object)$scale," offset = ", kpar(object)$offset,"\n")),
                   "polykernel" = cat(paste("Polynomial kernel function.", "\n","Hyperparameters :","degree = ",kpar(object)$degree," scale = ", kpar(object)$scale," offset = ", kpar(object)$offset,"\n")),
                   "vanillakernel" = cat(paste("Linear (vanilla) kernel function.", "\n"))
                     )
                 })
