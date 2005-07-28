

#kpca function

setGeneric("kpca",function(x, ...) standardGeneric("kpca"))
setMethod("kpca", signature(x = "formula"),
function(x, data = NULL, na.action = na.omit, ...)
{
    mt <- terms(x, data = data)
    if(attr(mt, "response") > 0) stop("response not allowed in formula")
    attr(mt, "intercept") <- 0
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    mf$formula <- mf$x
    mf$... <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    na.act <- attr(mf, "na.action")
    Terms <- attr(mf, "terms")
    x <- model.matrix(mt, mf)
    res <- kpca(x, ...)
    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1]] <- as.name("kpca")
    kcall(res) <- cl
    attr(Terms,"intercept") <- 0
    kterms(res) <- Terms
    if(!is.null(na.act)) 
        n.action(res) <- na.act
  
    return(res)
  })



setMethod("kpca",signature(x="matrix"),
          function(x, kernel = "rbfdot", kpar = list(sigma = 0.1), features = 0, th = 1e-4, na.action = na.omit, ...)
{
  x <- na.action(x)
  x <- as.matrix(x)
  m <- nrow(x)
  ret <- new("kpca")
  if(!is(kernel,"kernel"))
    {
      if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
      kernel <- do.call(kernel, kpar)
    }
  if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

  km <- kernelMatrix(kernel,x)
  #center kernel matrix
  i <- matrix((1/m), m, m)
  kc <- km - km %*% i - i %*% km + i %*% km %*% i
  #compute eigenvectors
  res <- eigen(kc/m,symmetric=TRUE)
  
  if(features == 0)
    features <- sum(res$values > th)
  else 
    if(res$values[features] < th)
      warning(paste("eigenvalues of the kernel matrix are below threshold!"))
 
  pcv(ret) <- t(t(res$vectors[,1:features])/sqrt(res$values[1:features]))
  eig(ret) <- res$values[1:features]
  names(eig(ret)) <- paste("Comp.", 1:features, sep = "")
  rotated(ret) <- kc %*% pcv(ret)
  kernelf(ret) <- kernel
  xmatrix(ret) <- x
  return(ret)
})

#project a new matrix into the feature space 
setMethod("predict",signature(object="kpca"),
function(object , x)
  {
    if (!is.null(kterms(object)))
      {
        if(!is.matrix(x))
          x <- model.matrix(delete.response(kterms(object)), as.data.frame(x), na.action = n.action(object))
      }
    else
      x  <- if (is.vector(x)) t(t(x)) else as.matrix(x)

    if (is.vector(x)||is.data.frame(x))
      x<-as.matrix(x)
    if (!is.matrix(x)) stop("x must be a matrix a vector or a data frame")
    n <- nrow(x)
    m <- nrow(xmatrix(object))
    knc <- kernelMatrix(kernelf(object),x,xmatrix(object))
    ka <- kernelMatrix(kernelf(object),xmatrix(object))
    #center
    yi <- matrix((1/m), m, m)
    xi <- matrix((1/m), n, m)
    ret <- knc - knc %*% yi - xi %*% ka + xi %*% ka %*% yi
    return(ret %*% pcv(object))
  })



  
