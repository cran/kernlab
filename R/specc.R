

#specc function

setGeneric("specc",function(x, ...) standardGeneric("specc"))
setMethod("specc", signature(x = "formula"),
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
    x <- model.matrix(mt, mf)
    res <- specc(x, ...)
   
    cl[[1]] <- as.name("specc")
    kcall(res) <- cl
    if(!is.null(na.act)) 
        n.action(res) <- na.action
   
    return(res)
  })

setMethod("specc",signature(x="matrix"),
          function(x, centers, kernel = "rbfdot", kpar = list(sigma = 0.1), iterations = 200, na.action = na.omit, ...)
{
  x <- na.action(x)
  x <- as.matrix(x)
  m <- nrow(x)
  if (missing(centers))
    stop("centers must be a number or a matrix")
  if (length(centers) == 1) {
    nc <-  centers
    if (m < centers)
      stop("more cluster centers than data points.")
  }
  else
    nc <- dim(centers)[2]
  
  ret <- new("specc")
  if(!is(kernel,"kernel"))
    {
      if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
      kernel <- do.call(kernel, kpar)
    }
  if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

  km <- kernelMatrix(kernel,x)
  if(is(kernel)[1] == "rbfkernel")
    diag(km) <- 0
  
  d <- 1/sqrt(rowSums(km))
  l <- d * km %*% diag(d)
  xi <- eigen(l)$vectors[,1:nc]
  yi <- xi/sqrt(rowSums(xi^2))
  res <- kmeans(yi, centers, iterations)

  cluster(ret) <- res$cluster 
  centers(ret) <- res$centers
  size(ret) <- res$size
  kernelf(ret) <- kernel
## use res$withinss for model selection !!
  return(ret)
})


  
