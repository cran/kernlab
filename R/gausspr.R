setGeneric("gausspr", function(x, ...) standardGeneric("gausspr"))
setMethod("gausspr",signature(x="formula"),
function (x, data=NULL, ..., subset, na.action = na.omit){
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- NULL
  m$formula <- m$x
  m$x <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  attr(Terms, "intercept") <- 0
  x <- model.matrix(Terms, m)
  y <- model.extract(m, response)
  ret <- gausspr(x, y, ...)
  kcall(ret) <- call
  kterms(ret) <- Terms
  if (!is.null(attr(m, "na.action")))
    n.action(ret) <- attr(m, "na.action")
  return (ret)
})

setMethod("gausspr",signature(x="vector"),
function(x,...)
  {
    x <- t(t(x))
    ret <- gausspr(x, ...)
    ret
  })
    
setMethod("gausspr",signature(x="matrix"),
function (x,
          y         = NULL,
          type      = NULL,
          kernel    = "rbfdot",
          kpar      = list(sigma = 0.1),
          var       = 1,
          tol       = 0.001,  
          cross     = 0,
          fit       = TRUE,
          ...
          ,subset 
         ,na.action = na.omit)
{

## subsetting and na-handling for matrices
  ret <- new("gausspr")
  if (!missing(subset)) x <- x[subset,]
  if (is.null(y))
    x <- na.action(x)
  else {
    df <- na.action(data.frame(y, x))
    y <- df[,1]
    x <- as.matrix(df[,-1])
  }
  ncols <- ncol(x)
  m <- nrows <- nrow(x)
  
 if (is.null (type)) type(ret) <-
   if (is.factor(y)) "classification"
    else "regression"
  else type(ret) <- type
  
  if (var < 10^-3)
    stop("Noise variance parameter var has to be greater than 10^-3")

 # in case of classification: transform factors into integers
  if (is.factor(y)) {
    lev(ret) <- levels (y)
    y <- as.integer (y)
  }
  else {
    if (type(ret) == "classification" && any(as.integer (y) != y))
        stop ("dependent variable has to be of factor or integer type for classification mode.")
    if(type(ret) == "classification")
      lev(ret) <- unique (y)
    }
 # initialize    
  nclass(ret) <- length (lev(ret))
  
  if(!is.null(type))
    type(ret) <- match.arg(type,c("classification", "regression"))
  
  if(!is(kernel,"kernel"))
    {
      if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
      kernel <- do.call(kernel, kpar)
    }
  if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

  p <- 0
  
  if (type(ret) == "classification")
    {
      indexes <- lapply(1:nclass(ret), function(kk) which(y == kk))
      for (i in 1:(nclass(ret)-1)) {
        jj <- i+1
        for(j in jj:nclass(ret)) {
          p <- p+1
          ##prepare data
          li <- length(indexes[[i]])
          lj <- length(indexes[[j]])
          xd <- matrix(0,(li+lj),dim(x)[2])
          xdi <- 1:(li+lj) <= li
          xd[xdi,rep(TRUE,dim(x)[2])] <- x[indexes[[i]],]
          xd[xdi == FALSE,rep(TRUE,dim(x)[2])] <- x[indexes[[j]],]
          if(y[indexes[[i]][1]] < y[indexes[[j]]][1])
            yd <- c(rep(1,li),rep(-1,lj))
          else
            yd <- c(rep(-1,li),rep(1,lj))
          K <- kernelMatrix(kernel,xd)
          gradnorm <- 1 
          alphag <- rep(0,li+lj)
          while (gradnorm > tol)
            {
              f <- crossprod(K,alphag)
              grad <- -yd/(1 + exp(yd*f))
              hess <- exp(yd*f)
              hess <- hess / ((1 + hess)^2)
              alphag <- alphag - solve(crossprod(K,diag(as.vector(hess))) + diag(rep(var,li+lj)))%*%(grad + alphag)
              gradnorm <- sqrt(sum((grad + alphag)^2))
            }
          alpha(ret)[[p]] <- alphag
          alphaindex(ret)[[p]] <- c(indexes[[i]],indexes[[j]])
        }
      }
    }
  
  if (type(ret) == "regression")
    {
      K <- kernelMatrix(kernel,x)

      #alpha(ret) <- solve(K + diag (rep(var, length=m)), y)
      alpha(ret) <- solve(K + diag(rep(var, length = m))) %*% y
    }

  kernelf(ret) <- kernel
  xmatrix(ret) <- x

  fit(ret)  <- if (fit)
    predict(ret, x) else NA

  if (fit){
    if(type(ret)=="classification")
      error(ret) <- 1 - .classAgreement(table(y,as.integer(fit(ret))))
    if(type(ret)=="regression")
      error(ret) <- drop(crossprod(fit(ret) - y)/m)
  }

  if(cross!=0)
    {
      cerror <- 0
      suppressWarnings(vgr<-split(sample(1:m,m),1:cross))
      for(i in 1:cross)
        {
          cind <- unsplit(vgr[-i],1:(m-length(vgr[[i]])))
          if(type(ret)=="classification")
            {
              cret <- gausspr(x[cind,],factor (lev(ret)[y[cind]], levels = lev(ret)),type=type(ret),kernel=kernel,C=C,var = var, cross = 0, fit = FALSE)
               cres <- predict(cret, x[vgr[[i]],])
            cerror <- (1 - .classAgreement(table(y[vgr[[i]]],as.integer(cres))))/cross + cerror
            }
          if(type(ret)=="regression")
            {
              cret <- gausspr(x[cind,],y[cind],type=type(ret),kernel=kernel,var = var,tol=tol, cross = 0, fit = FALSE)
              cres <- predict(cret, x[vgr[[i]],])
              cerror <- drop(crossprod(cres - y[vgr[[i]]])/m)/cross + cerror
            }
        }
      cross(ret) <- cerror
    }


  
  return(ret)
})


setMethod("predict", signature(object = "gausspr"),
function (object, newdata, type = "response", coupler = "minpair",...)
{
  if (missing(newdata))
    return(fit(object))
  ncols <- ncol(xmatrix(object))
  nrows <- nrow(xmatrix(object))
  oldco <- ncols

  if (!is.null(kterms(object)))
  {  
  newdata <- model.matrix(delete.response(kterms(object)), as.data.frame(newdata), na.action = na.action)
   }
  else
    newdata  <- if (is.vector (newdata)) t(t(newdata)) else as.matrix(newdata)

  newcols  <- 0
  newnrows <- nrow(newdata)
  newncols <- ncol(newdata)
  newco    <- newncols
    
  if (oldco != newco) stop ("test vector does not match model !")

  type <- match.arg(type,c("response","probabilities","votes"))
  p <- 0
  if(type == "response")
    {
  if(type(object)=="classification")
    {
      predres <- 1:newnrows
      votematrix <- matrix(0,nclass(object),nrows)
      for(i in 1:(nclass(object)-1))
        {
        jj <- i+1
        for(j in jj:nclass(object))
          {
            p <- p+1
            ret <- kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[p]],],alpha(object)[[p]])
            votematrix[i,ret>0] <- votematrix[i,ret>0] + 1
            votematrix[j,ret<0] <- votematrix[j,ret<0] + 1
          }
      }
      predres <- sapply(predres, function(x) which.max(votematrix[,x]))
    }

}

  if(type == "probabilities")
    {
      if(type(object)=="classification")
        {
          binprob <- matrix(0, newnrows, nclass(object)*(nclass(object) - 1)/2)
          for(i in 1:(nclass(object)-1))
            {
              jj <- i+1
              for(j in jj:nclass(object))
                {
                  p <- p+1
                  binprob[,p] <-  1/(1+exp(-kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[p]],],alpha(object)[[p]])))

                }
            }
          ## multiprob <- sapply(1:newnrows, function(x) couple(binprob[x ,],coupler = coupler))
          multiprob <- couple(binprob, coupler = coupler)
        }
    }
    
  
  if(type(object) == "regression")
    {
      predres <- kernelMult(kernelf(object),newdata,xmatrix(object),as.matrix(alpha(object)))
    }


 if (is.character(lev(object)))
    {
      ##classification & probabilities : return probabilitie matrix
      if(type == "probabilities")
        {
          colnames(multiprob) <- lev(object)
          return(multiprob)
        }
      ##classification & type response: return factors
      if(type == "response")
        return(factor (lev(object)[predres], levels = lev(object)))
      ##classification & votes : return votematrix
      if(type == "votes")
        return(votematrix)
    }
  else if (type(object) == "one-classification")
    ##one-class-classification: return TRUE/FALSE (probabilities ?)
    return(ret == 1)
  else
    ##else: return raw values
    return(predres)

})


setMethod("show","gausspr",
function(object){
  cat("Gaussian Processes object of class \"gausspr\"","\n")
  cat(paste("Problem type:", type(object),"\n"))
  cat("\n")
  show(kernelf(object))
  cat(paste("\nNumber of training instances learned :", length(alpha(object)),"\n"))
  if(!is.null(fit(object)))
    cat(paste("Train error :", round(error(object),9),"\n"))
  ##train error & loss
  if(!is.null(cross(object)))
    cat("Cross validation error :",round(cross(object),9),"\n")
})
