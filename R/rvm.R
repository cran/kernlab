## relevance vector machine
## author : alexandros

setGeneric("rvm", function(x, ...) standardGeneric("rvm"))
setMethod("rvm",signature(x="formula"),
function (x, data=NULL, ..., subset, na.action = na.omit){
  cl <- match.call()
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
   ret <- rvm(x, y, ...)
  kcall(ret) <- cl
  terms(ret) <- Terms
  if (!is.null(attr(m, "na.action")))
    n.action(ret) <- attr(m, "na.action")
  return (ret)
})

setMethod("rvm",signature(x="vector"),
function(x,...)
  {
    x <- t(t(x))
    ret <- rvm(x, ...)
    ret
  })
    
setMethod("rvm",signature(x="matrix"),
function (x,
          y         = NULL,
          type      = "regression",
          kernel    = "rbfdot",
          kpar      = list(sigma = 0.1),
          alpha     = ncol(as.matrix(x)),
          var = 0.1,        # variance 
          var.fix = FALSE,  # fixed variance?
          iterations = 100, # no. of iterations
          verbosity = 0,
          tol = .Machine$double.eps,
          minmaxdiff = 1e-3,
          cross     = 0,
          fit       = TRUE,
          ...
          ,subset 
         ,na.action = na.omit)
{

## subsetting and na-handling for matrices
  ret <- new("rvm")
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
  else
    type(ret) <- "regression"

 # in case of classification: transform factors into integers
  if (is.factor(y)) {
    lev(ret) <- levels (y)
    y <- as.integer (y)
    if (!is.null(class.weights)) {
      if (is.null(names (class.weights)))
        stop ("Weights have to be specified along with their according level names !")
      weightlabels <- match (names(class.weights),lev(ret))
      if (any(is.na(weightlabels)))
        stop ("At least one level name is missing or misspelled.")
    }
  } else {
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

  if(length(alpha) == m)
    thetavec <- 1/alpha
  else
    if (length(alpha) == 1)
      thetavec <- rep(1/alpha, m)
    else stop("length of initial alpha vector is wrong (has to be one or equal with number of train data")
  
  wvec <- rep(1, m)
  piter <- iterations*0.4
  
  if (type(ret) == "regression")
    {
      K <- kernelMatrix(kernel, x)
      Kml <- crossprod(K, y)
            
      for (i in 1:iterations) {
        nzindex <- thetavec > tol
        thetavec [!nzindex]  <- wvec [!nzindex] <- 0
        Kr <- K [ ,nzindex, drop = FALSE]
        thetatmp <- thetavec[nzindex]
        n <- sum (nzindex)

        Rinv <- backsolve(chol(crossprod(Kr)/var + diag(1/thetatmp)),diag(1,n))

        ## compute the new wvec coefficients
        wvec [nzindex] <- (Rinv %*% (crossprod(Rinv, Kml [nzindex])))/var

        diagSigma <- rowSums(Rinv^2)

        ## error
        err <- sum ((y - Kr %*% wvec [nzindex])^2)

        if(var < 2e-9)
          {
            warning("Model might be overfitted")
            break
          }
        ## log some information
        if (verbosity > 0) {
          log.det.Sigma.inv <- - 2 * sum (log (diag (Rinv)))
      
          ## compute the marginal likelihood to monitor convergence
          mlike <- -1/2 * (log.det.Sigma.inv +
                              sum (log (thetatmp)) +
                              m * log (var) + 1/var * err +
                              (wvec [nzindex]^2) %*% (1/thetatmp))

          cat ("Marg. Likelihood =", formatC (mlike), "\tnRV=", n, "\tvar=", var, "\n")
        }
        
        ## compute zeta
        zeta <- 1 - diagSigma / thetatmp
        ## compute logtheta for convergence checking
        logtheta <- - log(thetavec[nzindex])
        ## update thetavec
        if(i < piter){
          thetavec [nzindex] <- wvec [nzindex]^2 / zeta
        thetavec [thetavec <= 0] <- 0 }
        else{
          thetavec [nzindex] <- (wvec [nzindex]^2/zeta -  diagSigma)/zeta
          thetavec [thetavec <= 0] <- 0 
        }

        ## Stop if largest alpha change is too small
            
        maxdiff <- max(abs(logtheta[thetavec[which(nzindex)]!=0] + log(thetavec[thetavec!=0])))

        if(maxdiff < minmaxdiff)
          break;

        ## update variance
        if (!var.fix) {
          var <- err / (m - sum (zeta))
        }
      }

      if(verbosity == 0)
        mlike(ret) <- drop(-1/2 * (-2*sum(log(diag(Rinv))) +
                              sum (log (thetatmp)) +
                              m * log (var) + 1/var * err +
                              (wvec [nzindex]^2) %*% (1/thetatmp)))

      nvar(ret) <- var
      error(ret) <- sqrt(err/m)
      if(fit)
      fitted(ret) <- Kr %*% wvec [nzindex]
      
    }

  if(type(ret)=="classification")
    {
      stop("classification with the relevance vector machine not implemented yet")
    }
  kcall(ret) <- match.call()
  kernelf(ret) <- kernel
  alpha(ret) <- wvec[nzindex]
   tol(ret) <- tol
  xmatrix(ret) <- x
  ymatrix(ret) <- y
  RVindex(ret) <- which(nzindex)
  nRV(ret) <- length(RVindex(ret))

  
  fitted(ret)  <- if (fit)
    predict(ret, x) else NA
  
  if (fit){
    if(type(ret)=="classification")
      error(ret) <- 1 - .classAgreement(table(y,as.integer(fitted(ret))))
    if(type(ret)=="regression")
      error(ret) <- drop(crossprod(fitted(ret) - y)/m)
  }
  
  cross(ret) <- -1
  if(cross!=0)
    {
      cerror <- 0
      suppressWarnings(vgr<-split(sample(1:m,m),1:cross))
      for(i in 1:cross)
        {
          cind <- unsplit(vgr[-i],factor(rep((1:cross)[-i],unlist(lapply(vgr[-i],length)))))
          if(type(ret)=="classification")
            {
              cret <- rvm(x[cind,],factor (lev(ret)[y[cind]], levels = lev(ret)),type=type(ret),kernel=kernel,alpha = alpha,var = var, var.fix=var.fix, tol=tol, cross = 0, fit = FALSE)
               cres <- predict(cret, x[vgr[[i]],])
            cerror <- (1 - .classAgreement(table(y[vgr[[i]]],as.integer(cres))))/cross + cerror
            }
          if(type(ret)=="regression")
            {
              cret <- ksvm(x[cind,],y[cind],type=type(ret),kernel=kernel,C=C,nu=nu,epsilon=epsilon,tol=tol,alpha = alpha, var = var, var.fix=var.fix, cross = 0, fit = FALSE)
              cres <- predict(cret, x[vgr[[i]],])
              cerror <- drop(crossprod(cres - y[vgr[[i]]])/m)/cross + cerror
            }
        }
      cross(ret) <- cerror
    }

  return(ret)
})



setMethod("predict", signature(object = "rvm"),
function (object, newdata, ...)
{
  if (missing(newdata))
    return(fitted(object))
  ncols <- ncol(xmatrix(object))
  nrows <- nrow(xmatrix(object))
  oldco <- ncols

  if (!is.null(terms(object)))
  {  
  newdata <- model.matrix(delete.response(terms(object)), as.data.frame(newdata), na.action = na.action)
   }
  else
    newdata  <- if (is.vector (newdata)) t(t(newdata)) else as.matrix(newdata)

  newcols  <- 0
  newnrows <- nrow(newdata)
  newncols <- ncol(newdata)
  newco    <- newncols
    
  if (oldco != newco) stop ("test vector does not match model !")
  p<-0

  if(type(object) == "regression")
    {
      ret <- kernelMult(kernelf(object),newdata,as.matrix(xmatrix(object)[RVindex(object),]),alpha(object))
    }

  ret
})

setMethod("show","rvm",
function(object){
  cat("Relevance Vector Machine object of class \"rvm\"","\n")
   cat("Problem type: regression","\n","\n")
  show(kernelf(object))
 
  cat(paste("\nNumber of Relevance Vectors :", nRV(object),"\n"))
  cat("Variance : ",round(nvar(object),9))
  cat("\n")
    if(!is.null(fitted(object)))
  cat(paste("Training error :", round(error(object),9),"\n"))
  if(!is.null(cross(object)))
    cat("Cross validation error :",round(cross(object),9),"\n")
  ##train error & loss
})
