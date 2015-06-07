## kernel based on-line learning algorithms for classification, novelty detection and regression.
##
## created 15.09.04 alexandros
## updated 04.06.15 Stephen Tridgell

setGeneric("onlearn",function(obj, x, y = NULL, nu = 0.2, lambda = 0.99, eta = 0.001) standardGeneric("onlearn"))
setMethod("onlearn", signature(obj = "onlearn"),
          function(obj , x, y = NULL,  nu = 0.2, lambda = 0.99, eta = 0.001)
          {
            if(onstart(obj) == 1 && onstop(obj) < buffer(obj))
              buffernotfull  <- TRUE
            else
              buffernotfull <- FALSE

            if(is.vector(x))
              x <- matrix(x,,length(x))  
            d <- dim(x)[2]
            if (is.factor(y)) {
                lev(obj) <- levels(y)
                y <- as.integer (y)
                if (type(obj) == "classification" && length(lev(obj)) == 2)
                    y <- (2*(y - 1) - 1)
            }
            else {
                if ((type(obj) != "classification") && any(as.integer (y) != y))
                    stop ("dependent variable has to be of factor or integer type for classification mode.")
                if (type(obj) != "regression")
                    lev(obj) <- sort(unique (y))
            }
            if (type(obj) == "classification") {
                if (length(unique(y)) > 2)  stop ("Only two class classification is supported")
                if (length(unique(y)) == 1 && !(max(y) == 1 || min(y) == -1))
                    warning("y must be either 1 or -1")
            }
            for (i in 1:dim(x)[1])
              {
                xt <- x[i,,drop=FALSE]
                yt <- y[i]
            if(type(obj)=="novelty")
              {
                phi <- fit(obj)

                alpha(obj) <- (1 - eta) * alpha(obj)

                if(phi < 0)
                  {
                    if(buffernotfull)
                      onstop(obj) <- onstop(obj) + 1
                    else{
                      onstop(obj) <- onstop(obj)%%buffer(obj) + 1
                      onstart(obj) <- onstart(obj)%%buffer(obj) +1
                    }
                    alpha(obj)[onstop(obj)] <- eta
                    xmatrix(obj)[onstop(obj),] <- xt
                    rho(obj) <- rho(obj) - eta*(1 - nu)
                  }
                else
                  rho(obj) <- rho(obj) - eta*nu  # sign error in Online Learning with Kernels

                if(onstart(obj) == 1 && onstop(obj) < buffer(obj))
                  fit(obj) <- drop(kernelMult(kernelf(obj), xt, matrix(xmatrix(obj)[1:onstop(obj),],ncol=d),
                                              matrix(alpha(obj)[1:onstop(obj)],ncol=1)) - rho(obj))
                else
                  fit(obj) <- drop(kernelMult(kernelf(obj), xt, xmatrix(obj), matrix(alpha(obj),ncol=1)) - rho(obj))
              }
            if(type(obj)=="classification")
              { 
                phi <- fit(obj)
                
                alpha(obj) <- (1 - eta) * alpha(obj)

                if(yt*phi < rho(obj))
                  {
                    if(buffernotfull)
                      onstop(obj) <- onstop(obj) + 1
                    else{
                      onstop(obj) <- onstop(obj)%%buffer(obj) + 1
                      onstart(obj) <- onstart(obj)%%buffer(obj) +1
                    }
                    alpha(obj)[onstop(obj)] <- eta*yt
                    b(obj) <- b(obj)  + eta*yt
                    xmatrix(obj)[onstop(obj),] <- xt
                    rho(obj) <- rho(obj) - eta*(1 - nu)
                  }
                else
                  rho(obj) <- rho(obj) + eta*nu  # sign error in Online Learning with Kernels
                
                if(onstart(obj) == 1 && onstop(obj) < buffer(obj))
                  fit(obj) <- drop(kernelMult(kernelf(obj), xt, xmatrix(obj)[1:onstop(obj),,drop=FALSE],
                                              matrix(alpha(obj)[1:onstop(obj)],ncol=1)) + b(obj))
                else
                  fit(obj) <-drop(kernelMult(kernelf(obj), xt, xmatrix(obj), matrix(alpha(obj),ncol=1)) + b(obj))
          
              }
            if(type(obj)=="regression")
              {
                # Using epsilon insensitive loss function
                # name rho is used instead of epsilon to have one variable across algorithms
                phi <- fit(obj)
                
                alpha(obj) <- (1 - eta*lambda) * alpha(obj)

                if(abs(yt - phi) > rho(obj))
                  {
                    if(buffernotfull)
                      onstop(obj) <- onstop(obj) + 1
                    else{
                      onstop(obj) <- onstop(obj)%%buffer(obj) + 1
                      onstart(obj) <- onstart(obj)%% buffer(obj) +1
                    }
                    alpha(obj)[onstop(obj)] <- sign(yt - phi)*eta
                    xmatrix(obj)[onstop(obj),] <- xt
                    rho(obj) <- rho(obj) + eta*(1 - nu)
                  }
                else
                    rho(obj) <- rho(obj) - eta*nu

                if(onstart(obj) == 1 && onstop(obj) < buffer(obj))
                  fit(obj) <- drop(kernelMult(kernelf(obj), xt, matrix(xmatrix(obj)[1:onstop(obj),],ncol=d),
                                              matrix(alpha(obj)[1:onstop(obj)],ncol=1)))
                else
                  fit(obj) <- drop(kernelMult(kernelf(obj), xt, xmatrix(obj), matrix(alpha(obj),ncol=1)))
              }
              }
          return(obj)
          })


setGeneric("inlearn",function(d, kernel = "rbfdot", kpar = list(sigma=0.1), type = "novelty", buffersize = 1000) standardGeneric("inlearn"))
setMethod("inlearn", signature(d = "numeric"),
          function(d ,kernel = "rbfdot", kpar = list(sigma=0.1), type = "novelty", buffersize = 1000)
          {
            obj <- new("onlearn")

            if(!is(kernel,"kernel"))
              {
                if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
                kernel <- do.call(kernel, kpar)
              }
            if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

            type(obj) <- match.arg(type,c("novelty","classification","regression"))
            xmatrix(obj) <- matrix(0,buffersize,d)
            kernelf(obj) <- kernel
            onstart(obj) <- 1
            onstop(obj) <- 1
            fit(obj) <- 0 
            b(obj) <- 0
            alpha(obj) <- rep(0, buffersize)
            rho(obj) <- 0
            buffer(obj) <- buffersize
            return(obj)
          })


setMethod("show","onlearn",
function(object){
  cat("On-line learning object of class \"onlearn\"","\n")
  cat("\n")
  cat(paste("Learning problem :", type(object), "\n"))
  cat
  cat(paste("Data dimensions :", dim(xmatrix(object))[2], "\n"))
  cat(paste("Buffersize :", buffer(object), "\n"))
  cat("\n")
 show(kernelf(object))
})


setMethod("predict",signature(object="onlearn"),
function(object, x)
  {
    if(is.vector(x))
      x<- matrix(x,1)

    d <- dim(xmatrix(object))[2]
          
    if(type(object)=="novelty")
      {
        if(onstart(object) == 1 && onstop(object) < buffer(object))
          res <- drop(kernelMult(kernelf(object), x, matrix(xmatrix(object)[1:onstop(object),],ncol= d),
                                 matrix(alpha(object)[1:onstop(object)],ncol=1)) - rho(object))
        else
          res <- drop(kernelMult(kernelf(object), x, matrix(xmatrix(object),ncol=d), matrix(alpha(object)),ncol=1) - rho(object))
      }

    if(type(object)=="classification")
      {
        if(onstart(object) == 1 && onstop(object) < buffer(object))
          res <- drop(kernelMult(kernelf(object), x, matrix(xmatrix(object)[1:onstop(object),],ncol=d),
                                 matrix(alpha(object)[1:onstop(object)],ncol=1)))
        else
          res <- drop(kernelMult(kernelf(object), x, matrix(xmatrix(object),ncol=d), matrix(alpha(object)),ncol=1))
       
      }

    if(type(object)=="regression")
      {
        if(onstart(object) == 1 && onstop(object) < buffer(object))
          res <- drop(kernelMult(kernelf(object), x, matrix(xmatrix(object)[1:onstop(object),],ncol=d),
                                 matrix(alpha(object)[1:onstop(object)],ncol=1)))
        else
          res <- drop(kernelMult(kernelf(object), x, matrix(xmatrix(object),ncol=d), matrix(alpha(object)),ncol=1))
      }

    return(res)
    
  })

  
