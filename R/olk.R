## OLK algorithms for classification, novelty detection and regression. 
##
## derived from onlearn 2.6.15 Stephen Tridgell
## updated

setGeneric("olk",function(obj, x, y = NULL, olkC = 0.5, r = 1e-4, epsilon = 0.005) standardGeneric("olk"))
setMethod("olk", signature(obj = "olk"),
          function(obj , x, y = NULL, olkC = 0.5, r = 1e-4, epsilon = 0.005)
{
    if(olkstart(obj) == 1 && olkstop(obj) < bufferolk(obj))
        buffernotfull  <- TRUE
    else
        buffernotfull <- FALSE

    if(is.vector(x))
        x <- matrix(x,,length(x))
    d <- dim(x)[2]
        
    for (i in 1:dim(x)[1])
    {
        # for all available examples
        xt <- x[i,,drop=FALSE]
        yt <- y[i]
        if(type(obj)=="novelty")
        {
            ft <- 1 + r - fitolk(obj)
            if(ft > 0)
            {
                if(ft > olkC)
                    ft <- olkC
                alpha(obj) <- (1/(1+r)) * alpha(obj)
                if(buffernotfull)
                    olkstop(obj) <- olkstop(obj) + 1
                else{
                    olkstop(obj) <- olkstop(obj)%%bufferolk(obj) + 1
                    olkstart(obj) <- olkstart(obj)%%bufferolk(obj) +1
                }
                alpha(obj)[olkstop(obj)] <- ft/(1+r)
                xmatrix(obj)[olkstop(obj),] <- xt
            }

            if(olkstart(obj) == 1 && olkstop(obj) < bufferolk(obj))
                fitolk(obj) <- drop(kernelMult(kernelf(obj), xt, matrix(xmatrix(obj)[1:olkstop(obj),],ncol=d),
                                               matrix(alpha(obj)[1:olkstop(obj)],ncol=1))) 
            else
                fitolk(obj) <- drop(kernelMult(kernelf(obj), xt, xmatrix(obj), matrix(alpha(obj),ncol=1)))
        }
        if(type(obj)=="classification")
        { 
            if(is.null(pattern(obj)) && is.factor(y))
                pattern(obj) <- yt
            if(!is.null(pattern(obj)))
                if(pattern(obj) == yt)
                    yt <- 1
                else yt <-  -1
            
            ft <- 1 + r - yt*fitolk(obj)
            if (ft > 0)
            {
                if (ft > olkC)
                    ft <- olkC
                alpha(obj) <- (1/(1+r)) * alpha(obj)
                if(buffernotfull)
                    olkstop(obj) <- olkstop(obj) + 1
                else{
                    olkstop(obj) <- olkstop(obj)%%bufferolk(obj) + 1
                    olkstart(obj) <- olkstart(obj)%%bufferolk(obj) +1
                }
                alpha(obj)[olkstop(obj)] <- yt*ft/(1+r)
                xmatrix(obj)[olkstop(obj),] <- xt
            }

            if(olkstart(obj) == 1 && olkstop(obj) < bufferolk(obj))
                fitolk(obj) <- drop(kernelMult(kernelf(obj), xt, xmatrix(obj)[1:olkstop(obj),,drop=FALSE],
                                               matrix(alpha(obj)[1:olkstop(obj)],ncol=1)))
            else
                fitolk(obj) <-drop(kernelMult(kernelf(obj), xt, xmatrix(obj), matrix(alpha(obj),ncol=1)))
            
        }

        if(type(obj)=="regression")
        {
            ft <- fitolk(obj)
            alpha1 <- ft - (1+r)*(yt + epsilon)
            alpha2 <- (1+r)*(yt - epsilon) - ft
            if(alpha1 > 0 || alpha2 > 0)
            {
                if (alpha1 < 0)
                    alpha1 <- 0
                if (alpha2 < 0)
                    alpha2 <- 0
                if (alpha1 > olkC)
                    alpha1 <- olkC
                if (alpha2 > olkC)
                    alpha2 <- olkC
                alpha(obj) <- (1/(1+r)) * alpha(obj)
                if(buffernotfull)
                    olkstop(obj) <- olkstop(obj) + 1
                else{
                    olkstop(obj) <- olkstop(obj)%%bufferolk(obj) + 1
                    olkstart(obj) <- olkstart(obj)%% bufferolk(obj) +1
                }
                alpha(obj)[olkstop(obj)] <- (alpha1 - alpha2)/(1+r)
                xmatrix(obj)[olkstop(obj),] <- xt
            }

            if(olkstart(obj) == 1 && olkstop(obj) < bufferolk(obj))
                fitolk(obj) <- drop(kernelMult(kernelf(obj), xt, matrix(xmatrix(obj)[1:olkstop(obj),],ncol=d),
                                               matrix(alpha(obj)[1:olkstop(obj)],ncol=1)))
            else
                fitolk(obj) <- drop(kernelMult(kernelf(obj), xt, xmatrix(obj), matrix(alpha(obj),ncol=1))) 
        }
    }
    return(obj)
})


setGeneric("inlearn",function(d, kernel = "rbfdot", kpar = list(sigma=0.1), type = "novelty", buffersize = 1000) standardGeneric("inlearn"))
setMethod("inlearn", signature(d = "numeric"),
          function(d ,kernel = "rbfdot", kpar = list(sigma=0.1), type = "novelty", buffersize = 1000)
{
    obj <- new("olk")

    if(!is(kernel,"kernel"))
    {
        if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
        kernel <- do.call(kernel, kpar)
    }
    if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

    type(obj) <- match.arg(type,c("novelty","classification","regression"))
    xmatrix(obj) <- matrix(0,buffersize,d)
    kernelf(obj) <- kernel
    olkstart(obj) <- 1
    olkstop(obj) <- 1
    fitolk(obj) <- 0 
    alpha(obj) <- rep(0, buffersize)
    bufferolk(obj) <- buffersize
    return(obj)
})


setMethod("show","olk",
          function(object){
    cat("Online learning with kernels of class \"olk\"","\n")
    cat("\n")
    cat(paste("Learning problem :", type(object), "\n"))
    cat
    cat(paste("Data dimensions :", dim(xmatrix(object))[2], "\n"))
    cat(paste("Buffersize :", bufferolk(object), "\n"))
    cat("\n")
    show(kernelf(object))
})


setMethod("predict",signature(object="olk"),
          function(object, x)
{
    if(is.vector(x))
        x<- matrix(x,1)

    d <- dim(xmatrix(object))[2]
    
    if(type(object)=="novelty")
    {
        if(olkstart(object) == 1 && olkstop(object) < bufferolk(object))
            res <- drop(kernelMult(kernelf(object), x, matrix(xmatrix(object)[1:olkstop(object),],ncol= d),
                                   matrix(alpha(object)[1:olkstop(object)],ncol=1)) - 1) 
        else
            res <- drop(kernelMult(kernelf(object), x, matrix(xmatrix(object),ncol=d), matrix(alpha(object)),ncol=1) - 1)
    }

    if(type(object)=="classification")
    {
        if(olkstart(object) == 1 && olkstop(object) < bufferolk(object))
            res <- drop(kernelMult(kernelf(object), x, matrix(xmatrix(object)[1:olkstop(object),],ncol=d),
                                   matrix(alpha(object)[1:olkstop(object)],ncol=1)))
        else
            res <- drop(kernelMult(kernelf(object), x, matrix(xmatrix(object),ncol=d), matrix(alpha(object)),ncol=1))
        
    }

    if(type(object)=="regression")
    {
        if(olkstart(object) == 1 && olkstop(object) < bufferolk(object))
            res <- drop(kernelMult(kernelf(object), x, matrix(xmatrix(object)[1:olkstop(object),],ncol=d),
                                   matrix(alpha(object)[1:olkstop(object)],ncol=1))) 
        else
            res <- drop(kernelMult(kernelf(object), x, matrix(xmatrix(object),ncol=d), matrix(alpha(object)),ncol=1)) 
    }

    return(res)
    
})


