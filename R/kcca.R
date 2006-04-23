## Simple kernel canonical corelation analysis
## author: alexandros karatzoglou

setGeneric("kcca",function(x, y, kernel="rbfdot", kpar=list(sigma = 0.1), ...) standardGeneric("kcca"))
setMethod("kcca", signature(x = "matrix"),
          function(x,y,kernel="rbfdot",kpar=list(sigma=0.1), ...)
          {
            x <- as.matrix(x)
            y <- as.matrix(y)

            if(!(nrow(x)==nrow(y))) stop("Number of colums in x, y matrixes is not equall")
            if(!is(kernel,"kernel"))
              {
                if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
                kernel <- do.call(kernel, kpar)
              }
            if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")
            xpca <- kpca(x,kernel,...)
            ypca <- kpca(y,kernel,...)
            cca <- cancor(rotated(xpca), rotated(ypca))
            
            ret <- new("kcca")

            kcor(ret) <- cca$cor
            xcoef(ret) <- pcv(xpca) %*% cca$xcoef
            ycoef(ret) <- pcv(ypca) %*% cca$ycoef
            xvar(ret) <- rotated(xpca) %*% cca$xcoef
            yvar(ret) <- rotated(ypca) %*% cca$ycoef
            ret
          })
