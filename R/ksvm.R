setGeneric("ksvm", function(x, ...) standardGeneric("ksvm"))
setMethod("ksvm",signature(x="formula"),
function (x, data=NULL, ..., subset, na.action = na.omit, scaled = TRUE){
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- NULL
  m$formula <- m$x
  m$x <- NULL
  m$scaled <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  attr(Terms, "intercept") <- 0
  x <- model.matrix(Terms, m)
  y <- model.extract(m, response)
  if (length(scaled) == 1)
    scaled <- rep(scaled, ncol(x))
  if (any(scaled)) {
    remove <- unique(c(which(labels(Terms) %in% names(attr(x, "contrasts"))),
                       which(!scaled)
                       )
                     )
    scaled <- !attr(x, "assign") %in% remove
  }
   ret <- ksvm(x, y, scaled = scaled, ...)
  kcall(ret) <- call
  kterms(ret) <- Terms
  if (!is.null(attr(m, "na.action")))
    n.action(ret) <- attr(m, "na.action")
  return (ret)
})

setMethod("ksvm",signature(x="vector"),
function(x,...)
          { x <- t(t(x))
            ret <- ksvm(x, ...)
            return(ret)
          })
    
setMethod("ksvm",signature(x="matrix"),
function (x,
          y         = NULL,
          scaled    = TRUE,
          type      = NULL,
          kernel    = "rbfdot",
          kpar      = "automatic",
          C         = 1,
          nu        = 0.2,
          epsilon   = 0.1,
          prob.model = FALSE,
          class.weights = NULL,
          cross     = 0,
          fit       = TRUE,
          cache     = 40,
          tol       = 0.001,
          shrinking = TRUE,
          ...
          ,subset 
         ,na.action = na.omit)
{
  sparse  <- inherits(x, "matrix.csr")
  if (sparse) {
    if (!require(SparseM))
      stop("Need SparseM package for handling of sparse structures!")
  }
  ## subsetting and na-handling for matrices
  ret <- new("ksvm")
  if (!missing(subset)) x <- x[subset,]
  if (is.null(y))
    x <- na.action(x)
  else {
    df <- na.action(data.frame(y, x))
    y <- df[,1]
    x <- as.matrix(df[,-1])
  }
  n.action(ret) <- na.action
  
 if (is.null(type)) type(ret) <-
    if (is.null(y)) "one-svc"
    else if (is.factor(y)) "C-svc"
    else "eps-svr"
  
  if(!is.null(type))
  type(ret) <- match.arg(type,c("C-svc",
                                "nu-svc",
                                "kbb-svc",
                                "spoc-svc",
                                "C-bsvc",
                                "one-svc",
                                "eps-svr",
                                "eps-bsvr",
                                "nu-svr"))


  ## scaling, subsetting, and NA handling
  if (sparse) {
    scale <- rep(FALSE, ncol(x))
    if(!is.null(y)) na.fail(y)
    x <- t(t(x)) ## make shure that col-indices are sorted
  }
  unscaledx <- x  
  x.scale <- y.scale <- NULL
 ## scaling
  if (length(scaled) == 1)
    scaled <- rep(scaled, ncol(x))
  if (any(scaled)) {
    co <- !apply(x[,scaled, drop = FALSE], 2, var)
    if (any(co)) {
      scaled <- rep(FALSE, ncol(x))
      warning(paste("Variable(s)",
                    paste("`",colnames(x[,scaled, drop = FALSE])[co],
                          "'", sep="", collapse=" and "),
                    "constant. Cannot scale data.")
              )
    } else {
      xtmp <- scale(x[,scaled])
      x[,scaled] <- xtmp
      x.scale <- attributes(xtmp)[c("scaled:center","scaled:scale")]
      if (is.numeric(y)&&(type(ret)!="C-svc"&&type(ret)!="nu-svc"&&type(ret)!="C-bsvc"&&type(ret)!="spoc-svc"&&type(ret)!="kbb-svc")) {
        y <- scale(y)
        y.scale <- attributes(y)[c("scaled:center","scaled:scale")]
        y <- as.vector(y)
      }
      scaling(ret) <- list(scaled = scaled, x.scale = x.scale, y.scale = y.scale)
    }
  }
  ncols <- ncol(x)
  m <- nrows <- nrow(x)

  if(is.character(kernel)){
    kernel <- match.arg(kernel,c("rbfdot","polydot","tanhdot","vanilladot","laplacedot","besseldot","anovadot"))
  
  if((kernel == "tanhdot" || kernel == "vanilladot" || kernel == "polydot"|| kernel == "besseldot") &&  kpar=="automatic" )
    {
      cat (" Setting default kernel parameters ","\n")
      kpar <- list()
    }
  }
  
  if (!is.function(kernel))
  if (!is.list(kpar)&&is.character(kpar)&&(class(kernel)=="rbfkernel" || class(kernel) =="laplacedot" || kernel == "laplacedot"|| kernel=="rbfdot")){
    kp <- match.arg(kpar,"automatic")
    if(kp=="automatic")
      kpar <- list(sigma=sum(sigest(x,scaled=FALSE))/2)
   cat("Using automatic sigma estimation (sigest) for RBF or laplace kernel","\n")
   
  }
  if(!is(kernel,"kernel"))
    {
      if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
      kernel <- do.call(kernel, kpar)
    }

  if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

  if (!is.vector(y) && !is.factor (y) && !(type(ret)=="one-svc")) stop("y must be a vector or a factor.")
  if ((type(ret) != "one-svc") && nrows != nrow(x)) stop("x and y don't match.")
  if(nu > 1|| nu <0) stop("nu must be between 0 an 1.")
  
  weightlabels <- NULL
  nweights <- 0
  weight <- 0
  wl <- 0
  ## in case of classification: transform factors into integers
  if (type(ret) == "one-svc") # one class classification --> set dummy
    y <- 1
  else
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
    }
    else {
      if ((type(ret) =="C-svc" || type(ret) == "nu-svc" ||type(ret) == "C-bsvc" || type(ret) == "spoc-svc" || type(ret) == "kbb-svc") && any(as.integer (y) != y))
        stop ("dependent variable has to be of factor or integer type for classification mode.")

      if (type(ret) != "eps-svr" || type(ret) != "nu-svr"|| type(ret)!="eps-bsvr")
        lev(ret) <- sort(unique (y))
    }
 ## initialize    
  nclass(ret) <- length (unique(y))
  p <- 0
  K <- 0 
  svindex <- problem <- NULL
  sigma <- 0.1
  degree <- offset <- scale <- 1
  switch(is(kernel)[1],
         "rbfkernel" =
         {
           sigma <- kpar(kernel)$sigma
           ktype <- 2
         },
         "tanhkernel" =
         {
           sigma <- kpar(kernel)$scale
           offset <- kpar(kernel)$offset
           ktype <- 3
         },
         "polykernel" =
         {
           degree <- kpar(kernel)$degree
           sigma <- kpar(kernel)$scale
           offset <- kpar(kernel)$offset
           ktype <- 1
         },
         "vanillakernel" =
         {
           ktype <- 0
         },
	 "laplacekernel" =
	 {
	 ktype <- 5
	 sigma <- kpar(kernel)$sigma
	 },
         "besselkernel" =
         {
           ktype <- 6
           sigma <- kpar(kernel)$sigma
           degree <- kpar(kernel)$order
           offset <- kpar(kernel)$degree
         },
         "anovakernel" =
         {
           ktype <- 7
           sigma <- kpar(kernel)$sigma
           degree <-  kpar(kernel)$degree
         },
         {
           ktype <- 4
         }
         )
  prior(ret) <- list(NULL)


  if(type(ret) == "C-svc"){
    indexes <- lapply(sort(unique(y)), function(kk) which(y == kk))
    for (i in 1:(nclass(ret)-1)) {
      jj <- i+1
      for(j in jj:nclass(ret)) {
        p <- p+1
        ##prepare data
        li <- length(indexes[[i]])
        lj <- length(indexes[[j]])
        if(sparse)
          xd <- as.matrix.csr(0,(li+lj),dim(x)[2])
        else
          xd <- matrix(0,(li+lj),dim(x)[2])
        xdi <- 1:(li+lj) <= li
        xd[xdi,rep(TRUE,dim(x)[2])] <- x[indexes[[i]],]
        xd[xdi == FALSE,rep(TRUE,dim(x)[2])] <- x[indexes[[j]],]
        if(y[indexes[[i]][1]] < y[indexes[[j]]][1])
          {
            yd <- c(rep(1,li),rep(-1,lj))
            if(!is.null(class.weights)){
            weight <- class.weights[weightlabels[c(i,j)]]
            wl <- c(1,0)
            nweights <- 2
          }
          }
        else
          {
            yd <- c(rep(-1,li),rep(1,lj))
            if(!is.null(class.weights)){
            weight <- class.weights[weightlabels[c(j,i)]]
            wl <- c(0,1)
            nweigths <- 2
          }
          }
       
        boolabel <- yd >= 0
        prior1 <- sum(boolabel)

        md <- length(yd)
        prior0 <- md - prior1
        prior(ret)[[p]] <- list(prior1 = prior1, prior0 = prior0) 

        if(ktype==4)
          K <- kernelMatrix(kernel,xd)
        
        resv <- .Call("smo_optim",
                      as.double(t(xd)),
                      as.integer(nrow(xd)),
                      as.integer(ncol(xd)),
                      as.double(yd),
                      as.double(K),
                      
                      as.integer(if (sparse) xd@ia else 0),
                      as.integer(if (sparse) xd@ja else 0),
                      as.integer(sparse),
                      
                      as.double(matrix(rep(-1,m))), ##linear term
                      as.integer(ktype),
                      as.integer(0), 
                      as.double(C),
                      as.double(nu),
                      as.double(epsilon),
                      as.double(sigma),
                      as.double(degree),
                      as.double(offset),
                      as.integer(wl), ##weightlabel
                      as.double(weight),
                      as.integer(nweights),
                      as.double(cache), 
                      as.double(tol),
                      as.integer(shrinking),
                      PACKAGE="kernlab")
       
         alpha(ret)[p] <- list(resv[-(li+lj+1)])
        ## nonzero alpha*y
        coeff(ret)[p] <- list(alpha(ret)[[p]][alpha(ret)[[p]]>0]*yd[alpha(ret)[[p]]>0])
        ## store SV indexes from current problem for later use in predict
        alphaindex(ret)[p] <- list(c(indexes[[i]],indexes[[j]])[which(alpha(ret)[[p]]>0)])
        ## save the indexes from all the SV in a vector (use unique?)
        svindex <- c(svindex,alphaindex(ret)[[p]])
        ## store betas in a vector 
        b(ret) <- c(b(ret), resv[li+lj+1])
        ## used to reconstruct indexes for the patterns matrix x from "indexes" (really usefull ?)
        problem[p] <- list(c(i,j))
        ##store C  in return object
        param(ret)$C <- C
##        margin(ret)[p] <- (min(kernelMult(kernel,xd[1:li,],,alpha(ret)[[p]][1:li])) - max(kernelMult(kernel,xd[li:(li+lj),],,alpha(ret)[[p]][li:(li+lj)])))/2
      }
    }
  } 

if(type(ret) == "nu-svc"){
  indexes <- lapply(sort(unique(y)), function(kk) which(y == kk))
    for (i in 1:(nclass(ret)-1)) {
      jj <- i+1
      for(j in jj:nclass(ret)) {
        p <- p+1
       ##prepare data
        li <- length(indexes[[i]])
        lj <- length(indexes[[j]])
         if(sparse)
          xd <- as.matrix.csr(0,(li+lj),dim(x)[2])
        else
          xd <- matrix(0,(li+lj),dim(x)[2])
        xdi <- 1:(li+lj) <= li
        xd[xdi,rep(TRUE,dim(x)[2])] <- x[indexes[[i]],]
        xd[xdi == FALSE,rep(TRUE,dim(x)[2])] <- x[indexes[[j]],]
        if(y[indexes[[i]][1]] < y[indexes[[j]]][1])
          {
            yd <- c(rep(1,li),rep(-1,lj))
            if(!is.null(class.weights)){
            weight <- class.weights[weightlabels[c(i,j)]]
            wl <- c(1,0)
            nweights <- 2
          }
          }
        else
          {
            yd <- c(rep(-1,li),rep(1,lj))
            if(!is.null(class.weights)){
            weight <- class.weights[weightlabels[c(j,i)]]
            wl <- c(0,1)
            nweigths <- 2
          }
          }
       
        boolabel <- yd >= 0
        prior1 <- sum(boolabel)
        md <- length(yd)
        prior0 <- md - prior1
        prior(ret)[[p]] <- list(prior1 = prior1, prior0 = prior0)

        if(ktype==4)
             K <- kernelMatrix(kernel,xd)
        
        resv <- .Call("smo_optim",
                      as.double(t(xd)),
                      as.integer(nrow(xd)),
                      as.integer(ncol(xd)),
                      as.double(yd),
                      as.double(K),
                      as.integer(if (sparse) xd@ia else 0),
                      as.integer(if (sparse) xd@ja else 0),
                      as.integer(sparse),
                      
                      as.double(matrix(rep(-1,m))), #linear term
                      as.integer(ktype),
                      as.integer(1),
                      as.double(C),
                      as.double(nu),
                      as.double(epsilon),
                      as.double(sigma),
                      as.double(degree),
                      as.double(offset),
                      as.integer(wl), #weightlabl.
                      as.double(weight),
                      as.integer(nweights),
                      as.double(cache),
                      as.double(tol), 
                      as.integer(shrinking),
                      PACKAGE="kernlab")
        
        alpha(ret)[p] <- list(resv[-(li+lj+1)])
        ## alpha*y whithout including zeros (smaller kernel matrixes)
        coeff(ret)[p] <- list(alpha(ret)[[p]][alpha(ret)[[p]]!= 0])
        ##store SV indexes from current problem for later use in predict
         alphaindex(ret)[p] <- list(c(indexes[[i]],indexes[[j]])[which(resv[-(li+lj+1)]!= 0)])
        ##alphaindex(ret)[p] <- list(c(which(alpha(ret)[[p]][1:li]!=0)+li*(i-1),which(alpha(ret)[[p]][-(1:li)]!=0)+li+lj*(j-2)))
        ##save the indexes from all the SV in a vector (use unique!)
               svindex <- c(svindex,alphaindex(ret)[[p]])
        ## store betas in a vector 
               b(ret) <- c(b(ret), resv[li+lj+1])
        ## used to reconstruct indexes for the patterns matrix x from "indexes"
               problem[p] <- list(c(i,j))
               param(ret)$nu <- nu
        
##        margin(ret)[p] <- (min(kernelMult(kernel,xd[1:li,],,alpha(ret)[[p]][1:li])) - max(kernelMult(kernel,xd[li:(li+lj),],,alpha(ret)[[p]][li:(li+lj)])))/2
    
      }
    }
  } 
  


  if(type(ret) == "C-bsvc"){
     if(!is.null(class.weights))
     weightedC <- class.weights[weightlabels] * rep(C,nclass(ret))
    else
      weightedC <- rep(C,nclass(ret)) 
    indexes <- lapply(sort(unique(y)), function(kk) which(y == kk))
    for (i in 1:(nclass(ret)-1)) {
      jj <- i+1
      for(j in jj:nclass(ret)) {
        p <- p+1
        ##prepare data
        li <- length(indexes[[i]])
        lj <- length(indexes[[j]])
        if(sparse)
          xd <- as.matrix.csr(0,(li+lj),dim(x)[2])
        else
          xd <- matrix(0,(li+lj),dim(x)[2])
        xdi <- 1:(li+lj) <= li
        xd[xdi,rep(TRUE,dim(x)[2])] <- x[indexes[[i]],]
        xd[xdi == FALSE,rep(TRUE,dim(x)[2])] <- x[indexes[[j]],]
        if(y[indexes[[i]][1]] < y[indexes[[j]]][1])
          {
            yd <- c(rep(1,li),rep(-1,lj))
            if(!is.null(class.weights)){
            weight <- class.weights[weightlabels[c(i,j)]]
            wl <- c(1,0)
            nweights <- 2
          }
          }
        else
          {
            yd <- c(rep(-1,li),rep(1,lj))
            if(!is.null(class.weights)){
            weight <- class.weights[weightlabels[c(j,i)]]
            wl <- c(0,1)
            nweigths <- 2
          }
          }
       
        boolabel <- yd >= 0
        prior1 <- sum(boolabel)
        md <- length(yd)
        prior0 <- md - prior1
        prior(ret)[[p]] <- list(prior1 = prior1, prior0 = prior0) 

           if(ktype==4)
          K <- kernelMatrix(kernel,xd)
        
        resv <- .Call("tron_optim",
                      as.double(t(xd)),
                      as.integer(nrow(xd)),
                      as.integer(ncol(xd)),
                      as.double(yd),
                      as.double(K),
                      as.integer(if (sparse) xd@ia else 0),
                      as.integer(if (sparse) xd@ja else 0),
                      as.integer(sparse),
                      as.integer(2),
                      as.double(0), ##countc
                      as.integer(ktype),
                      as.integer(5), 
                      as.double(C),
                      as.double(epsilon),
                      as.double(sigma),
                      as.double(degree),
                      as.double(offset),
                      as.double(1),  ##  cost value of alpha seeding
                      as.double(2),  ## step value of alpha seeding
                      as.integer(wl), ##weightlabel
                      as.double(weight),
                      as.integer(nweights),
                      as.double(weightedC),
                      as.double(cache), 
                      as.double(tol),
                      as.integer(10), ##qpsize
                      as.integer(shrinking),
                      PACKAGE="kernlab")
        alpha(ret)[p] <- list(resv)
        ## nonzero alpha*y
        coeff(ret)[p] <- list(alpha(ret)[[p]][alpha(ret)[[p]]>0]*yd[alpha(ret)[[p]]>0])
        ## store SV indexes from current problem for later use in predict
        alphaindex(ret)[p] <- list(c(indexes[[i]],indexes[[j]])[which(resv>0)])
        ## save the indexes from all the SV in a vector (use unique?)
        svindex <- c(svindex,alphaindex(ret)[[p]])
        ## store betas in a vector 
        b(ret) <- - sapply(coeff(ret),sum) 
        ## used to reconstruct indexes for the patterns matrix x from "indexes" (really usefull ?)
        problem[p] <- list(c(i,j))
        ##store C  in return object
        param(ret)$C <- C
##        margin(ret)[p] <- (min(kernelMult(kernel,xd[1:li,],,alpha(ret)[[p]][1:li])) - max(kernelMult(kernel,xd[li:(li+lj),],,alpha(ret)[[p]][li:(li+lj)])))/2
      }
    }
  } 

if(type(ret) =="spoc-svc")
  {
    if(!is.null(class.weights))
     weightedC <- class.weights[weightlabels] * rep(C,nclass(ret))
    else
      weightedC <- rep(C,nclass(ret)) 
    yd <- sort(y,method="quick", index.return = TRUE)
    x<-x[yd$ix,]
    count <- 0

       if(ktype==4)
          K <- kernelMatrix(kernel,x)
    
    resv <- .Call("tron_optim",
                  as.double(t(x)),
                  as.integer(nrow(x)),
                  as.integer(ncol(x)),
                  as.double(rep(yd$x-1,2)),
                  as.double(K),
                  as.integer(if (sparse) x@ia else 0),
                  as.integer(if (sparse) x@ja else 0),
                  as.integer(sparse),
                  as.integer(nclass(ret)),
                  as.integer(count),
                  as.integer(ktype),
                  as.integer(7), 
                  as.double(C),
                  as.double(epsilon),
                  as.double(sigma),
                  as.double(degree),
                  as.double(offset),
                  as.double(C), 
                  as.double(2), #Cstep
                  as.integer(0), #weightlabel
                  as.double(0),
                  as.integer(0),
                  as.double(weightedC),
                  as.double(cache), 
                  as.double(tol),
                  as.integer(10), #qpsize
                  as.integer(shrinking),
                  PACKAGE="kernlab")

    alpha(ret) <- t(matrix(resv,nclass(ret)))
    coeff(ret) <- lapply(1:nclass(ret), function(x) alpha(ret)[,x][alpha(ret)[,x]!=0])
    names(coeff(ret)) <- lev(ret)
    alphaindex(ret) <-  lapply(1:nclass(ret), function(x) which(alpha(ret)[,x]!=0))
    names(alphaindex(ret)) <- lev(ret)
    svindex <- which(alpha(ret)!=0)
    b(ret) <- 0
    param(ret)$C <- C
  }

if(type(ret) =="kbb-svc")
  {
    if(!is.null(class.weights))
      weightedC <- weightlabels * rep(C,nclass(ret))
    else
      weightedC <- rep(C,nclass(ret)) 
    yd <- sort(y,method="quick", index.return = TRUE)
    x<-x[yd$ix,]
    count <-  sapply(unique(yd$x), function(c) length(yd$x[yd$x==c]))

       if(ktype==4)
          K <- kernelMatrix(kernel,x)
    resv <- .Call("tron_optim",
                  as.double(t(x)),
                  as.integer(nrow(x)),
                  as.integer(ncol(x)),
                  as.double(yd$x-1),
                  as.double(K),
                  as.integer(if (sparse) x@ia else 0),
                  as.integer(if (sparse) x@ja else 0),
                  as.integer(sparse),
                  as.integer(nclass(ret)),
                  as.integer(count),
                  as.integer(ktype),
                  as.integer(8),
                  as.double(C),
                  as.double(epsilon),
                  as.double(sigma),
                  as.double(degree),
                  as.double(offset),
                  as.double(1), #Cbegin
                  as.double(2), #Cstep
                  as.integer(0), #weightlabl.
                  as.double(0),
                  as.integer(0),
                  as.double(weightedC),
                  as.double(cache),
                  as.double(tol),
                  as.integer(10), #qpsize
                  as.integer(shrinking),
                  PACKAGE="kernlab")
    alpha(ret) <- matrix(resv,nrow(x))

    coeff(ret) <-  lapply(1:(nclass(ret)-1), function(x) alpha(ret)[,x][alpha(ret)[,x]!=0])
    alphaindex(ret) <-  lapply(1:(nclass(ret)-1), function(x) which(alpha(ret)[,x]!=0))
    svindex <- which(resv !=0)  ## have to figure out what to do with this...!
    b(ret) <- - sapply(coeff(ret),sum) 
  }

  if(type(ret) =="one-svc")
  {
       if(ktype==4)
          K <- kernelMatrix(kernel,x)
       
    resv <- .Call("smo_optim",
                  as.double(t(x)),
                  as.integer(nrow(x)),
                  as.integer(ncol(x)),
                  as.double(matrix(rep(1,m))),
                  as.double(K),
                  as.integer(if (sparse) x@ia else 0),
                  as.integer(if (sparse) x@ja else 0),
                  as.integer(sparse),
                  as.double(matrix(rep(-1,m))),
                  as.integer(ktype),
                  as.integer(2),
                  as.double(C),
                  as.double(nu),
                  as.double(epsilon),
                  as.double(sigma),
                  as.double(degree),
                  as.double(offset),
                  as.integer(0), #weightlabl.
                  as.double(0),
                  as.integer(0),
                  as.double(cache),
                  as.double(tol),
                  as.integer(shrinking),
                  PACKAGE="kernlab")
    alpha(ret) <- resv[-(m+1)]
    coeff(ret) <- alpha(ret)[alpha(ret)!=0]
    alphaindex(ret) <- which(alpha(ret)!=0) ## in this case and in regr. the same with svindex
    svindex <- which(alpha(ret) !=0) 
    b(ret) <- resv[(m+1)]
    param(ret)$nu <- nu
  }

  if(type(ret) =="eps-svr")
  {
       if(ktype==4)
          K <- kernelMatrix(kernel,x)
       
     resv <- .Call("smo_optim",
                  as.double(t(x)),
                  as.integer(nrow(x)),
                  as.integer(ncol(x)),
                  as.double(y),
                   as.double(K),
                  as.integer(if (sparse) x@ia else 0),
                  as.integer(if (sparse) x@ja else 0),
                  as.integer(sparse),
                  as.double(matrix(rep(-1,m))),
                  as.integer(ktype),
                  as.integer(3),
                  as.double(C),
                  as.double(nu),
                  as.double(epsilon),
                  as.double(sigma),
                  as.double(degree),
                  as.double(offset),
                  as.integer(0), #weightlabl.
                  as.double(0),
                  as.integer(0),
                  as.double(cache), 
                  as.double(tol), 
                  as.integer(shrinking), 
                   PACKAGE="kernlab")
    alpha(ret) <- resv[-(m+1)]
    coeff(ret) <- alpha(ret)[alpha(ret)!=0]
    alphaindex(ret) <- which(alpha(ret)!=0)
    svindex <- which(alpha(ret) !=0) 
    b(ret) <- resv[(m+1)]
    param(ret)$epsilon <- epsilon
  }

 if(type(ret) =="nu-svr")
  {
       if(ktype==4)
          K <- kernelMatrix(kernel,x)
    
    resv <- .Call("smo_optim",
                  as.double(t(x)),
                  as.integer(nrow(x)),
                  as.integer(ncol(x)),
                  as.double(y),
                  as.double(K),
                  as.integer(if (sparse) x@ia else 0),
                  as.integer(if (sparse) x@ja else 0),
                  as.integer(sparse),
                  as.double(matrix(rep(-1,m))),
                  as.integer(ktype),
                  as.integer(4),
                  as.double(C),
                  as.double(nu),
                  as.double(epsilon),
                  as.double(sigma),
                  as.integer(degree),
                  as.double(offset),
                  as.integer(0),
                  as.double(0),
                  as.integer(0),
                  as.double(cache), 
                  as.double(tol), 
                  as.integer(shrinking), 
                  PACKAGE="kernlab")
    alpha(ret) <- resv[-(m+1)]
    coeff(ret) <- alpha(ret)[alpha(ret)!=0]
    alphaindex(ret) <- which(alpha(ret)!=0)
    svindex <- which(alpha(ret) !=0) 
    b(ret) <- resv[(m+1)]
    param(ret)$epsilon <- epsilon
    param(ret)$nu <- nu
  }

   if(type(ret) =="eps-bsvr")
  {
       if(ktype==4)
          K <- kernelMatrix(kernel,x)
       
     resv <- .Call("tron_optim",
                  as.double(t(x)),
                  as.integer(nrow(x)),
                  as.integer(ncol(x)),
                  as.double(y),
                   as.double(K),
                  as.integer(if (sparse) x@ia else 0),
                  as.integer(if (sparse) x@ja else 0),
                  as.integer(sparse),
                  as.double(2),
                   as.double(0),
                  as.integer(ktype),
                  as.integer(6),
                   as.double(C),
                   as.double(epsilon),
                   as.double(sigma),
                  as.double(degree),
                  as.double(offset),
                   as.double(1),  #Cbegin
                   as.double(2), #Cstep
                  as.integer(0), #weightlabl.
                  as.double(0),
                  as.integer(0),
                   as.double(0),
                  as.double(cache), 
                  as.double(tol),
                   as.integer(10), #qpsize
                  as.integer(shrinking), 
                   PACKAGE="kernlab")
    alpha(ret) <- resv
    coeff(ret) <- alpha(ret)[alpha(ret)!=0]
    alphaindex(ret) <- which(alpha(ret)!=0)
    svindex <- which(alpha(ret) !=0) 
    b(ret) <- -sum(alpha(ret))
    param(ret)$epsilon <- epsilon
  }

  
  kcall(ret) <- match.call()
  kernelf(ret) <- kernel
  ## param(ret) <- list(C=C, nu = nu, epsilon = epsilon)
  xmatrix(ret) <- x
  ymatrix(ret) <- y
  SVindex(ret) <- unique(svindex)
  nSV(ret)  <- length(unique(svindex))
  if(nSV(ret)==0)
    stop("No Support Vectors found. You may want to change your parameters")
  fit(ret)  <- if (fit)
    predict(ret, unscaledx) else NA

  if (fit){
    if(type(ret)=="C-svc"||type(ret)=="nu-svc"||type(ret)=="spoc-svc"||type(ret)=="kbb-svc"||type(ret)=="C-bsvc")
      error(ret) <- 1 - .classAgreement(table(y,as.integer(fit(ret))))
    if(type(ret)=="eps-svr"||type(ret)=="nu-svr"||type(ret)=="eps-bsvr")
      error(ret) <- drop(crossprod(fit(ret) - y)/m)
  }

  cross(ret) <- -1
  if(cross == 1)
    cat("\n","cross should be >1 no cross-validation done!","\n","\n")
  else if (cross > 1)
    {
     
      cerror <- 0
      suppressWarnings(vgr<-split(sample(1:m,m),1:cross))
      for(i in 1:cross)
        {
         
          cind <-  unsplit(vgr[-i],factor(rep((1:cross)[-i],unlist(lapply(vgr[-i],length)))))
          if(type(ret)=="C-svc"||type(ret)=="nu-svc"||type(ret)=="spoc-svc"||type(ret)=="kbb-svc"||type(ret)=="C-bsvc")
            {
              if(is.null(class.weights))
                cret <- ksvm(x[cind,],y[cind],type = type(ret),kernel=kernel,kpar = NULL, C=C, nu=nu, tol=tol, scaled=FALSE, cross = 0, fit = FALSE ,cache = cache)
              else
                cret <- ksvm(x[cind,],lev(ret)[y[cind]],type = type(ret),kernel=kernel,kpar = NULL, C=C, nu=nu, tol=tol, scaled=FALSE, cross = 0, fit = FALSE, class.weights = class.weights,cache = cache)
               cres <- predict(cret, x[vgr[[i]],])
            cerror <- (1 - .classAgreement(table(y[vgr[[i]]],as.integer(cres))))/cross + cerror
            }
          if(type(ret)=="eps-svr"||type(ret)=="nu-svr"||type(ret)=="eps-bsvr")
            {
              cret <- ksvm(x[cind,],y[cind],type=type(ret),kernel=kernel,kpar = NULL,C=C,nu=nu,epsilon=epsilon,tol=tol,scaled=FALSE, cross = 0, fit = FALSE, cache = cache, prob.model = FALSE)
              cres <- predict(cret, x[vgr[[i]],])
              cerror <- drop(crossprod(cres - y[vgr[[i]]])/m)/cross + cerror
            }
         }
      cross(ret) <- cerror
    }

  prob.model(ret) <- list(NULL)
  
  if(prob.model)
    {

      pres <- NULL
      if(type(ret)=="C-svc"||type(ret)=="nu-svc"||type=="C-bsvc")
        {
          p <- 0
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
                {
                  yd <- c(rep(1,li),rep(-1,lj))
                  if(!is.null(class.weights)){
                    weight <- weightlabels[c(i,j)]
                    wl <- c(1,0)
                    nweights <- 2
                  }
                }
              else
                {
                  yd <- c(rep(-1,li),rep(1,lj))
                  if(!is.null(class.weights)){
                    weight <- weightlabels[c(j,i)]
                    wl <- c(0,1)
                    nweigths <- 2
                  }
                }
              m <- li+lj
              suppressWarnings(vgr<-split(sample(1:m,m),1:3))
              for(k in 1:3)
                {
                  cind <- unsplit(vgr[-k],factor(rep((1:3)[-k],unlist(lapply(vgr[-k],length)))))
                  cret <- ksvm(xd[cind,], yd[cind], type = type(ret),kernel=kernel,kpar = NULL, C=C, nu=nu, tol=tol, scaled=FALSE, cross = 0, fit = FALSE,cache = cache, prob.model=FALSE)
                  pres <- rbind(pres,predict(cret, xd[vgr[[k]],],type="decision"))
                  
                }
              prob.model(ret)[[p]] <- .probPlatt(pres)
            }
          }
        }
      if(type(ret) == "eps-svr"||type(ret) == "nu-svr"||type(ret)=="eps-bsvr"){
        suppressWarnings(vgr<-split(sample(1:m,m),1:3))
        for(i in 1:3)
          {
            cind <- unsplit(vgr[-i],factor(rep((1:cross)[-i],unlist(lapply(vgr[-i],length)))))

            cret <- ksvm(x[cind,],y[cind],type=type(ret),kernel=kernel,kpar = NULL,C=C,nu=nu,epsilon=epsilon,tol=tol,scaled=FALSE, cross = 0, fit = FALSE, cache = cache, prob.model = FALSE)
            cres <- predict(cret, x[vgr[[i]],])
            pres <- rbind(pres,predict(cret, x[vgr[[i]],],type="decision"))
          }
        pres[abs(pres) > (5*sd(pres))] <- 0
        prob.model(ret) <- list(sum(abs(pres))/dim(pres)[1])
      }
    }
 
  xmatrix(ret) <- x
  ## loss(ret) <- sum((1 - y * fitted(ret))[(1 - y * fitted(ret))>0]/m)
  return(ret)
})

.classAgreement <- function (tab) {
  n <- sum(tab)
  if (!is.null(dimnames(tab))) {
    lev <- intersect(colnames(tab), rownames(tab))
    p0 <- sum(diag(tab[lev, lev])) / n
  } else {
    m <- min(dim(tab))
    p0 <- sum(diag(tab[1:m, 1:m])) / n
  }
  return(p0)
}


#**************************************************************#


setMethod("predict", signature(object = "ksvm"),
function (object, newdata, type = "response", coupler = "minpair")
{
  sc <- 0
  type <- match.arg(type,c("response","probabilities","votes","decision"))
  if (missing(newdata) && type!="response")
    return(fit(object))
  else if(missing(newdata))
    {
      newdata <- xmatrix(object)
      sc <- 1
    }
  
  ncols <- ncol(xmatrix(object))
  nrows <- nrow(xmatrix(object))
  oldco <- ncols

  if (!is.null(kterms(object)))
    {
      if(!is.matrix(newdata))
        newdata <- model.matrix(delete.response(kterms(object)), as.data.frame(newdata), na.action = n.action(object))
    }
  else
    newdata  <- if (is.vector(newdata)) t(t(newdata)) else as.matrix(newdata)

  newcols  <- 0
  newnrows <- nrow(newdata)
  newncols <- ncol(newdata)
  newco    <- newncols
    
  if (oldco != newco) stop ("test vector does not match model !")
  p<-0

  if (is.list(scaling(object)) && sc != 1)
    newdata[,scaling(object)$scaled] <-
      scale(newdata[,scaling(object)$scaled, drop = FALSE],
            center = scaling(object)$x.scale$"scaled:center",
            scale  = scaling(object)$x.scale$"scaled:scale"
            )
 
  if(type == "response" || type =="decision" || type=="votes")
    {
  if(type(object)=="C-svc"||type(object)=="nu-svc"||type(object)=="C-bsvc")
    {
      predres <- 1:newnrows
      votematrix <- matrix(0,nclass(object),newnrows)
      for(i in 1:(nclass(object)-1))
        {
        jj <- i+1
        for(j in jj:nclass(object))
          {
            p <- p+1
            ret <- kernelMult(kernelf(object),newdata,as.matrix(xmatrix(object)[alphaindex(object)[[p]],]),coeff(object)[[p]]) - b(object)[p]
            if(type=="decision")
              votematrix[p,] <- ret
            else{
              votematrix[i,ret>0] <- votematrix[i,ret>0] + 1
              votematrix[j,ret<0] <- votematrix[j,ret<0] + 1
            }
          }
      }
      if(type == "decision")
        {
          if (nclass(object) == 2)
            predres <- t(votematrix)[,1,drop = FALSE]
          else
            predres <-  t(votematrix)
        }
      else 
        predres <- sapply(predres, function(x) which.max(votematrix[,x]))
    }

  if(type(object) == "spoc-svc")
    {
      predres <- 1:newnrows
      votematrix <- matrix(0,nclass(object),newnrows)
      for(i in 1:nclass(object))
        votematrix[i,] <- kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[i]],],coeff(object)[[i]])
      predres <- sapply(predres, function(x) which.max(votematrix[,x]))
    }

  if(type(object) == "kbb-svc")
    { 
      predres <- 1:newnrows
      votematrix <- matrix(0,nclass(object),newnrows)
      se <- 1:(nclass(object)-1)
      A <- rowSums(alpha(object))
      Aind <- unique(unlist(alphaindex(object)))
      sA <- sum(A)
      Aindex <- which(A!=0)
      for(i in 1:nclass(object))
        {
          for(k in se[se <= i])
            votematrix[k,] <- votematrix[k,] - kernelMult(kernelf(object), newdata, xmatrix(object)[alphaindex(object)[[k]],],coeff(object)[[k]]) + b(object)[k]

          votematrix[i,] <- votematrix[i,] + kernelMult(kernelf(object), newdata, xmatrix(object)[Aind,],A[Aind]) + sA

          for(k in se[!se<i])
            votematrix[k+1,] <- votematrix[k+1,] - kernelMult(kernelf(object), newdata, xmatrix(object)[alphaindex(object)[[k]],],coeff(object)[[k]]) + b(object)[k]
        }
      
      predres <- sapply(predres, function(x) which.max(votematrix[,x]))
    }
}

  if(type == "probabilities")
    {
      if(is.null(prob.model(object)[[1]]))
        stop("ksvm object contains no probability model. Make sure you set the paramater prob.model in ksvm during training.")
      
      if(type(object)=="C-svc"||type(object)=="nu-svc"||type(object)=="C-bsvc")
        {
          binprob <- matrix(0, newnrows, nclass(object)*(nclass(object) - 1)/2)
          for(i in 1:(nclass(object)-1))
            {
              jj <- i+1
              for(j in jj:nclass(object))
                {
                  p <- p+1
                  binprob[,p] <- .SigmoidPredict(as.vector(kernelMult(kernelf(object),newdata,as.matrix(xmatrix(object)[alphaindex(object)[[p]],]),coeff(object)[[p]]) - b(object)[p]), prob.model(object)[[p]]$A, prob.model(object)[[p]]$B)
                }
            }
          multiprob <- couple(binprob, coupler = coupler)
        }
      else
        stop("probability estimates only supported for C-svc and nu-svc")
    }
  
  if(type(object) == "one-svc")
    {
      ret <- kernelMult(kernelf(object),newdata,as.matrix(xmatrix(object)[alphaindex(object),]),coeff(object)) - b(object)
      ret[ret>0]<-1
      ##one-class-classification: return TRUE/FALSE (probabilities ?)
      return(ret == 1)      
    }
  else {
    if(type(object)=="eps-svr"||type(object)=="nu-svr"||type(object)=="eps-bsvr")
      {
        predres <- kernelMult(kernelf(object),newdata,as.matrix(xmatrix(object)[alphaindex(object),]),coeff(object)) - b(object)
      }
    else {
      ##classification & votes : return votematrix
      if(type == "votes")
        return(votematrix)
      
      ##classification & probabilities : return probability matrix
      if(type == "probabilities")
        {
          colnames(multiprob) <- lev(object)
          return(multiprob)
        }

      if(is.numeric(lev(object)) && type == "response")
         return(lev(object)[predres])
      
      if (is.character(lev(object)) && type!="decision")
        {
          ##classification & type response: return factors
          if(type == "response")
            return(factor (lev(object)[predres], levels = lev(object)))
        }
    }
  }
 
  if (!is.null(scaling(object)$y.scale))
    ## return raw values, possibly scaled back
    return(predres * scaling(object)$y.scale$"scaled:scale" + scaling(object)$y.scale$"scaled:center")
  else
    ##else: return raw values
    return(predres)
})

#****************************************************************************************#

setMethod("show","ksvm",
function(object){
  cat("Support Vector Machine object of class \"ksvm\"","\n")
  cat("\n")
   cat(paste("SV type:", type(object)))
  
  switch(type(object),
         "C-svc" = cat(paste("  (classification)", "\n")),
         "nu-svc" = cat(paste("  (classification)", "\n")),
         "C-bsvc" = cat(paste("  (classification)", "\n")),
         "one-svc" = cat(paste("  (novelty detection)", "\n")),
         "spoc-svc" = cat(paste("  (classification)", "\n")),
         "kbb-svc" = cat(paste("  (classification)", "\n")),
         "eps-svr" = cat(paste("  (regression)","\n")),
         "nu-svr" = cat(paste("  (regression)","\n"))
         )
  
  switch(type(object),
         "C-svc" = cat(paste(" parameter : cost C =",param(object)$C, "\n")),
         "nu-svc" = cat(paste(" parameter : nu =", param(object)$nu, "\n")),
         "C-bsvc" = cat(paste(" parameter : cost C =",param(object)$C, "\n")),
         "one-svc" = cat(paste(" parameter : nu =", param(object)$nu, "\n")),
         "spoc-svc" = cat(paste(" parameter : cost C =",param(object)$C, "\n")),
         "kbb-svc" = cat(paste(" parameter : cost C =",param(object)$C, "\n")),
         "eps-svr" = cat(paste(" parameter : epsilon =",param(object)$epsilon,"\n")),
         "nu-svr" = cat(paste(" parameter : epsilon =", param(object)$epsilon, "nu =", param(object)$nu,"\n"))
         )
  cat("\n")
 show(kernelf(object))
  cat(paste("\nNumber of Support Vectors :", nSV(object),"\n"))


##    if(type(object)=="C-svc" || type(object) == "nu-svc")
##      cat(paste("Margin width :",margin(object),"\n"))
  if(!is.null(fit(object)))
    cat(paste("Training error :", round(error(object),6),"\n"))
  if(cross(object)!= -1)
    cat("Cross validation error :",round(cross(object),6),"\n")
  if(!is.null(prob.model(object)[[1]])&&(type(object)=="eps-svr" ||type(object)=="nu-svr"||type(object)=="eps-bsvr"))
    cat("Laplace distr. width :",round(prob.model(object)[[1]],6),"\n")
  ##train error & loss
})


setMethod("plot", signature(x = "ksvm", y ="missing"),
function(x, data = NULL, grid = 50, slice = list(), ...) {

  if (type(x) =="C-svc" || type(x) == "nu-svc") {
    if(nclass(x) > 2)
      stop("plot function only supports binary classification")
  
    if (!is.null(kterms(x))&&!is.null(data))
      {
        if(!is.matrix(data))
          sub <- model.matrix(delete.response(kterms(x)), as.data.frame(data), na.action = n.action(x))
      }
    else if(!is.null(data))
      sub <-  as.matrix(data)
    else
      sub <- xmatrix(x)
      
    sub <- sub[,!colnames(xmatrix(x))%in%names(slice)]
    xr <- seq(min(sub[,2]), max(sub[,2]), length = grid)
    yr <- seq(min(sub[,1]), max(sub[,1]), length = grid)
    sc <- 0
    if(is.null(data))
      {
        sc  <- 1
        data <- xmatrix(x)
      }
    if(is.data.frame(data) || !is.null(kterms(x))){
      lis <- c(list(yr), list(xr), slice)
      names(lis)[1:2] <- colnames(sub)
      new <- expand.grid(lis)[,labels(kterms(x))]
    }
    else
      new <- expand.grid(xr,yr)
    
    if(sc== 1) 
      scaling(x) <- NULL

    preds <- predict(x, new ,type = "decision")
    
    if(is.null(kterms(x)))
      xylb <- colnames(sub)
    else
      xylb <- names(lis)
    lvl <- 37
    
    mymax <- max(abs(preds))
    mylevels <- pretty(c(0, mymax), 15)
    ##Z$ mycols <- rev(heat.colors(length(mylevels)-1))
    ##Z# mycols <- c(rev(mycols), mycols)
    nl <- length(mylevels)-2
    mycols <- c(hsv(0, (nl:0/nl)^1.3, 1), hsv(2/3, (0:nl/nl)^1.3, 1))
    mylevels <- c(-rev(mylevels[-1]), mylevels)

    index <- max(which(mylevels < min(preds))):min(which(mylevels > max(preds)))
    mycols <- mycols[index]
    mylevels <- mylevels[index]
    
    
    filled.contour(xr, yr, matrix(as.numeric(preds), nr = length(xr), byrow = TRUE), col = mycols, levels = mylevels
                    ,plot.axes = {
                     axis(1)
                     axis(2)
                     if(!is.null(data)){
                       points(sub[-SVindex(x),2], sub[-SVindex(x),1],col= (ymatrix(x)[-SVindex(x)]+3))
                       points(sub[SVindex(x),2], sub[SVindex(x),1], pch = "x",col=(ymatrix(x)[SVindex(x)]+3))}
                     else{
                       points(sub[-SVindex(x),],col=(ymatrix(x)[-SVindex(x)]+3))
                       points(sub[SVindex(x),], pch ="x",col=(ymatrix(x)[SVindex(x)]+3))
                     }},
                   nlevels = lvl,
                   plot.title = title(main = "SVM classification plot", xlab = xylb[2], ylab = xylb[1]),
                   ...
                   )
  } else {
    stop("Only plots of classification ksvm objects supported")
  }
})


setGeneric(".probPlatt", function(deci) standardGeneric(".probPlatt"))
setMethod(".probPlatt",signature(deci="ANY"),
function(deci)
  {
    if (is.matrix(deci))
    deci <- as.vector(deci)
    if (!is.vector(deci))
      stop("input should be matrix or vector")
    ## Create label and count priors
    boolabel <- deci >= 0
    prior1 <- sum(boolabel)
    m <- length(deci)
    prior0 <- m - prior1

    ## set parameters (should be on the interface I guess)
    maxiter <- 100
    minstep <- 1e-10
    sigma <- 1e-3
    eps <- 1e-5

   
    
    ## Construct target support
    hiTarget <- (prior1 + 1)/(prior1 + 2)
    loTarget <- 1/(prior0 + 2)
    length <- prior1 + prior0
    t <- rep(loTarget, m)
    t[boolabel] <- hiTarget

    ##Initial Point & Initial Fun Value
    A <- 0
    B <- log((prior0 + 1)/(prior1 + 1))
    fval <- 0

    fApB <- deci*A + B
    bindex <- fApB >= 0

    fval <- sum(t[bindex]*fApB[bindex] + log(1 + exp(-fApB[bindex])))
    fval <- fval + sum((t[!bindex] - 1)*fApB[!bindex] + log(1+exp(fApB[!bindex])))


    for (it in 1:maxiter)
      {
        h11 <- h22 <- sigma
        h21 <- g1 <- g2 <- 0
        fApB <- deci*A + B
        
        for (l in 1:2)
          {
          if (l == 1)
            {
              bindex <- fApB >= 0
              p <- exp(-fApB[bindex])/(1 + exp(-fApB[bindex]))
              q <- 1/(1+exp(-fApB[bindex]))
            }
          if (l == 2)
            {
              bindex <- fApB < 0
              p <- 1/(1 + exp(fApB[bindex]))
              q <- exp(fApB[bindex])/(1 + exp(fApB[bindex]))
            }
          
          d2 <- p*q
          h11 <- h11 + sum(d2*deci[bindex]^2)
          h22 <- h22 + sum(d2)
          h21 <- h21 + sum(deci[bindex]*d2)
          d1 <- t[bindex] - p
          g1 <- g1 + sum(deci[bindex]*d1)
          g2 <- g2 + sum(d1)
        }
        ## Stopping Criteria
        if (abs(g1) < eps && abs(g2) < eps)
          break

        ## Finding Newton Direction -inv(t(H))%*%g
        det <- h11*h22 - h21^2
        dA <- -(h22*g1 - h21*g2) / det
        dB <- -(-h21*g1 + h11*g2) / det
        gd <- g1*dA + g2*dB

        ## Line Search
        stepsize <- 1

        while(stepsize >= minstep)
          {
            newA <- A + stepsize * dA
            newB <- B + stepsize * dB

            ## New function value
            newf <- 0
            fApB <- deci * newA + newB
            bindex <- fApB >= 0 
            newf <- sum(t[bindex] * fApB[bindex] + log(1 + exp(-fApB[bindex])))
            newf <- newf + sum((t[!bindex] - 1)*fApB[!bindex] + log(1 + exp(fApB[!bindex])))

           ## Check decrease
            if (newf < (fval + 0.0001 * stepsize * gd))
              {
                A <- newA
                B <- newB
                fval <- newf
                break
              }
            else
              stepsize <- stepsize/2
          }
        if (stepsize < minstep)
          {
            cat("line search fails", A, B, g1, g2, dA, dB, gd)
            ret <- .SigmoidPredict(deci, A, B)
            return(ret)
          }
      }
    if(it >= maxiter -1)
      cat("maximum number of iterations reached",g1,g2)

    ret <- list(A=A, B=B)
    return(ret)
  })

 ## Sigmoid predict function

.SigmoidPredict <- function(deci, A, B)
  {
fApB <- deci*A +B
k<-length(deci)
ret <- rep(0,k)
bindex <- fApB >= 0
ret[bindex] <- exp(-fApB[bindex])/(1 + exp(-fApB[bindex]))
ret[!bindex] <- 1/(1 + exp(fApB[!bindex]))
return(ret)
}

