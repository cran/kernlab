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
  {
    x <- t(t(x))
    ret <- ksvm(x, ...)
    return(ret)
  })
    
setMethod("ksvm",signature(x="matrix"),
function (x,
          y         = NULL,
          scaled    = TRUE,
          type      = NULL,
          kernel    = "rbfdot",
          kpar      = list(sigma = 0.1),
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
    if (is.null(y)) "one-classification"
    else if (is.factor(y)) "C-classification"
    else "eps-regression"
  
  if(!is.null(type))
  type(ret) <- match.arg(type,c("C-classification",
                          "nu-classification",
                         "spoc-classification",
                          "kbb-classification",
                          "one-classification",
                          "eps-regression",
                          "nu-regression"))


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
      if (is.numeric(y)&&(type(ret)!="C-classification"&&type(ret)!="nu-classification"&&type(ret)!="spoc-classification")) {
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
  
  if((kernel == "tanhdot" || kernel == "vanilladot" || kernel == "polydot") && any(names(kpar) == "sigma") )
    {
      cat("no sigma parameter in this kernel. Did you set a hyperparameter for this kernel ?")
      cat ("\n", " Setting default kernel parameters")
      kpar <- NULL
    }
  }
  if (!is.list(kpar)&&is.character(kpar)&&( class(kernel)=="rbfkernel" || kernel=="rbfdot")){
    kp <- match.arg(kpar,"automatic")
    if(kp=="automatic")
      kpar <- list(sigma=sum(sigest(x,scaled=FALSE))/2)
   cat("Using automatic sigma estimation (sigest) for RBF kernel","\n")
  }
  if(!is(kernel,"kernel"))
    {
      if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
      kernel <- do.call(kernel, kpar)
    }

  if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")



  ##internal function used by smo_optim                                   
  .kernelin <- function(v)
    {
      return(kernelMatrix(kernel,xd,xd[v,,drop=FALSE]))
    }

  if (!is.vector(y) && !is.factor (y) && !(type(ret)=="one-classification")) stop("y must be a vector or a factor.")
  if ((type(ret) != "one-classification") && nrows != nrow(x)) stop("x and y don't match.")
  if(nu > 1|| nu <0) stop("nu must be between 0 an 1.")
  
  weightlabels <- NULL
  nweights <- 0
  weight <- 0
  wl <- 0
  ## in case of classification: transform factors into integers
  if (type(ret) == "one-classification") # one class classification --> set dummy
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
      if ((type(ret) =="C-classification" || type(ret) == "nu-classification" || type(ret) == "spoc-classification" || type(ret) == "kbb-classification") && any(as.integer (y) != y))
        stop ("dependent variable has to be of factor or integer type for classification mode.")

      if (type(ret) != "eps-regression" || type(ret) != "nu-regression")
        lev(ret) <- unique (y)
    }
 ## initialize    
  nclass(ret) <- length (lev(ret))
  p <- 0
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
         ktype <- 4
         )
  prior(ret) <- list(NULL)
  if(type(ret) == "C-classification"){
    indexes <- lapply(sort(unique(y)), function(kk) which(y == kk))
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
       
        boolabel <- yd >= 0
        prior1 <- sum(boolabel)
        md <- length(yd)
        prior0 <- md - prior1
        prior(ret)[[p]] <- list(prior1 = prior1, prior0 = prior0) 
        
        resv <- .Call("smo_optim",
                      as.double(t(xd)),
                      as.integer(nrow(xd)),
                      as.integer(ncol(xd)),
                      as.double(yd),

                      
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
                      .kernelin,
                      environment(.kernelin), PACKAGE="kernlab")
        alpha(ret)[p] <- list(resv[-(li+lj+1)])
        ## nonzero alpha*y
        coeff(ret)[p] <- list(alpha(ret)[[p]][alpha(ret)[[p]]>0]*yd[alpha(ret)[[p]]>0])
        ## store SV indexes from current problem for later use in predict
        alphaindex(ret)[p] <- list(c(indexes[[i]],indexes[[j]])[which(resv[-(li+lj+1)]>0)])
        ## save the indexes from all the SV in a vector (use unique?)
        svindex <- c(svindex,alphaindex(ret)[[p]])
        ## store betas in a vector 
        b(ret) <- c(b(ret), resv[li+lj+1])
        ## used to reconstruct indexes for the patterns matrix x from "indexes" (really usefull ?)
        problem[p] <- list(c(i,j))
        ##store C  in return object
        param(ret)$C <- C
      }
    }
  } 

if(type(ret) == "nu-classification"){
  indexes <- lapply(sort(unique(y)), function(kk) which(y == kk))
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

        boolabel <- yd >= 0
        prior1 <- sum(boolabel)
        md <- length(yd)
        prior0 <- md - prior1
        prior(ret)[[p]] <- list(prior1 = prior1, prior0 = prior0)
        
        resv <- .Call("smo_optim",
                      as.double(t(xd)),
                      as.integer(nrow(xd)),
                      as.integer(ncol(xd)),
                      as.double(yd),
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
                      as.integer(0), #weightlabl.
                      as.double(0),
                      as.integer(0),
                      as.double(cache),
                      as.double(tol), 
                      as.integer(shrinking),
                      .kernelin,
                      environment(.kernelin), PACKAGE="kernlab")
        
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
      }
    }
  } 
  

if(type(ret) =="spoc-classification")
  {
    if(!is.null(class.weights))
     weightedC <- weightlabels * rep(C,nclass(ret))
    else
      weightedC <- rep(C,nclass(ret)) 
    yd <- sort(y,method="quick", index.return = TRUE)
    x<-x[yd$ix,]
    count <- 0
    resv <- .Call("tron_optim",
                  as.double(t(x)),
                  as.integer(nrow(x)),
                  as.integer(ncol(x)),
                  as.double(rep(yd$x-1,2)),
                  as.integer(if (sparse) x@ia else 0),
                  as.integer(if (sparse) x@ja else 0),
                  as.integer(sparse),
                  as.integer(nclass(ret)),
                  as.integer(count),
                  as.integer(ktype),
                  as.integer(2), 
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
                  .kernelin,
                  environment(.kernelin), PACKAGE="kernlab")

    alpha(ret) <- t(matrix(resv,nclass(ret)))
    coeff(ret) <- lapply(1:nclass(ret), function(x) alpha(ret)[,x][alpha(ret)[,x]!=0])
    names(coeff(ret)) <- lev(ret)
    alphaindex(ret) <-  lapply(1:nclass(ret), function(x) which(alpha(ret)[,x]!=0))
    names(alphaindex(ret)) <- lev(ret)
    svindex <- which(alpha(ret)!=0)
    b(ret) <- 0
    param(ret)$C <- C
  }

if(type(ret) =="kbb-classification")
  {
    if(!is.null(class.weights))
      weightedC <- weightlabels * rep(C,nclass(ret))
    else
      weightedC <- rep(C,nclass(ret)) 
    yd <- sort(y,method="quick", index.return = TRUE)
    x<-x[yd$ix,]
    count <-  sapply(unique(yd$x), function(c) length(yd$x[yd$x==c]))
    resv <- .Call("tron_optim",
                  as.double(t(x)),
                  as.integer(nrow(x)),
                  as.integer(ncol(x)),
                  as.double(yd$x-1),
                  as.integer(if (sparse) x@ia else 0),
                  as.integer(if (sparse) x@ja else 0),
                  as.integer(sparse),
                  as.integer(nclass(ret)),
                  as.integer(count),
                  as.integer(ktype),
                  as.integer(1),
                  as.double(C),
                  as.double(epsilon),
                  as.double(sigma),
                  as.double(degree),
                  as.double(offset),
                  as.double(C),
                  as.double(2), #Cstep
                  as.integer(0), #weightlabl.
                  as.double(0),
                  as.integer(0),
                  as.double(weightedC),
                  as.double(cache),
                  as.double(tol),
                  as.integer(10), #qpsize
                  as.integer(shrinking),
                  .kernelin,
                  environment(.kernelin), PACKAGE="kernlab")
    start <-rep(0,nclass(ret))
    start[1]<-0
    start2<-rep(0,nclass(ret))
    alpha(ret)<-matrix(0,nrow(x),nclass(ret)-1)
    start[2:nclass(ret)]<-cumsum(count)[1:nclass(ret)-1]
    for (i in 2:nclass(ret))
      start2[i] <- start2[i-1] + nrow(x) - count[i]
    p<-1
    se<-1:nclass(ret)
    for(i in se){
      for(j in (start[i]+1):(start[i]+count[i]))
        {
          for(k in se[se<i])
            alpha(ret)[p,k] <- resv[start2[k]+j-count[k]]
          for(k in se[se>i])
            alpha(ret)[p,k-1] <- resv[start2[k]+j]
          p <- p+1
        }
    }

    coeff(ret) <-  lapply(1:(nclass(ret)-1), function(x) alpha(ret)[,x][alpha(ret)[,x]!=0])
    alphaindex(ret) <-  lapply(1:(nclass(ret)-1), function(x) which(alpha(ret)[,x]!=0))
    svindex <- which(resv !=0)  ## have to figure out what to do with this...!
    b(ret) <- 0
  }

  if(type(ret) =="one-classification")
  {
    resv <- .Call("smo_optim",
                  as.double(t(x)),
                  as.integer(nrow(x)),
                  as.integer(ncol(x)),
                  as.double(matrix(rep(1,m))),
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
                  .kernelin,
                  environment(.kernelin), PACKAGE="kernlab")
    alpha(ret) <- resv[-(m+1)]
    coeff(ret) <- alpha(ret)[alpha(ret)!=0]
    alphaindex(ret) <- which(alpha(ret)!=0) ## in this case and in regr. the same with svindex
    svindex <- which(alpha(ret) !=0) 
    b(ret) <- resv[(m+1)]
    param(ret)$nu <- nu
  }

  if(type(ret) =="eps-regression")
  {
     resv <- .Call("smo_optim",
                  as.double(t(x)),
                  as.integer(nrow(x)),
                  as.integer(ncol(x)),
                  as.double(y),
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
                  .kernelin,
                  environment(.kernelin),PACKAGE="kernlab")
    alpha(ret) <- resv[-(m+1)]
    coeff(ret) <- alpha(ret)[alpha(ret)!=0]
    alphaindex(ret) <- which(alpha(ret)!=0)
    svindex <- which(alpha(ret) !=0) 
    b(ret) <- resv[(m+1)]
    param(ret)$epsilon <- epsilon
  }

 if(type(ret) =="nu-regression")
  {
    resv <- .Call("smo_optim",
                  as.double(t(x)),
                  as.integer(nrow(x)),
                  as.integer(ncol(x)),
                  as.double(y),
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
                  .kernelin,
                  environment(.kernelin), PACKAGE="kernlab")
    alpha(ret) <- resv[-(m+1)]
    coeff(ret) <- alpha(ret)[alpha(ret)!=0]
    alphaindex(ret) <- which(alpha(ret)!=0)
    svindex <- which(alpha(ret) !=0) 
    b(ret) <- resv[(m+1)]
    param(ret)$epsilon <- epsilon
    param(ret)$nu <- nu
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
    if(type(ret)=="C-classification"||type(ret)=="nu-classification"||type(ret)=="spoc-classification"||type(ret)=="kbb-classification")
      error(ret) <- 1 - .classAgreement(table(y,as.integer(fit(ret))))
    if(type(ret)=="eps-regression"||type(ret)=="nu-regression")
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
           cind <- unsplit(vgr[-i],1:(m-length(vgr[[i]])))
          if(type(ret)=="C-classification"||type(ret)=="nu-classification"||type(ret)=="spoc-classification"||type(ret)=="kbb-classification")
            { 
              cret <- ksvm(x[cind,],factor (lev(ret)[y[cind]], levels = lev(ret)),type = type(ret),kernel=kernel,kpar = NULL, C=C, nu=nu, tol=tol, scaled=FALSE, cross = 0, fit = FALSE, class.weights = class.weights,cache = cache)
               cres <- predict(cret, x[vgr[[i]],])
            cerror <- (1 - .classAgreement(table(y[vgr[[i]]],as.integer(cres))))/cross + cerror
            }
          if(type(ret)=="eps-regression"||type(ret)=="nu-regression")
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
          pres <- NULL
          m <- li+lj
          suppressWarnings(vgr<-split(sample(1:m,m),1:3))
          for(k in 1:3)
            {
              cind <- unsplit(vgr[-k],1:(m-length(vgr[[k]])))
              if(type(ret)=="C-classification"||type(ret)=="nu-classification")
                { 
                  cret <- ksvm(xd[cind,], yd[cind], type = type(ret),kernel=kernel,kpar = NULL, C=C, nu=nu, tol=tol, scaled=FALSE, cross = 0, fit = FALSE,cache = cache, prob.model=FALSE)
                  pres <- rbind(pres,predict(cret, xd[vgr[[k]],],type="decision"))
                }
            }
          prob.model(ret)[[p]] <- .probPlatt(pres)
        }
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
  if (missing(newdata))
    return(fit(object))
  ncols <- ncol(xmatrix(object))
  nrows <- nrow(xmatrix(object))
  oldco <- ncols

  if (!is.null(kterms(object)))
    {  
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

  if (is.list(scaling(object)))
    newdata[,scaling(object)$scaled] <-
      scale(newdata[,scaling(object)$scaled, drop = FALSE],
            center = scaling(object)$x.scale$"scaled:center",
            scale  = scaling(object)$x.scale$"scaled:scale"
            )

  type <- match.arg(type,c("response","probabilities","votes","decision"))
 
  if(type == "response" || type =="decision")
    {
  if(type(object)=="C-classification"||type(object)=="nu-classification")
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
        predres <- t(votematrix)
      else 
        predres <- sapply(predres, function(x) which.max(votematrix[,x]))
    }

  if(type(object) == "spoc-classification")
    {
      predres <- 1:newnrows
      votematrix <- matrix(0,nclass(object),newnrows)
      for(i in 1:nclass(object))
        votematrix[i,] <- kernelMult(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[i]],],coeff(object)[[i]])
      predres <- sapply(predres, function(x) which.max(votematrix[,x]))
    }

  if(type(object) == "kbb-classification")
    {
      predres <- 1:newnrows
      votematrix <- matrix(0,nclass(object),newnrows)
      se <- 1:(nclass(object)-1)
      A <- rowSums(alpha(object))
      Aindex <- which(A!=0)
      for(i in 1:nclass(object))
        {
          for(k in se[se<i])
            votematrix[k,] <- votematrix[k,] - (kernelMatrix(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[k]],]) +1)%*%coeff(object)[[k]]

          votematrix[i,] <- votematrix[i,] + (kernelMatrix(kernelf(object),newdata,xmatrix(object)[Aindex,])+1)%*%A[Aindex]      

          for(k in se[!se<i])
            votematrix[k+1,] <- votematrix[k+1,] - (kernelMatrix(kernelf(object),newdata,xmatrix(object)[alphaindex(object)[[k]],]) +1)%*%coeff(object)[[k]]
        }
      
      predres <- sapply(predres, function(x) which.max(votematrix[,x]))
    }
}

  if(type == "probabilities")
    {
      if(is.null(prob.model(object)[[1]]))
        stop("ksvm object contains no probability model. Make sure you set the paramater prob.model in ksvm during training.")
      
      if(type(object)=="C-classification"||type(object)=="nu-classification")
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
        stop("probability estimates only supported for C-classification and nu-classification")
    }
  
  if(type(object) == "one-classification")
    {
      ret <- kernelMult(kernelf(object),newdata,as.matrix(xmatrix(object)[alphaindex(object),]),coeff(object)) - b(object)
      ret[ret>0]<-1
    }
  
  if(type(object)=="eps-regression"||type(object)=="nu-regression")
    {
      predres <- kernelMult(kernelf(object),newdata,as.matrix(xmatrix(object)[alphaindex(object),]),coeff(object)) - b(object)
    }
  
  if (is.character(lev(object)) && type!="decision")
    {
      ##classification & probabilities : return probability matrix
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
  else if (!is.null(scaling(object)$y.scale))
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
   cat(paste("SV type:", type(object), "\n"))
  switch(type(object),
         "C-classification" = cat(paste(" parameter : cost C =",param(object)$C, "\n")),
         "nu-classification" = cat(paste(" parameter : nu =", param(object)$nu, "\n")),
         "one-classification" = cat(paste(" parameter : nu =", param(object)$nu, "\n")),
         "spoc-classification" = cat(paste(" parameter : cost C =",param(object)$C, "\n")),
         "eps-regression" = cat(paste(" parameter : epsilon =",param(object)$epsilon,"\n")),
         "nu-regression" = cat(paste(" parameter : epsilon =", param(object)$epsilon, "nu =", param(object)$nu,"\n"))
         )
  cat("\n")
 show(kernelf(object))
  cat(paste("\nNumber of Support Vectors :", nSV(object),"\n"))
  if(!is.null(fit(object)))
  cat(paste("Training error :", round(error(object),6),"\n"))
  if(cross(object)!= -1)
    cat("Cross validation error :",round(cross(object),6),"\n")
  ##train error & loss
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

