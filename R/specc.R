

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

setMethod("specc",signature(x="matrix"),function(x, centers, kernel = "rbfdot", kpar = "local", iterations = 200, mod.sample =  1, na.action = na.omit, ...)
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

  
  if(is.character(kpar)) {
   kpar <- match.arg(kpar,c("automatic","local"))
   


    if(kpar == "automatic")
      {
        sam <- sample(1:m, floor(mod.sample*m))
        sx <- x[sam,]
        dota <- rowSums(sx*sx)/2
        ktmp <- crossprod(t(sx))
        for (i in 1:length(sam))
          ktmp[i,]<- 2*(-ktmp[i,] + dota + rep(dota[i],m))

        ## fix numerical prob.
        ktmp[ktmp<0] <- 0
        ktmp <- sqrt(ktmp)

        kmax <- max(ktmp)
        kmin <- min(ktmp + diag(rep(Inf,dim(ktmp)[1])))
        kmea <- mean(ktmp)
        lsmin <- log2(kmin)
        lsmax <- log2(kmax)
        midmax <- min(c(2*kmea, kmax))
        midmin <- max(c(kmea/2,kmin))
        rtmp <- c(seq(midmin,0.9*kmea,0.05*kmea), seq(kmea,midmax,0.08*kmea))
        if ((lsmax - (Re(log2(midmax))+0.5)) < 0.5) step <- (lsmax - (Re(log2(midmax))+0.5))
        else step <- 0.5
        if (((Re(log2(midmin))-0.5)-lsmin) < 0.5 ) stepm <-  ((Re(log2(midmin))-0.5) - lsmin)
        else stepm <- 0.5
        
        tmpsig <- c(2^(seq(lsmin,(Re(log2(midmin))-0.5), stepm)), rtmp, 2^(seq(Re(log2(midmax))+0.5, lsmax,step)))
        diss <- matrix(rep(Inf,length(tmpsig)*nc),ncol=nc)

        for (i in 1:length(tmpsig)){
          ka <- exp((-(ktmp^2))/(2*(tmpsig[i]^2)))
          diag(ka) <- 0
          
          d <- 1/sqrt(rowSums(ka))
     
          if(!any(d==Inf) && !any(is.na(d))&& (max(d)[1]-min(d)[1] < 10^4))
            {
              l <- d * ka %*% diag(d)
              xi <- eigen(l,symmetric=TRUE)$vectors[,1:nc]
              yi <- xi/sqrt(rowSums(xi^2))
              res <- kmeans(yi, centers, iterations)
              diss[i,] <- res$withinss
            }
        }

        ms <- which.min(rowSums(diss))
        kernel <- rbfdot((tmpsig[ms]^(-2))/2)
        ## Compute Affinity Matrix
        km <- kernelMatrix(kernel, x)
      }
 }
  if (kpar=="local")
    {
      s <- rep(0,m)
      dota <- rowSums(x*x)/2
      dis <- crossprod(t(x))
      for (i in 1:m)
        dis[i,]<- 2*(-dis[i,] + dota + rep(dota[i],m))



      ## fix numerical prob.
      dis[dis<0] <- 0
      
      for (i in 1:m)
        s[i] <- median(sort(sqrt(dis[i,]))[1:5])


      ## Compute Affinity Matrix
      km <- exp(-dis / s%*%t(s))
      kernel <- rbfdot(1)
    }
  else
    {
      if(!is(kernel,"kernel"))
        {
          if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
          kernel <- do.call(kernel, kpar)
        }
      if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

      ## Compute Affinity Matrix
      km <- kernelMatrix(kernel, x)
    }

    

  if(is(kernel)[1] == "rbfkernel")
    diag(km) <- 0
  
  d <- 1/sqrt(rowSums(km))
  l <- d * km %*% diag(d)
  xi <- eigen(l)$vectors[,1:nc]
  yi <- xi/sqrt(rowSums(xi^2))
  res <- kmeans(yi, centers, iterations)
  cent <- matrix(unlist(lapply(1:nc,ll<- function(l){colMeans(x[which(res$cluster==l),])})),ncol=dim(x)[2], byrow=TRUE)

  withss <- unlist(lapply(1:nc,ll<- function(l){sum((x[which(res$cluster==l),] - cent[l,])^2)}))
  
  return(new("specc", .Data=res$cluster, size = res$size, centers=cent, withinss=withss, kernelf= kernel))

})



setMethod("show","specc",
function(object){
 
  cat("Spectral Clustering object of class \"specc\"","\n")
  cat("\n","Cluster memberships:","\n","\n")
   cat(object@.Data,"\n","\n")
  show(kernelf(object))
  cat("\n")
  cat(paste("Centers: ","\n"))
  show(centers(object))
  cat("\n")
  cat(paste("Cluster size: ","\n"))
  show(size(object))
  cat("\n")
  cat(paste("Within-cluster sum of squares: ", "\n"))
  show(withinss(object))
  cat("\n")
})


  
