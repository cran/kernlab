

#kkmeans function

setGeneric("kkmeans",function(x, ...) standardGeneric("kkmeans"))
setMethod("kkmeans", signature(x = "formula"),
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
    res <- kkmeans(x, ...)
   
    cl[[1]] <- as.name("kkmeans")
    kcall(res) <- cl
    if(!is.null(na.act)) 
        n.action(res) <- na.action
   
    return(res)
  })

setMethod("kkmeans",signature(x="matrix"),function(x, centers, kernel
                                = "rbfdot", kpar = list(sigma=0.1),
                                alg ="kkmeans", p = 1,
                               max.iter = 200, mod.sample =  1, na.action = na.omit, ...)
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

  if(!is(kernel,"kernel"))
    {
      if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
      kernel <- do.call(kernel, kpar)
    }
  if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")


  if(length(centers) == 1){
    suppressWarnings(vgr<- vgr2 <- split(sample(1:m,m),1:centers))
    ncenters <- centers
  }
  else
    ncenters <- dim(centers)[1]


  if(is.character(kpar)) 
    alg <- match.arg(alg,c("kkmeans","kerninghan", "normcut"))

    if(alg == "kkmeans")
      {
        p <- NULL
        D <- NULL
        D1 <- NULL
        w <- rep(1,m)
      }
    if(alg=="kerninghan")
      {
        p <- p
        D <- kernelMult(kernel,x, , rep(1,m))
        w <- rep(1,m)
        D1 <- NULL
      }
    if(alg=="normcut")
      {
        p <- p
        D1 <- 1
        w <- kernelMult(kernel,x, , rep(1,m))
      }
    
    
  
  a <- rowSums(x^2)
  dismat <- matrix(0,m,ncenters)
  
  kdiag <- drop(kernel(1,1))
  dc  <- secsum <- rep(1,ncenters)
  mindis <- rep(0,m)
  cind <- 1:ncenters

    for ( i in 1:ncenters)
        {
          ## compute last sum eq. 1 once per cluster 
          ## for(k in 1:length(vgr[[i]]))
          ##     {
          ##       c <- rowSums(x[vgr[[i]],]^2)
          ##       secsum <- secsum + sum(kernelFast(kernel,x[vgr[[i]],,drop=FALSE],x[vgr[[i]][k],,drop=FALSE], c)*w[vgr[[i]]]*w[vgr[[i]][k]])/sum(w[[vgr[[i]]]])^2
          ## }
          
          secsum[i] <- sum(affinMult(kernel, x[vgr[[i]],,drop=FALSE],,
          w[vgr[[i]]], p , D, D1) * w[vgr[[i]]])/sum(w[vgr[[i]]])^2
          
        ## compute second sum eq. 1
          dismat[,i] <- - 2 *
          affinMult(kernel,x,x[vgr[[i]],,drop=FALSE], w[vgr[[i]]], p ,
          D, D1)/sum(w[vgr[[i]]]) + secsum[i] + kdiag 
            
          ##        for (j in 1:m)
          ##        {
          ##          dismat[j,i] <-
          ##        2*sum(kernelFast(kernel,x[vgr[[i]],,drop=FALSE],
          ##        x[j,,drop=FALSE], a) * w[[vgr[[i]]]])/sum(w[[vgr[[i]]]]) +
          ##        secsum
          ##        }
       
         
      }

  cluserm <- max.col(-dismat)
  for(i in 1:ncenters)
    vgr2[[i]] <- which(cluserm==i)
    
  while(sum(dc) < 0.00001){

    for(h in 1:ncenters)
          dc[h] <- -2*sum(affinMult(kernel,x[vgr2[[h]],,drop=FALSE],x[vgr[[h]],,drop=FALSE],w[vgr[[h]]], p, D, D1)*w[vgr2[[h]]])/(sum(w[vgr[[h]]])*sum(w[vgr2[[h]]]))+sum(affinMult(kernel, x[vgr[[h]],,drop=FALSE],,w[vgr[[h]]],p,D,D1) * w[vgr[[h]]])/sum(w[[vgr[[h]]]])^2 +sum(affinMult(kernel, x[vgr2[[h]],,drop=FALSE], ,w[vgr2[[h]]], p , D, D1)* w[vgr2[[h]]])/sum(w[vgr2[[h]]])^2

    vgr <- vgr2
    dismat <- t(t(dismat) - dc)
    cl2 <- max.col(-dismat)

    for (u in 1:ncenters){
      secsum[u] <- sum(affinMult(kernel, x[vgr[[u]],,drop=FALSE], ,w[vgr[[u]]], p, D, D1) * w[vgr[[u]]])/sum(w[[vgr[[u]]]])^2
      mindis[cl2==u] <- - 2*affinMult(kernel,x[cl2==u,,drop=FALSE],x[vgr[[u]],,drop=FALSE], w[vgr[[u]]],p,D,D1)/sum(w[[vgr[[u]]]]) + secsum[u] + kdiag
  }

    dismat[,cl2] <- mindis
    cl3 <- max.col(-dismat)
    recal <- cl3!=cl2

     for (u in 1:ncenters){
       dismat[recal,cint!=cl2] <- - 2 *affinMult(kernel,x[recal[cl2!=u],,drop=FALSE],x[vgr[[u]],,drop=FALSE], w[vgr[[u]]],p,D,D1)/sum(w[[vgr[[u]]]]) + secsum[u] + kdiag
             
     }
  }

cluster <- max.col(-dismat)
size <- unlist(lapply(1:ncenters, ll <- function(l){length(which(cluster==l))}))
cent <- matrix(unlist(lapply(1:ncenters,ll<- function(l){colMeans(x[which(cluster==l),])})),ncol=dim(x)[2], byrow=TRUE)
withss <- unlist(lapply(1:ncenters,ll<- function(l){sum((x[which(cluster==l),] - cent[l,])^2)}))
  
  return(new("specc", .Data=cluster, size = size, centers=cent, withinss=withss, kernelf= kernel))
  
})




setGeneric("affinMult",function(kernel, x, y = NULL, z, p, D, D1, blocksize = 256) standardGeneric("affinMult"))

affinMult.rbfkernel <- function(kernel, x, y=NULL, z, p, D, D1,blocksize = 256)
{
  if(is.null(p)&is.null(D)&is.null(D1))
    res <- kernelMult(kernel,x,y,z)
  else{
  if(!is.matrix(y)&&!is.null(y)) stop("y must be a matrix")
  if(!is.matrix(z)&&!is.vector(z)) stop("z must be a matrix or a vector")
  sigma <- kpar(kernel)$sigma
  n <- dim(x)[1]
  m <- dim(x)[2]
  nblocks <- floor(n/blocksize)
  lowerl <- 1
  upperl <- 0
  dota <- as.matrix(rowSums(x^2))

  if (is.null(y) & is.null(D1))
    {
      if(is.vector(z))
        {
          if(!length(z) == n) stop("vector z length must be equal to x rows")
          z <- matrix(z,n,1)
        }
      if(!dim(z)[1]==n)
        stop("z rows must equal x rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      if(nblocks > 0)
        {
          dotab <- rep(1,blocksize)%*%t(dota)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            res[lowerl:upperl,] <- exp(sigma*(2*x[lowerl:upperl,]%*%t(x) - dotab - dota[lowerl:upperl]%*%t(rep.int(1,n))))%*%z - z[lowerl:upperl,]*(1-p)
            lowerl <- upperl + 1
          
          }
      }
      if(lowerl <= n)
        res[lowerl:n,] <- exp(sigma*(2*x[lowerl:n,]%*%t(x) - rep.int(1,n+1-lowerl)%*%t(dota) - dota[lowerl:n]%*%t(rep.int(1,n))))%*%z- z[lowerl:upperl,]*(1-p)

    }
  if(is.matrix(y) & is.null(D1))
    {
      n2 <- dim(y)[1]
      if(is.vector(z))
        {
          if(!length(z) == n2) stop("vector z length must be equal to y rows")
          z <- matrix(z,n2,1)
        }
      if(!dim(z)[1]==n2)
        stop("z length must equal y rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      dotb <- as.matrix(rowSums(y*y))
      
       if(nblocks > 0)
         {
           dotbb <-  rep(1,blocksize)%*%t(dotb)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            if(upperl < n2)
            res[lowerl:upperl,] <- exp(sigma*(2*x[lowerl:upperl,]%*%t(y) - dotbb - dota[lowerl:upperl]%*%t(rep.int(1,n2))))%*%z-z[lowerl:upperl,]*(1-p) - z[lowerl:upperl,]*D[lowerl:upperl] 
            if(upperl >n2 & lowerl <n2){
              res[lowerl:upperl,] <- exp(sigma*(2*x[lowerl:upperl,]%*%t(y) - dotbb - dota[lowerl:upperl]%*%t(rep.int(1,n2))))%*%z
              res[lowerl:n2,] <- res[lowerl:n2,] - z[lowerl:n2,]*(1-p) - z[lowerl:n2,]*D[lowerl:n2]
            }
            else
              res[lowerl:upperl,] <- exp(sigma*(2*x[lowerl:upperl,]%*%t(y) - dotbb - dota[lowerl:upperl]%*%t(rep.int(1,n2))))%*%z
            lowerl <- upperl + 1
          }
      }
      if(lowerl <= n){
        if(lowerl >n2 & n>=n2){
          res[lowerl:n,] <- exp(sigma*(2*x[lowerl:n,]%*%t(y) - rep.int(1,n+1-lowerl)%*%t(dotb) -dota[lowerl:n]%*%t(rep.int(1,n2))))%*%z
          res[lowerl:n2,] <- res[lowerl:n2,] - z[lowerl:n2,]*(1-p) - z[lowerl:n2,]*D[lowerl:n2]
          
        }
      else
          res[lowerl:n,] <- exp(sigma*(2*x[lowerl:n,]%*%t(y) - rep.int(1,n+1-lowerl)%*%t(dotb) - dota[lowerl:n]%*%t(rep.int(1,n2))))%*%z
      }
    }
   
   if (is.null(y) & !is.null(D1))
    {
      if(is.vector(z))
        {
          if(!length(z) == n) stop("vector z length must be equal to x rows")
          z <- matrix(z,n,1)
        }
      if(!dim(z)[1]==n)
        stop("z rows must equal x rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      if(nblocks > 0)
        {
          dotab <- rep(1,blocksize)%*%t(dota)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            tmp <- exp(sigma*(2*x[lowerl:upperl,]%*%t(x) - dotab - dota[lowerl:upperl]%*%t(rep.int(1,n))))
            D1 <- 1/colSums(tmp)
            res[lowerl:upperl,] <- D1*tmp%*%diag(D1)%*%z - z[lowerl:upperl,]*(1-D1)
            lowerl <- upperl + 1
          
          }
      }
      if(lowerl <= n){
        tmp <- exp(sigma*(2*x[lowerl:n,]%*%t(x) - rep.int(1,n+1-lowerl)%*%t(dota) - dota[lowerl:n]%*%t(rep.int(1,n))))
        res[lowerl:n,] <- D1*tmp%*%diag(D1)%*%z- z[lowerl:upperl,]*(1-D1)
      }
    }
  if(is.matrix(y) &!is.null(D1))
    {
      n2 <- dim(y)[1]
      if(is.vector(z))
        {
          if(!length(z) == n2) stop("vector z length must be equal to y rows")
          z <- matrix(z,n2,1)
        }
      if(!dim(z)[1]==n2)
        stop("z length must equal y rows")
      res <- matrix(rep(0,dim(z)[2]*n), ncol = dim(z)[2])
      dotb <- as.matrix(rowSums(y*y))
      ones <- rep(1,blocksize)
       if(nblocks > 0)
         {
           dotbb <-  rep(1,blocksize)%*%t(dotb)
        for(i in 1:nblocks)
          {
            upperl = upperl + blocksize
            if(upperl < n2)
            tmp <-   exp(sigma*(2*x[lowerl:upperl,]%*%t(y) - dotbb - dota[lowerl:upperl]%*%t(rep.int(1,n2))))
            D1 <- 1/colSums(tmp)
            res[lowerl:upperl,] <- D1*tmp%*%diag(D1)%*%z-z[lowerl:upperl,]*(1-D1)  
            if(upperl >n2 & lowerl <n2){
              tmp <- exp(sigma*(2*x[lowerl:upperl,]%*%t(y) - dotbb -dota[lowerl:upperl]%*%t(rep.int(1,n2))))
              D1 <- 1/colSums(tmp)
              res[lowerl:upperl,] <- D1*tmp%*%diag(D1)%*%z
              res[lowerl:n2,] <- res[lowerl:n2,] - z[lowerl:n2,]*(1-D1)
            }
            else{
              tmp <- exp(sigma*(2*x[lowerl:upperl,]%*%t(y) - dotbb - dota[lowerl:upperl]%*%t(rep.int(1,n2))))
              D1 <- 1/colSums(tmp)
              res[lowerl:upperl,] <- D1*tmp%*%Diag(D1)%*%z
            }
            lowerl <- upperl + 1
          }
      }
      if(lowerl <= n){
        if(lowerl >n2 & n>=n2){
          tmp <- exp(sigma*(2*x[lowerl:n,]%*%t(y) -rep.int(1,n+1-lowerl)%*%t(dotb) -dota[lowerl:n]%*%t(rep.int(1,n2))))
          D1 <- 1/colSums(tmp)
          res[lowerl:n,] <- D1*tmp%*%diag(D1)%*%z
          res[lowerl:n2,] <- res[lowerl:n2,] - z[lowerl:n2,]*(1-D1)
          
        }
        else{
          tmp <- exp(sigma*(2*x[lowerl:n,]%*%t(y) -rep.int(1,n+1-lowerl)%*%t(dotb) -dota[lowerl:n]%*%t(rep.int(1,n2))))
          D1 <- 1/colSums(tmp)
          res[lowerl:n,] <- D1*tmp%*%diag(D1)%*%z
        }
      }
    }
}
   
  return(res)
}  
setMethod("affinMult",signature(kernel="kernel", x="matrix"),affinMult.rbfkernel)


  
