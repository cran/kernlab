setGeneric("ipop",function(c, H, A, b, l, u, r, sigf=7, maxiter=40, margin=0.05, bound=10, verb=0) standardGeneric("ipop"))
setMethod("ipop",signature(H="matrix"),
function(c, H, A, b, l, u, r, sigf=7, maxiter=40, margin=0.05, bound=10, verb=0)
  {

##ipop solves the quadratic programming problem
##minimize   c' * primal + 1/2 primal' * H * primal
##subject to b <= A*primal <= b + r
##           l <= x <= u
##           d is the optimizer itself
##returns primal and dual variables (i.e. x and the Lagrange
##multipliers for b <= A * primal <= b + r)
##for additional documentation see
##     R. Vanderbei
##     LOQO: an Interior Point Code for Quadratic Programming, 1992
## Author:      R version Alexandros Karatzoglou, orig. matlab Alex J. Smola
## Created:     12/12/97
## R Version:   12/08/03
## Updated:     13/08/03
## This code is released under the GNU Public License

    if(!is.matrix(H)) stop("H must be a matrix")
    if(!is.matrix(A)&&!is.vector(A)) stop("A must be a matrix or a vector")
    if(!is.matrix(c)&&!is.vector(c)) stop("c must be a matrix or a vector")
    if(!is.matrix(l)&&!is.vector(l)) stop("l must be a matrix or a vector")
    if(!is.matrix(u)&&!is.vector(u)) stop("u must be a matrix or a vector")

    n <- dim(H)[2]
    m <- dim(A)[1]
    primal<-rep(0,n)
    if (missing(b))
      bvec <- rep(0, m)
    if(n !=nrow(H))
      stop("H matrix is not symmetric")
    if (n != length(c))
      stop("H and c are incompatible!")
    if (n != ncol(A))
      stop("A and c are incompatible!")
    if (m != length(b))
      stop("A and b are incompatible!")
    if(n !=length(u))
      stop("u is incopatible with H")
    if(n !=length(l))
      stop("l is incopatible with H")

    m <- nrow(A)
    n <- ncol(A)
    H.diag <- diag(H)
    H.x <- H
    b.plus.1 <- max(svd(b)$d) + 1
    c.plus.1 <- max(svd(c)$d) + 1
    one.x <- -matrix(1,n,1)
    one.y <- -matrix(1,m,1)
                                        # starting point
    diag(H.x) <- H.diag + 1
    H.y <- diag(1,m)
    c.x <- c
    c.y <- b
  ## solve the system [-H.x A' A H.y] [x, y] = [c.x c.y]

    AP <- matrix(0,m+n,m+n)
    xp <- 1:(m+n) <= n
    AP[xp,xp] <- -H.x
    AP[xp == FALSE,xp] <- A
    AP[xp,xp == FALSE] <- t(A)
    AP[xp == FALSE, xp== FALSE] <- H.y
    s.tmp <- solve(AP,c(c.x,c.y))
    x<-s.tmp[1:n]
    y<-s.tmp[-(1:n)]

    g <- pmax(abs(x - l), bound)
    z <- pmax(abs(x), bound)
    t <- pmax(abs(u - x), bound)
    s <- pmax(abs(x), bound)
    v <- pmax(abs(y), bound)
    w <- pmax(abs(y), bound)
    p <- pmax(abs(r - w), bound)
    q <- pmax(abs(y), bound)
    mu <- as.vector(crossprod(z,g) + crossprod(v,w) + crossprod(s,t) + crossprod(p,q))/(2 * (m + n))
    sigfig <- 0
    counter <- 0
    alfa <- 1
    if (verb > 0)	                       # print at least one status report
      cat("Iter    PrimalInf  DualInf  SigFigs  Rescale  PrimalObj  DualObj","\n")

    while (counter < maxiter)
      {
                                       # #update the iteration counter
        counter <- counter + 1
                                        ##central path (predictor)
        H.dot.x <- H %*% x
        rho <- b - A %*% x + w
        nu <- l - x + g
        tau <- u - x - t
        alpha <- r - w - p
        sigma <- c  - crossprod(A, y) - z + s + H.dot.x
        beta <- y + q - v
        gamma.z <- - z
        gamma.w <- - w
        gamma.s <- - s
        gamma.q <- - q
         ## instrumentation
        x.dot.H.dot.x <-  crossprod(x, H.dot.x)
        primal.infeasibility <- max(svd(rbind(rho, tau, alpha, nu))$d) / b.plus.1
        dual.infeasibility <- max(svd(rbind(sigma,beta))$d) / c.plus.1
        primal.obj <- crossprod(c,x) + 0.5 * x.dot.H.dot.x
        dual.obj <- crossprod(b,y) - 0.5 * x.dot.H.dot.x + crossprod(l, z) - crossprod(u,s) - crossprod(r,q)
        old.sigfig <- sigfig
        sigfig <- max(-log10(abs(primal.obj - dual.obj)/(abs(primal.obj) + 1)), 0)
        if (sigfig >= sigf) break
        if (verb > 0)		      	# final report
          cat( counter, "\t", signif(primal.infeasibility,6), signif(dual.infeasibility,6), sigfig, alfa, primal.obj, dual.obj,"\n")
                                       ## some more intermediate variables (the hat section)
        hat.beta <- beta - v * gamma.w / w
        hat.alpha <- alpha - p * gamma.q / q
        hat.nu <- nu + g * gamma.z / z
        hat.tau <- tau - t * gamma.s / s
                                        ##the diagonal terms
        d <- z / g + s / t
        e <- 1 / (v / w + q / p)
                                        ## initialization before the big cholesky
        diag(H.x) <- H.diag + d
        diag(H.y) <- e
        c.x <- sigma - z * hat.nu / g - s * hat.tau / t
        c.y <- rho - e * (hat.beta - q * hat.alpha / p)
                                        ## and solve the system [-H.x A' A H.y] [delta.x, delta.y] <- [c.x c.y]
        AP[xp,xp] <- -H.x
        AP[xp == FALSE, xp== FALSE] <- H.y
        s1.tmp <- solve(AP,c(c.x,c.y))
        delta.x<-s1.tmp[1:n] ; delta.y <- s1.tmp[-(1:n)]
                                        ##backsubstitution
        delta.w <- - e * (hat.beta - q * hat.alpha / p + delta.y)
        delta.s <- s * (delta.x - hat.tau) / t
        delta.z <- z * (hat.nu - delta.x) / g
        delta.q <- q * (delta.w - hat.alpha) / p
        delta.v <- v * (gamma.w - delta.w) / w
        delta.p <- p * (gamma.q - delta.q) / q
        delta.g <- g * (gamma.z - delta.z) / z
        delta.t <- t * (gamma.s - delta.s) / s
                                        ##compute update step now (sebastian's trick)
        alfa <- - (1 - margin) / min(c(delta.g / g, delta.w / w, delta.t / t, delta.p / p, delta.z / z, delta.v / v, delta.s / s, delta.q / q, -1))
        newmu <- (crossprod(z,g) + crossprod(v,w) + crossprod(s,t) + crossprod(p,q))/(2 * (m + n))
        newmu <- mu * ((alfa - 1) / (alfa + 10))^2
        gamma.z <- mu / g - z - delta.z * delta.g / g
        gamma.w <- mu / v - w - delta.w * delta.v / v
        gamma.s <- mu / t - s - delta.s * delta.t / t
        gamma.q <- mu / p - q - delta.q * delta.p / p
                                        ## some more intermediate variables (the hat section)
        hat.beta <- beta - v * gamma.w / w
        hat.alpha <- alpha - p * gamma.q / q
        hat.nu <- nu + g * gamma.z / z
        hat.tau <- tau - t * gamma.s / s
                                        ## initialization before the big cholesky
                                        ##for (  i  in  1 : n H.x(i,i) <- H.diag(i) + d(i) ) {
                                        ##H.y <- diag(e)
        c.x <- sigma - z * hat.nu / g - s * hat.tau / t
        c.y <- rho - e * (hat.beta - q * hat.alpha / p)

                                        ## and solve the system [-H.x A' A H.y] [delta.x, delta.y] <- [c.x c.y]
        AP[xp,xp] <- -H.x
        AP[xp == FALSE, xp== FALSE] <- H.y
        s1.tmp <- solve(AP,c(c.x,c.y))
        delta.x<-s1.tmp[1:n] ; delta.y<-s1.tmp[-(1:n)]
                                              ##backsubstitution
        delta.w <- - e * (hat.beta - q * hat.alpha / p + delta.y)
        delta.s <- s * (delta.x - hat.tau) / t
        delta.z <- z * (hat.nu - delta.x) / g
        delta.q <- q * (delta.w - hat.alpha) / p
        delta.v <- v * (gamma.w - delta.w) / w
        delta.p <- p * (gamma.q - delta.q) / q
        delta.g <- g * (gamma.z - delta.z) / z
        delta.t <- t * (gamma.s - delta.s) / s
                                        ##compute the updates
        alfa <- - (1 - margin) / min(c(delta.g / g, delta.w / w, delta.t / t, delta.p / p, delta.z / z, delta.v / v, delta.s / s, delta.q / q, -1))
        x <- x + delta.x * alfa
        g <- g + delta.g * alfa
        w <- w + delta.w * alfa
        t <- t + delta.t * alfa
        p <- p + delta.p * alfa
        y <- y + delta.y * alfa
        z <- z + delta.z * alfa
        v <- v + delta.v * alfa
        s <- s + delta.s * alfa
        q <- q + delta.q * alfa
                                        ## these two lines put back in
        mu <- (crossprod(z,g) + crossprod(v,w) + crossprod(s,t) + crossprod(p,q))/(2 * (m + n))
        mu <- mu * ((alfa - 1) / (alfa + 10))^2
        mu <- newmu
      }
    if (verb > 0)		      	## final report
      cat( counter, primal.infeasibility, dual.infeasibility, sigfig, alfa, primal.obj, dual.obj)

    ret <- new("ipop")               ## repackage the results
    primal(ret) <- x
    dual(ret)   <- y
    if ((sigfig > sigf) & (counter < maxiter))
      how(ret)    <- 'converged'
    else
      {					## must have run out of counts
        if ((primal.infeasibility > 10e5) & (dual.infeasibility > 10e5))
          how(ret)    <- 'primal and dual infeasible'
        if (primal.infeasibility > 10e5)
          how(ret)    <- 'primal infeasible'
        if (dual.infeasibility > 10e5)
          how(ret)    <- 'dual infeasible'
        else					## don't really know
          how(ret)    <- 'slow convergence, change bound?'
      }
    ret
})
