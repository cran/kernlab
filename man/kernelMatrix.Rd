\name{kernelMatrix}
\alias{kernelMatrix}
\alias{kernelMult}
\alias{kernelPol}
\alias{kernelFast}
\alias{kernelPol,kernel-method}
\alias{kernelMatrix,kernel-method}
\alias{kernelMult,kernel-method}
\alias{kernelFast,kernel-method}
\alias{kernelMatrix,rbfkernel-method}
\alias{kernelMatrix,polykernel-method}
\alias{kernelMatrix,vanillakernel-method}
\alias{kernelMatrix,tanhkernel-method}
\alias{kernelMatrix,laplacekernel-method}
\alias{kernelMatrix,anovakernel-method}
\alias{kernelMatrix,splinekernel-method}
\alias{kernelMatrix,besselkernel-method}
\alias{kernelMatrix,stringkernel-method}
\alias{kernelMult,rbfkernel,ANY-method}
\alias{kernelMult,splinekernel,ANY-method}
\alias{kernelMult,polykernel,ANY-method}
\alias{kernelMult,tanhkernel,ANY-method}
\alias{kernelMult,laplacekernel,ANY-method}
\alias{kernelMult,besselkernel,ANY-method}
\alias{kernelMult,anovakernel,ANY-method}
\alias{kernelMult,vanillakernel,ANY-method}
\alias{kernelMult,character,kernelMatrix-method}
\alias{kernelMult,stringkernel,ANY-method}
\alias{kernelPol,rbfkernel-method}
\alias{kernelPol,splinekernel-method}
\alias{kernelPol,polykernel-method}
\alias{kernelPol,tanhkernel-method}
\alias{kernelPol,vanillakernel-method}
\alias{kernelPol,anovakernel-method}
\alias{kernelPol,besselkernel-method}
\alias{kernelPol,laplacekernel-method}
\alias{kernelPol,stringkernel-method}
\alias{kernelFast,rbfkernel-method}
\alias{kernelFast,splinekernel-method}
\alias{kernelFast,polykernel-method}
\alias{kernelFast,tanhkernel-method}
\alias{kernelFast,vanillakernel-method}
\alias{kernelFast,anovakernel-method}
\alias{kernelFast,besselkernel-method}
\alias{kernelFast,laplacekernel-method}
\alias{kernelFast,stringkernel-method}
\alias{kernelFast,splinekernel-method}

\title{Kernel Matrix functions}
\description{
  \code{kernelMatrix} calculates the kernel matrix \eqn{K_{ij} = k(x_i,x_j)} or \eqn{K_{ij} =
    k(x_i,y_j)}.\cr
  \code{kernelPol} computes the quadratic kernel expression  \eqn{H = z_i z_j
    k(x_i,x_j)}, \eqn{H = z_i k_j k(x_i,y_j)}.\cr
    \code{kernelMult} calculates the kernel expansion \eqn{f(x_i) =
      \sum_{i=1}^m z_i  k(x_i,x_j)}\cr
    \code{kernelFast} computes the kernel matrix, identical
    to \code{kernelMatrix}, except that it also requires the squared
    norm of the first argument as additional input, useful in iterative
    kernel matrix calculations.
 }
\usage{
\S4method{kernelMatrix}{kernel}(kernel, x, y = NULL)

\S4method{kernelPol}{kernel}(kernel, x, y = NULL, z, k = NULL)

\S4method{kernelMult}{kernel}(kernel, x, y = NULL, z, blocksize = 256)

\S4method{kernelFast}{kernel}(kernel, x, y, a)
}

\arguments{
  \item{kernel}{the kernel function to be used to calculate the kernel
    matrix.
  This has to be a function of class \code{kernel}, i.e. which can be
  generated either one of the build in 
  kernel generating functions (e.g., \code{rbfdot} etc.) or a user defined
  function of class \code{kernel} taking two vector arguments and returning a scalar.}
  \item{x}{a data matrix to be used to calculate the kernel matrix, or a
  list of vector when a \code{stringkernel} is used}
  \item{y}{second data matrix to calculate the kernel matrix, or a
  list of vector when a \code{stringkernel} is used}
  \item{z}{a suitable vector or matrix}
  \item{k}{a suitable vector or matrix}
  \item{a}{the squared norm of \code{x}, e.g., \code{rowSums(x^2)}}
  \item{blocksize}{the kernel expansion computations are done block wise
    to avoid storing the kernel matrix into memory. \code{blocksize}
    defines the size of the computational blocks.}
}
\details{
  Common functions used during kernel based computations.\cr
  The \code{kernel} parameter can be set to any function, of class
  kernel, which computes the inner product in feature space between two
  vector arguments. \pkg{kernlab} provides the most popular kernel functions
  which can be initialized by using the following
  functions:
    \itemize{
      \item \code{rbfdot} Radial Basis kernel function
      \item \code{polydot} Polynomial kernel function
      \item \code{vanilladot} Linear kernel function
      \item \code{tanhdot} Hyperbolic tangent kernel function
      \item \code{laplacedot} Laplacian kernel function
      \item \code{besseldot} Bessel kernel function
      \item \code{anovadot} ANOVA RBF kernel function
      \item \code{splinedot} the Spline kernel 
    }  (see example.)

    \code{kernelFast} is mainly used in situations where columns of the
    kernel matrix are computed per invocation. In these cases,
    evaluating the norm of each row-entry over and over again would
    cause significant computational overhead.
  }
  
\value{
  \code{kernelMatrix} returns a symmetric diagonal semi-definite matrix.\cr
  \code{kernelPol} returns a matrix.\cr
  \code{kernelMult} usually returns a one-column matrix.
}
\author{Alexandros Karatzoglou \cr
\email{alexandros.karatzoglou@ci.tuwien.ac.at}}


\seealso{\code{\link{rbfdot}}, \code{\link{polydot}},
  \code{\link{tanhdot}}, \code{\link{vanilladot}}}

\examples{
## use the spam data
data(spam)
dt <- as.matrix(spam[c(10:20,3000:3010),-58])

## initialize kernel function 
rbf <- rbfdot(sigma = 0.05)
rbf

## calculate kernel matrix
kernelMatrix(rbf, dt)

yt <- as.matrix(as.integer(spam[c(10:20,3000:3010),58]))
yt[yt==2] <- -1

## calculate the quadratic kernel expression
kernelPol(rbf, dt, ,yt)

## calculate the kernel expansion
kernelMult(rbf, dt, ,yt)
}
\keyword{algebra}
\keyword{array}
