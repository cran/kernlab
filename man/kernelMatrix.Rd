\name{kernelMatrix}
\alias{kernelMatrix}
\alias{kernelMult}
\alias{kernelPol}
\alias{kernelPol,kernel-method}
\alias{kernelMatrix,kernel-method}
\alias{kernelMult,kernel-method}
\alias{kernelMatrix,rbfkernel,matrix-method}
\alias{kernelMatrix,polykernel,matrix-method}
\alias{kernelMatrix,vanillakernel,matrix-method}
\alias{kernelMatrix,tanhkernel,matrix-method}
\alias{kernelMult,rbfkernel,matrix-method}
\alias{kernelMult,polykernel,matrix-method}
\alias{kernelMult,tanhkernel,matrix-method}
\alias{kernelMult,vanillakernel,matrix-method}
\alias{kernelPol,rbfkernel,matrix-method}
\alias{kernelPol,polykernel,matrix-method}
\alias{kernelPol,tanhkernel,matrix-method}
\alias{kernelPol,vanillakernel,matrix-method}

\title{Kernel Matrix functions}
\description{
  \code{kernelMatrix} calculates the kernel matrix \eqn{K_{ij} = k(x_i,x_j)} or \eqn{K_{ij} =
    k(x_i,y_j)}.\cr
  \code{kernelPol} computes the quadratic kernel expression  \eqn{H = z_i z_j
    k(x_i,x_j)}, \eqn{H = z_i k_j k(x_i,y_j)}.\cr
    \code{kernelMult} calculates the kernel expansion \eqn{f(x_i) =
      sum_{i=1}^m z k(x_i,x_j)} 
 }
\usage{
\S4method{kernelMatrix}{kernel}(kernel, x, y = NULL)

\S4method{kernelPol}{kernel}(kernel, x, y = NULL, z, k = NULL)

\S4method{kernelMult}{kernel}(kernel, x, y = NULL, z, blocksize = 256)
}

\arguments{
  \item{kernel}{the kernel function to be used to calculate the kernel
    matrix.
  This has to be a function of class \code{kernel}}
  \item{x}{a data matrix to be used to calculate the kernel matrix}
  \item{y}{second data matrix to calculate the kernel matrix}
  \item{z}{a suitable vector or matrix}
  \item{k}{a suitable vector or matrix}
  \item{blocksize}{the kernel expansion computations are done block wise
    to avoid storing the kernel matrix into memory. \code{blocksize}
    defines the size of the computational blocks.}
}
\details{
  Common functions used during kernel based computations.\cr
  This \code{kernel} parameter can be set to any function, of class
  kernel, which computes a dot product between two
  vector arguments. kernlab provides the most popular kernel functions
  which can be initialized by using the following
  functions:
    \itemize{
      \item \code{rbfdot} (Radial Basis kernel function)
      \item \code{polydot} (Polynomial kernel function)
      \item \code{vanilladot} (Linear kernel function)
      \item \code{tanhdot} (Hyperbolic tangent kernel function)


    }  (see example.)
  }
  
\value{
  \code{kernelMatrix} returns a symmetric diagonal semi-definite matrix.
  \code{kernelPol} returns a matrix.
  \code{kernelMult}{usually returns a one-column matrix}
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
