\name{dots}
\alias{dots}
\alias{kernels}
\alias{rbfdot}
\alias{polydot}
\alias{tanhdot}
\alias{vanilladot}
\alias{laplacedot}
\alias{besseldot}
\alias{anovadot}
\alias{fourierdot}
\alias{splinedot}
\alias{kpar}
\alias{kfunction}
\alias{show,kernel-method}
\title{Kernel Functions}
\description{
  The kernel generating functions provided in kernlab. \cr
  The Gaussian RBF kernel \eqn{k(x,x') = \exp(-\sigma \|x - x'\|^2)} \cr
  The Polynomial kernel \eqn{k(x,x') = (scale <x, x'> + offset)^{degree}}\cr
  The Linear kernel \eqn{k(x,x') = <x, x'>}\cr
  The Hyperbolic tangent kernel \eqn{k(x, x') = \tanh(scale <x, x'> +  offset)}\cr
  The Laplacian kernel \eqn{k(x,x') = \exp(-\sigma \|x - x'\|)} \cr
  The Bessel kernel \eqn{k(x,x') = (- Bessel_{(\nu+1)}^n \sigma \|x - x'\|^2)} \cr
  The ANOVA RBF kernel \eqn{k(x,x') = \sum_{1\leq i_1 \ldots < i_D \leq
      N} \prod_{d=1}^D k(x_{id}, {x'}_{id})} where k(x,x) is a Gaussian
  RBF kernel. \cr
  The Spline kernel \eqn{ \prod_{d=1}^D 1 + x_i x_j + x_i x_j min(x_i,
    x_j)  - \frac{x_i + x_j}{2} min(x_i,x_j)^2 +
    \frac{min(x_i,x_j)^3}{3}} \\
  The String kernels (see \code{stringdot}.
}
\usage{
rbfdot(sigma = 1)

polydot(degree = 1, scale = 1, offset = 1)

tanhdot(scale = 1, offset = 1)

vanilladot()

laplacedot(sigma = 1)

besseldot(sigma = 1, order = 1, degree = 1)

anovadot(sigma = 1, degree = 1)

splinedot()
}

\arguments{
  \item{sigma}{The inverse kernel width used by the Gaussian the
    Laplacian, the Bessel and the ANOVA kernel }
  \item{degree}{The degree of the polynomial, bessel or ANOVA
    kernel function. This has to be an positive integer.}
  \item{scale}{The scaling parameter of the polynomial and tangent
    kernel is a convenient way of normalizing
    patterns without the need to modify the data itself}
  \item{offset}{The offset used in a polynomial or hyperbolic tangent
    kernel}
  \item{order}{The order of the Bessel function to be used as a kernel}
}
\details{
  The kernel generating functions are used to initialize a kernel
  function
  which calculates the dot (inner) product between two feature vectors in a
  Hilbert Space. These functions can be passed as a \code{kernel} argument on almost all
  functions in \pkg{kernlab}(e.g., \code{ksvm}, \code{kpca}  etc).
  
  Although using one of the existing kernel functions as a
  \code{kernel} argument in various functions in \pkg{kernlab} has the
  advantage that optimized code is used to calculate various kernel expressions,
  any other function implementing a dot product of class \code{kernel} can also be used as a kernel
  argument. This allows the user to use, test and develop special kernels
  for a given data set or algorithm.
  For details on the string kernels see \code{stringdot}.
  }
\value{
 Return an S4 object of class \code{kernel} which extents the
 \code{function} class. The resulting function implements the given
 kernel calculating the inner (dot) product between two vectors.
 \item{kpar}{a list containing the kernel parameters (hyperparameters)
   used.}
 The kernel parameters can be accessed by the \code{kpar} function.
 }

\author{Alexandros Karatzoglou\cr
  \email{alexandros.karatzoglou@ci.tuwien.ac.at}}

\note{If the offset in the Polynomial kernel is set to $0$, we obtain homogeneous polynomial
  kernels, for positive values, we have inhomogeneous
  kernels. Note that for negative values the kernel does not satisfy Mercer's
  condition and thus the optimizers may fail. \cr

  In the Hyperbolic tangent kernel if the offset is negative the likelihood of obtaining a kernel
  matrix that is not positive definite is much higher (since then even some
  diagonal elements may be negative), hence if this kernel has to be used, the
  offset should always be positive. Note, however, that this is no guarantee
  that the kernel will be positive.
}




\seealso{\code{stringdot}, \code{\link{kernelMatrix} }, \code{\link{kernelMult}}, \code{\link{kernelPol}}}
\examples{
rbfkernel <- rbfdot(sigma = 0.1)
rbfkernel

kpar(rbfkernel)

## create two vectors
x <- rnorm(10)
y <- rnorm(10)

## calculate dot product
rbfkernel(x,y)

}
\keyword{symbolmath}

