\name{inlearn}
\alias{inlearn}
\alias{inlearn,numeric-method}
\title{Onlearn object initialization}
\description{
 Online Kernel Algorithm object \code{onlearn} initialization function.
}
\usage{

\S4method{inlearn}{numeric}(d, kernel = "rbfdot", kpar = list(sigma = 0.1),
        type = "novelty", buffersize = 1000)
}
\arguments{
  \item{d}{the dimensionality of the data to be learned}
  
  \item{kernel}{the kernel function used in training and predicting.
    This parameter can be set to any function, of class kernel, which computes a dot product between two
    vector arguments. kernlab provides the most popular kernel functions
    which can be used by setting the kernel parameter to the following
    strings:
    \itemize{
      \item \code{rbfdot} Radial Basis kernel function "Gaussian"
      \item \code{polydot} Polynomial kernel function
      \item \code{vanilladot} Linear kernel function
      \item \code{tanhdot} Hyperbolic tangent kernel function
      \item \code{laplacedot} Laplacian kernel function
      \item \code{besseldot} Bessel kernel function
      \item \code{anovadot} ANOVA RBF kernel function
    }
    The kernel parameter can also be set to a user defined function of
    class kernel by passing the function name as an argument.
  }

  \item{kpar}{the list of hyper-parameters (kernel parameters).
    This is a list which contains the parameters to be used with the
    kernel function. For valid parameters for existing kernels are :
    \itemize{
      \item \code{sigma} inverse kernel width for the Radial Basis
      kernel function "rbfdot" and the Laplacian kernel "laplacedot".
      \item \code{degree, scale, offset} for the Polynomial kernel "polydot"
      \item \code{scale, offset} for the Hyperbolic tangent kernel
      function "tanhdot"
      \item \code{sigma, order, degree} for the Bessel kernel "besseldot".
      \item \code{sigma, degree} for the ANOVA kernel "anovadot".
    }
    Hyper-parameters for user defined kernels can be passed through the
    \code{kpar} parameter as well.}

   \item{type}{the type of problem to be learned by the online algorithm
     :
   \code{classification}, \code{regression}, \code{novelty}}
 \item{buffersize}{the size of the buffer to be used}
}
\details{
The \code{inlearn} is used to initialize a blank \code{onlearn} object.
}
\value{
 The function returns an \code{S4} object of class \code{onlearn} that
 can be used by the \code{onlearn} function.
}
\author{Alexandros Karatzoglou\cr
\email{alexandros.karatzoglou@ci.tuwien.ac.at}}

\seealso{ \code{\link{onlearn}}, \code{\link{onlearn-class}} }
\examples{

## create toy data set
x <- rbind(matrix(rnorm(100),,2),matrix(rnorm(100)+3,,2))
y <- matrix(c(rep(1,50),rep(-1,50)),,1)

## initialize onlearn object
on <- inlearn(2, kernel = "rbfdot", kpar = list(sigma = 0.2),
              type = "classification")

## learn one data point at the time
for(i in sample(1:100,100))
on <- onlearn(on,x[i,],y[i],nu=0.03,lambda=0.1)

sign(predict(on,x))

}
\keyword{classif}
\keyword{neural}
\keyword{regression}
\keyword{ts}
