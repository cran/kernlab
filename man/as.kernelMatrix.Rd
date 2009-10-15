\name{as.kernelMatrix}
\docType{methods}
\alias{kernelMatrix-class}
\alias{as.kernelMatrix}
\alias{as.kernelMatrix-methods}
\alias{as.kernelMatrix,matrix-method}
\title{Assing kernelMatrix class to matrix objects}

\description{\code{as.kernelMatrix} in package \pkg{kernlab} can be used 
  to coerce the kernelMatrix class to matrix objects representing a
  kernel matrix.  These matrices can then be used with the kernelMatrix
  interfaces which most of the functions in \pkg{kernlab} support.}

\usage{
\S4method{as.kernelMatrix}{matrix}(x, center = FALSE)
}
\arguments{
  \item{x}{matrix to be assigned the \code{kernelMatrix} class }
  \item{center}{center the kernel matrix in feature space (default: FALSE) }
}

\author{
  Alexandros Karatzoglou \cr
    \email{alexandros.karatzoglou@ci.tuwien.ac.at}
    }

\seealso{\code{\link{kernelMatrix}}, \code{\link{dots}}}
    
\keyword{methods}


\examples{
## Create toy data
x <- rbind(matrix(rnorm(10),,2),matrix(rnorm(10,mean=3),,2))
y <- matrix(c(rep(1,5),rep(-1,5)))

### Use as.kernelMatrix to label the cov. matrix as a kernel matrix
### which is eq. to using a linear kernel 

K <- as.kernelMatrix(crossprod(t(x)))

K

svp2 <- ksvm(K, y, type="C-svc")

svp2

}
