\name{inchol}
\alias{inchol}
\alias{inchol,matrix-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Incomplete Cholesky decomposition}
\description{
  \code{inchol} computes the incomplete Cholesky decomposition
  of the kernel matrix from a data matrix. 
}
\usage{
inchol(x, kernel="rbfdot", kpar=list(sigma=0.1), tol = 0.001, 
            maxiter = dim(x)[1], blocksize = 50, verbose = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The data matrix indexed by row}
  \item{kernel}{the kernel function used in training and predicting.
    This parameter can be set to any function, of class \code{kernel},
    which computes the inner product in feature space between two
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
       \item \code{splinedot} Spline kernel 

    }
    The kernel parameter can also be set to a user defined function of
    class kernel by passing the function name as an argument.
  }

  \item{kpar}{the list of hyper-parameters (kernel parameters).
    This is a list which contains the parameters to be used with the
    kernel function. Valid parameters for existing kernels are :
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
    kpar parameter as well.
  }
  
  \item{tol}{algorithm stops when remaining pivots bring less accuracy
    then \code{tol} (default: 0.001)}
  \item{maxiter}{maximum number of iterations and columns in \eqn{Z}}
  \item{blocksize}{add this many columns to matrix per iteration}
  \item{verbose}{print info on algorithm convergence}
}
\details{An incomplete cholesky decomposition calculates
  \eqn{Z} where \eqn{K= ZZ'} \eqn{K} being the kernel matrix.
  Since the rank of a kernel matrix is usually low, \eqn{Z} tends to be smaller
  then the complete kernel matrix. The decomposed matrix can be
  used to create memory efficient kernel-based algorithms without the
  need to compute and store a complete kernel matrix in memory.}
\value{
  An S4 object of class "inchol" which is an extension of the class
  "matrix". The object is the decomposed kernel matrix along with 
  the slots :
  \item{pivots}{Indices on which pivots where done}
  \item{diagresidues}{Residuals left on the diagonal}
  \item{maxresiduals}{Residuals picked for pivoting}

  slots can be accessed either by \code{object@slot}
or by accessor functions with the same name (e.g., \code{pivots(object))}}

\references{
      Francis R. Bach, Michael I. Jordan\cr
      \emph{Kernel Independent Component Analysis}\cr
      Journal of Machine Learning Research  3, 1-48\cr
      \url{http://www.jmlr.org/papers/volume3/bach02a/bach02a.pdf}
    }
    
\author{Alexandros Karatzoglou (based on Matlab code by 
  S.V.N. (Vishy) Vishwanathan and Alex Smola)\cr
\email{alexandros.karatzoglou@ci.tuwien.ac.at}}

\seealso{\code{\link{csi}}, \code{\link{inchol-class}}, \code{\link{chol}}}
\examples{

data(iris)
datamatrix <- as.matrix(iris[,-5])
# initialize kernel function
rbf <- rbfdot(sigma=0.1)
rbf
Z <- inchol(datamatrix,kernel=rbf)
dim(Z)
pivots(Z)
# calculate kernel matrix
K <- crossprod(t(Z))
# difference between approximated and real kernel matrix
(K - kernelMatrix(kernel=rbf, datamatrix))[6,]

}
\keyword{methods}
\keyword{algebra}
\keyword{array}
