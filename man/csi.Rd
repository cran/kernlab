\name{csi}
\docType{methods}
\alias{csi}
\alias{csi-methods}
\alias{csi,matrix-method}
\title{Cholesky decomposition with Side Information}
\description{
 The \code{csi} function in \pkg{kernlab} is an implementation of an
 incomplete Cholesky decomposition algorithm which exploits side
 information (e.g., classification labels, regression responses) to
 compute a low rank decomposition of a kernel matrix from the data.
}
\usage{
\S4method{csi}{matrix}(x, y, kernel="rbfdot", kpar=list(sigma=0.1), rank,
centering = TRUE, kappa = 0.99 ,delta = 40 ,tol = 1e-5)
}

\arguments{
  \item{x}{The data matrix indexed by row}

  \item{y}{the classification labels or regression responses. In
    classification y is a \eqn{m \times n} matrix where \eqn{m}
    the number of data and \eqn{n} the number of classes \eqn{y} and \eqn{y_i} is 1 if
    the corresponding x belongs to class i.}
  
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
      \item \code{stringdot} String kernel
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
  \item{rank}{maximal rank of the computed kernel matrix}

  \item{centering}{if \code{TRUE} centering is performed (default: TRUE)}

  \item{kappa}{trade-off between approximation of K and prediction of Y (default: 0.99)}
  
  \item{delta}{number of columns of cholesky performed in advance (default: 40)}

  \item{tol}{minimum gain at each iteration (default: 1e-4)}
}


\details{An incomplete cholesky decomposition calculates
  \eqn{Z} where \eqn{K= ZZ'} \eqn{K} being the kernel matrix.
  Since the rank of a kernel matrix is usually low, \eqn{Z} tends to
  be smaller then the complete kernel matrix. The decomposed matrix can be
  used to create memory efficient kernel-based algorithms without the
  need to compute and store a complete kernel matrix in memory. \cr
  \code{csi} uses the class labels, or regression responses to compute a
  more appropriate approximation for the problem at hand considering the
  additional information from the response variable. }

\value{
  An S4 object of class "csi" which is an extension of the class
  "matrix". The object is the decomposed kernel matrix along with 
  the slots :
  \item{pivots}{Indices on which pivots where done}
  \item{diagresidues}{Residuals left on the diagonal}
  \item{maxresiduals}{Residuals picked for pivoting}
  \item{predgain}{predicted gain before adding each column}
  \item{truegain}{actual gain after adding each column}
  \item{Q}{QR decomposition of the kernel matrix}
  \item{R}{QR decomposition of the kernel matrix}
  
  slots can be accessed either by \code{object@slot}
  or by accessor functions with the same name
  (e.g., \code{pivots(object))}}

\references{ 
      Francis R. Bach, Michael I. Jordan\cr
      \emph{Predictive low-rank decomposition for kernel methods.}\cr
      Proceedings of the Twenty-second International Conference on Machine Learning (ICML) 2005\cr
      \url{http://cmm.ensmp.fr/~bach/bach_jordan_csi.pdf}
    }

\author{Alexandros Karatzoglou (based on Matlab code by 
  Francis Bach)\cr
\email{alexandros.karatzoglou@ci.tuwien.ac.at}}

\seealso{\code{\link{inchol}}, \code{\link{chol}}, \code{\link{csi-class}}}
\examples{

data(iris)

## create multidimensional y matrix
yind <- t(matrix(1:3,3,150))
ymat <- matrix(0, 150, 3)
ymat[yind==as.integer(iris[,5])] <- 1

datamatrix <- as.matrix(iris[,-5])
# initialize kernel function
rbf <- rbfdot(sigma=0.1)
rbf
Z <- csi(datamatrix,ymat, kernel=rbf, rank = 30)
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

