\name{kpca}
\alias{kpca}
\alias{kpca,formula-method}
\alias{kpca,matrix-method}
\alias{kpca,kernelMatrix-method}
\alias{kpca,list-method}
\alias{predict,kpca-method}
\title{Kernel Principal Components Analysis}
\description{
Kernel Principal Components Analysis is a nonlinear form of principal
component analysis.}
\usage{
\S4method{kpca}{formula}(x, data = NULL, na.action, ...)

\S4method{kpca}{matrix}(x, kernel = "rbfdot", kpar = list(sigma = 0.1),
    features = 0, th = 1e-4, na.action = na.omit, ...)

\S4method{kpca}{kernelMatrix}(x, features = 0, th = 1e-4, ...)

\S4method{kpca}{list}(x, kernel = "stringdot", kpar = list(length = 4, lambda = 0.5),
    features = 0, th = 1e-4, na.action = na.omit, ...)
}

\arguments{
  \item{x}{the data matrix indexed by row or a formula describing the
    model, or a kernel Matrix of class \code{kernelMatrix}, or a list of character vectors} 
\item{data}{an optional data frame containing the variables in
	  the model (when using a formula).}
	
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
    kpar parameter as well.}
  
  \item{features}{Number of features (principal components) to
    return. (default: 0 , all)}

    \item{th}{the value of the eigenvalue under which principal
      components are ignored (only valid when features =  0). (default : 0.0001) }

  \item{na.action}{A function to specify the action to be taken if \code{NA}s are
          found. The default action is \code{na.omit}, which leads to rejection of cases
          with missing values on any required variable. An alternative
	  is \code{na.fail}, which causes an error if \code{NA} cases
	  are found. (NOTE: If given, this argument must be named.)}
    
  \item{\dots}{ additional parameters}

}
\details{Using kernel functions one can efficiently compute
  principal components in high-dimensional 
  feature spaces, related to input space by some non-linear map.\cr
  The data can be passed to the \code{kpca} function in a \code{matrix} or a
\code{data.frame}, in addition \code{kpca} also supports input in the form of a
kernel matrix of class \code{kernelMatrix} or as a list of character
vectors where a string kernel has to be used.
}
\value{
 An S4 object containing the principal component vectors along with the
 corresponding eigenvalues. 
  \item{pcv}{a matrix containing the principal component vectors (column
  wise)}
\item{eig}{The corresponding eigenvalues}
\item{rotated}{The original data projected (rotated) on the principal components}
\item{xmatrix}{The original data matrix}

all the slots of the object can be accessed by accessor functions.
}
\note{The predict function can be used to embed new data on the new space}
\references{
  Schoelkopf B., A. Smola, K.-R. Mueller :\cr
  \emph{Nonlinear component analysis as a kernel eigenvalue problem}\cr
  Neural Computation 10, 1299-1319\cr
  \url{http://mlg.anu.edu.au/~smola/papers/SchSmoMul98.pdf}
}
\author{Alexandros Karatzoglou \cr
\email{alexandros.karatzoglou@ci.tuwien.ac.at}}


\seealso{\code{\link{kcca}}, \code{pca}}
\examples{
# another example using the iris
data(iris)
test <- sample(1:150,20)

kpc <- kpca(~.,data=iris[-test,-5],kernel="rbfdot",
            kpar=list(sigma=0.2),features=2)

#print the principal component vectors
pcv(kpc)

#plot the data projection on the components
plot(rotated(kpc),col=as.integer(iris[-test,5]),
     xlab="1st Principal Component",ylab="2nd Principal Component")

#embed remaining points 
emb <- predict(kpc,iris[test,-5])
points(emb,col=as.integer(iris[test,5]))
}
\keyword{cluster}

