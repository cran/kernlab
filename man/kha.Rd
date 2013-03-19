\name{kha}
\alias{kha}
\alias{kha,formula-method}
\alias{kha,matrix-method}
\alias{predict,kha-method}
\encoding{latin1}
\title{Kernel Principal Components Analysis}
\description{
Kernel Hebbian Algorithm is a nonlinear iterative algorithm for principal
component analysis.}
\usage{
\S4method{kha}{formula}(x, data = NULL, na.action, ...)

\S4method{kha}{matrix}(x, kernel = "rbfdot", kpar = list(sigma = 0.1), features = 5, 
         eta = 0.005, th = 1e-4, maxiter = 10000, verbose = FALSE,
        na.action = na.omit, ...)
}

\arguments{
  \item{x}{ The data matrix indexed by row
    or a formula describing the model. Note, that an
          intercept is always included, whether given in the formula or
          not.} 
	\item{data}{an optional data frame containing the variables in
	  the model
	  (when using a formula).}
	
	\item{kernel}{the kernel function used in training and predicting.
    This parameter can be set to any function, of class kernel, which
    computes the inner product in feature space between two
    vector arguments (see \code{\link{kernels}}).
    \pkg{kernlab} provides the most popular kernel functions
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
    return. (default: 5)}
  \item{eta}{The hebbian learning rate (default : 0.005)}
    \item{th}{the smallest value of the convergence step (default : 0.0001) }
    \item{maxiter}{the maximum number of iterations.}
    \item{verbose}{print convergence every 100 iterations. (default : FALSE)}
  \item{na.action}{A function to specify the action to be taken if \code{NA}s are
          found. The default action is \code{na.omit}, which leads to rejection of cases
          with missing values on any required variable. An alternative
	  is \code{na.fail}, which causes an error if \code{NA} cases
	  are found. (NOTE: If given, this argument must be named.)}
    
  \item{\dots}{ additional parameters}
}



\details{The original form of KPCA can only be used on small data sets
  since it requires the estimation of the eigenvectors of a full kernel
  matrix. The Kernel Hebbian Algorithm iteratively estimates the Kernel
  Principal Components with only linear order memory complexity.
  (see ref. for more details)
}

\value{
 An S4 object containing the principal component vectors along with the
 corresponding normalization values. 
  \item{pcv}{a matrix containing the principal component vectors (column
  wise)}
\item{eig}{The normalization values}
\item{xmatrix}{The original data matrix}

all the slots of the object can be accessed by accessor functions.
}
\note{The predict function can be used to embed new data on the new space}
\references{Kwang In Kim, M.O. Franz and B. Schölkopf\cr
  \emph{Kernel Hebbian Algorithm for Iterative Kernel Principal Component Analysis}\cr
  Max-Planck-Institut für biologische Kybernetik, Tübingen (109)\cr
  \url{http://www.kyb.tuebingen.mpg.de/publications/pdfs/pdf2302.pdf}
}
\author{Alexandros Karatzoglou \cr
\email{alexandros.karatzoglou@ci.tuwien.ac.at}}


\seealso{\code{\link{kpca}}, \code{\link{kfa}}, \code{\link{kcca}}, \code{pca}}
\examples{
# another example using the iris
data(iris)
test <- sample(1:150,70)

kpc <- kha(~.,data=iris[-test,-5],kernel="rbfdot",
           kpar=list(sigma=0.2),features=2, eta=0.001, maxiter=65)

#print the principal component vectors
pcv(kpc)

#plot the data projection on the components
plot(predict(kpc,iris[,-5]),col=as.integer(iris[,5]),
     xlab="1st Principal Component",ylab="2nd Principal Component")

}
\keyword{cluster}

