\name{kfa}
\alias{kfa}
\alias{kfa,formula-method}
\alias{kfa,matrix-method}
\alias{show,kfa-method}
\alias{coef,kfa-method}
\title{Kernel Feature Analysis}
\description{
The Kernel Feature Analysis algorithm is an algorithm for extracting
structure from possibly high-dimensional data sets.
Similar to \code{kpca} a new basis for the data is found.
The data can then be projected on the new basis. 
}
\usage{
\S4method{kfa}{formula}(x, data = NULL, na.action = na.omit, ...)

\S4method{kfa}{matrix}(x, kernel = "rbfdot", kpar = list(sigma = 0.1),
   features = 0, subset = 59, normalize = TRUE, na.action = na.omit)
}

\arguments{

  \item{x}{ The data matrix indexed by row or a formula
    describing the model. Note, that an intercept is always
    included, whether given in the formula or not.}

  \item{data}{an optional data frame containing the variables in
          the model
          (when using a formula).}

  \item{kernel}{the kernel function used in training and predicting.
    This parameter can be set to any function, of class kernel, which
    computes an inner product in feature space between two
    vector arguments. \pkg{kernlab} provides the most popular kernel functions
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

  \item{subset}{the number of features sampled (used) from the data set}

 \item{normalize}{normalize the feature selected (default: TRUE)}  
  
  \item{na.action}{A function to specify the action to be taken if \code{NA}s are
          found. The default action is \code{na.omit}, which leads to rejection of cases
          with missing values on any required variable. An alternative
          is \code{na.fail}, which causes an error if \code{NA} cases
          are found. (NOTE: If given, this argument must be named.)}

  \item{\dots}{ additional parameters}
}
\details{
  Kernel Feature analysis is similar to Kernel PCA, but instead of
extracting eigenvectors of the training dataset in feature space, it
approximates the eigenvectors by selecting training patterns which are good
basis vectors for the training set. It works by choosing a fixed size
subset of the data set and scaling it to unit length (under the kernel).
It then chooses the features that maximize the value of the inner
product (kernel function) with the rest of the patterns.
}
\value{
  \code{kfa} returns an object of class \code{kfa} containing the
  features selected by the algorithm. 
  \item{xmatrix}{contains the features selected}
  \item{alpha}{contains the sparse alpha vector}

  The \code{predict} function can be used to embed new data points into to the
  selected feature base.
}
\references{Alex J. Smola, Olvi L. Mangasarian and Bernhard Schoelkopf\cr
  \emph{Sparse Kernel Feature Analysis}\cr
  Data Mining Institute Technical Report 99-04, October 1999\cr
  \url{ftp://ftp.cs.wisc.edu/pub/dmi/tech-reports/99-04.ps}
}

\author{Alexandros Karatzoglou\cr
\email{alexandros.karatzoglou@ci.tuwien.ac.at}}

\seealso{\code{\link{kpca}}, \code{\link{kfa-class}}}
\examples{
data(promotergene)
f <- kfa(~.,data=promotergene,features=2,kernel="rbfdot",
         kpar=list(sigma=0.01))
plot(predict(f,promotergene),col=as.numeric(promotergene[,1]))
}
\keyword{cluster}

