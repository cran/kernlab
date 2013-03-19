\name{lssvm}
\docType{methods}
\alias{lssvm}
\alias{lssvm-methods}
\alias{lssvm,formula-method}
\alias{lssvm,vector-method}
\alias{lssvm,matrix-method}
\alias{lssvm,list-method}
\alias{lssvm,kernelMatrix-method}
\alias{show,lssvm-method}
\alias{coef,lssvm-method}
\alias{predict,lssvm-method}
\title{Least Squares Support Vector Machine}
\description{
  The \code{lssvm} function is an
  implementation of the Least Squares SVM. \code{lssvm} includes a
  reduced version of Least Squares SVM using a decomposition of the
  kernel matrix which is calculated by the \code{csi} function.
}

\usage{

\S4method{lssvm}{formula}(x, data=NULL, ..., subset, na.action = na.omit, scaled = TRUE)

\S4method{lssvm}{vector}(x, ...)

\S4method{lssvm}{matrix}(x, y, scaled = TRUE, kernel = "rbfdot", kpar = "automatic",
      type = NULL, tau = 0.01, reduced = TRUE, tol = 0.0001,
      rank = floor(dim(x)[1]/3), delta = 40, cross = 0, fit = TRUE,
      ..., subset, na.action = na.omit)

\S4method{lssvm}{kernelMatrix}(x, y, type = NULL, tau = 0.01,
      tol = 0.0001, rank = floor(dim(x)[1]/3), delta = 40, cross = 0,
      fit = TRUE, ...)

\S4method{lssvm}{list}(x, y, scaled = TRUE,
      kernel = "stringdot", kpar = list(length=4, lambda = 0.5),
      type = NULL, tau = 0.01, reduced = TRUE, tol = 0.0001,
      rank = floor(dim(x)[1]/3), delta = 40, cross = 0, fit = TRUE,
      ..., subset)
}

\arguments{

\item{x}{a symbolic description of the model to be fit, a matrix or
  vector containing the training data when a formula interface is not
  used or a \code{kernelMatrix} or a list of character vectors.}
	   
  \item{data}{an optional data frame containing the variables in the model.
          By default the variables are taken from the environment which
          `lssvm' is called from.}
	  
  \item{y}{a response vector with one label for each row/component of \code{x}. Can be either
    a factor (for classification tasks) or a numeric vector (for
    classification or regression - currently nor supported -).}

   \item{scaled}{A logical vector indicating the variables to be
    scaled. If \code{scaled} is of length 1, the value is recycled as
    many times as needed and all non-binary variables are scaled.
    Per default, data are scaled internally to zero mean and unit
    variance. The center and scale values are returned and used for later predictions.}
  
  \item{type}{Type of problem. Either "classification" or "regression".
    Depending on whether \code{y} is a factor or not, the default
    setting for \code{type} is "classification" or "regression" respectively,
    but can be overwritten by setting an explicit value. (regression is
currently not supported)\cr}

 \item{kernel}{the kernel function used in training and predicting.
    This parameter can be set to any function, of class kernel, which computes a dot product between two
    vector arguments. kernlab provides the most popular kernel functions
    which can be used by setting the kernel parameter to the following
    strings:
    \itemize{

      \item \code{rbfdot} Radial Basis kernel "Gaussian"

      \item \code{polydot} Polynomial kernel 

      \item \code{vanilladot} Linear kernel 

      \item \code{tanhdot} Hyperbolic tangent kernel 

      \item \code{laplacedot} Laplacian kernel 

      \item \code{besseldot} Bessel kernel 

      \item \code{anovadot} ANOVA RBF kernel 

      \item \code{splinedot} Spline kernel

      \item \code{stringdot} String kernel
    }
    Setting the kernel parameter to "matrix" treats \code{x} as a kernel
    matrix calling the \code{kernelMatrix} interface.\cr
    
    
    The kernel parameter can also be set to a user defined function of
    class kernel by passing the function name as an argument.
  }

  \item{kpar}{

    the list of hyper-parameters (kernel parameters).
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

      \item \code{length, lambda, normalized} for the "stringdot" kernel
      where length is the length of the strings considered, lambda the
      decay factor and normalized a logical parameter determining if the
      kernel evaluations should be normalized.
    }
   Hyper-parameters for user defined kernels can be passed through the
   kpar parameter as well.\cr

    \code{kpar} can also be set to the string "automatic" which uses the heuristics in 
    \code{\link{sigest}} to calculate a good \code{sigma} value for the
    Gaussian RBF or Laplace kernel, from the data. (default = "automatic").
  }

  \item{tau}{the regularization parameter (default 0.01) }

  
  \item{reduced}{if set to \code{FALSE} the full linear problem of the
    lssvm is solved, when \code{TRUE} a reduced method using \code{csi} is used.}
  
  \item{rank}{the maximal rank of the decomposed kernel matrix, see
    \code{csi}}

  \item{delta}{number of columns of cholesky performed in advance, see
    \code{csi} (default 40)}
 
  \item{tol}{tolerance of termination criterion for the \code{csi}
    function, lower tolerance leads to more precise approximation but
    may increase the training time and the decomposed matrix size (default: 0.0001)}

  \item{fit}{indicates whether the fitted values should be computed and
    included in the model or not (default: 'TRUE')}

  \item{cross}{if a integer value k>0 is specified, a k-fold cross
    validation on the training data is performed to assess the
    quality of the model: the Mean Squared Error for regression}

  \item{subset}{An index vector specifying the cases to be used in the
    training sample.  (NOTE: If given, this argument must be
    named.)}

  \item{na.action}{A function to specify the action to be taken if \code{NA}s are
    found. The default action is \code{na.omit}, which leads to rejection of cases
    with missing values on any required variable. An alternative
    is \code{na.fail}, which causes an error if \code{NA} cases
    are found. (NOTE: If given, this argument must be named.)}
	
  \item{\dots}{ additional parameters}

}
\details{Least Squares Support Vector Machines are reformulation to the
  standard SVMs that lead to solving linear KKT systems.
  The algorithm is based on the minimization of a classical penalized
  least-squares cost function. The current implementation approximates
  the kernel matrix by an incomplete Cholesky factorization obtained by
  the \code{\link{csi}} function, thus the solution is an approximation
  to the exact solution of the lssvm optimization problem. The quality
  of the solution depends on the approximation and can be influenced by
  the "rank" , "delta", and "tol" parameters.
}


\value{
    An S4 object of class \code{"lssvm"} containing the fitted model,
  Accessor functions can be used to access the slots of the object (see
  examples) which include:
 \item{alpha}{the parameters of the \code{"lssvm"}}
\item{coef}{the model coefficients (identical to alpha)}
\item{b}{the model offset.}
\item{xmatrix}{the training data used by the model}
}

  \references{
     J. A. K. Suykens and  J. Vandewalle\cr
     \emph{Least Squares Support Vector Machine Classifiers}\cr
     Neural Processing Letters vol. 9, issue 3, June 1999\cr
     }
\author{Alexandros Karatzoglou \cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}}


\seealso{\code{\link{ksvm}}, \code{\link{gausspr}}, \code{\link{csi}} }

\examples{
## simple example
data(iris)

lir <- lssvm(Species~.,data=iris)

lir

lirr <- lssvm(Species~.,data= iris, reduced = FALSE)

lirr

## Using the kernelMatrix interface

iris <- unique(iris)

rbf <- rbfdot(0.5)

k <- kernelMatrix(rbf, as.matrix(iris[,-5]))

klir <- lssvm(k, iris[, 5])

klir

pre <- predict(klir, k)
}
\keyword{classif}
\keyword{nonlinear}
\keyword{methods}

