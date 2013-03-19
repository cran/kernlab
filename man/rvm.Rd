\name{rvm}
\alias{rvm}
\alias{rvm-methods}
\alias{rvm,formula-method}
\alias{rvm,list-method}
\alias{rvm,vector-method}
\alias{rvm,kernelMatrix-method}
\alias{rvm,matrix-method}
\alias{show,rvm-method}
\alias{predict,rvm-method}
\alias{coef,rvm-method}
\title{Relevance Vector Machine}

\description{
  The Relevance Vector Machine is a Bayesian model for regression and
  classification of identical functional form to the support vector
  machine.
  The \code{rvm} function currently supports only regression.
}
\usage{
\S4method{rvm}{formula}(x, data=NULL, ..., subset, na.action = na.omit)

\S4method{rvm}{vector}(x, ...)

\S4method{rvm}{matrix}(x, y, type="regression",
    kernel="rbfdot", kpar="automatic",
    alpha= ncol(as.matrix(x)), var=0.1, var.fix=FALSE, iterations=100,
    verbosity = 0, tol = .Machine$double.eps, minmaxdiff = 1e-3,
    cross = 0, fit = TRUE, ... , subset, na.action = na.omit) 

\S4method{rvm}{list}(x, y, type = "regression",
    kernel = "stringdot", kpar = list(length = 4, lambda = 0.5),
    alpha = 5, var = 0.1, var.fix = FALSE, iterations = 100,
    verbosity = 0, tol = .Machine$double.eps, minmaxdiff = 1e-3,
    cross = 0, fit = TRUE, ..., subset, na.action = na.omit)
}



\arguments{

  \item{x}{a symbolic description of the model to be fit.
    When not using a formula x can be a matrix or vector containing the training
    data or a kernel matrix of class \code{kernelMatrix} of the training data
    or a list of character vectors (for use with the string
    kernel). Note, that the intercept is always excluded, whether
    given in the formula or not.}
  \item{data}{an optional data frame containing the variables in the model.
          By default the variables are taken from the environment which
          `rvm' is called from.}
	  
  \item{y}{a response vector with one label for each row/component of \code{x}. Can be either
    a factor (for classification tasks) or a numeric vector (for
    regression).}
  \item{type}{\code{rvm} can only be used for regression at the moment.}


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
      
      \item \code{length, lambda, normalized} for the "stringdot" kernel
      where length is the length of the strings considered, lambda the
      decay factor and normalized a logical parameter determining if the
      kernel evaluations should be normalized.
    }

    Hyper-parameters for user defined kernels can be passed through the
    kpar parameter as well. In the case of a Radial Basis kernel function (Gaussian)
     kpar can also be set to the string "automatic" which uses the heuristics in 
    \code{\link{sigest}} to calculate a good \code{sigma} value for the
    Gaussian RBF or Laplace kernel, from the data.
    (default = "automatic").}

  \item{alpha}{The initial alpha vector. Can be either a vector of
    length equal to the number of data points or a single number.}

  \item{var}{the initial noise variance}

  \item{var.fix}{Keep noise variance fix during iterations (default: FALSE)}

  \item{iterations}{Number of iterations allowed (default: 100)} 
  \item{tol}{tolerance of termination criterion}
  \item{minmaxdiff}{termination criteria. Stop when max difference is
    equal to this parameter (default:1e-3) }
  \item{verbosity}{print information on algorithm convergence (default
    = FALSE)}
  \item{fit}{indicates whether the fitted values should be computed and
          included in the model or not (default: TRUE)}
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
\details{The Relevance Vector Machine typically leads to sparser models
  then the SVM. It also performs better in many cases (specially in
  regression). 
  }
\value{
  An S4 object of class "rvm" containing the fitted model.
  Accessor functions can be used to access the slots of the
  object which include :
  
  
  \item{alpha}{The resulting relevance vectors}
  \item{alphaindex}{ The index of the resulting relevance vectors in the data
    matrix}
  \item{nRV}{Number of relevance vectors}
  \item{RVindex}{The indexes of the relevance vectors}
  \item{error}{Training error (if \code{fit = TRUE})}
  
  ...
}
\references{
      Tipping, M. E.\cr
      \emph{Sparse Bayesian learning and the relevance vector machine}\cr
      Journal of Machine Learning Research  1, 211-244\cr
      \url{http://www.jmlr.org/papers/volume1/tipping01a/tipping01a.pdf}
    }
\author{ Alexandros Karatzoglou \cr
\email{alexandros.karatzoglou@ci.tuwien.ac.at}}

\seealso{ \code{\link{ksvm}}}
\examples{
# create data
x <- seq(-20,20,0.1)
y <- sin(x)/x + rnorm(401,sd=0.05)

# train relevance vector machine
foo <- rvm(x, y)
foo
# print relevance vectors
alpha(foo)
RVindex(foo)

# predict and plot
ytest <- predict(foo, x)
plot(x, y, type ="l")
lines(x, ytest, col="red")
}
\keyword{regression}
\keyword{nonlinear}
