\name{gausspr}
\alias{gausspr}
\alias{gausspr,formula-method}
\alias{gausspr,vector-method}
\alias{gausspr,matrix-method}
\alias{show,gausspr-method}
\alias{predict,gausspr-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gaussian processes for regression and classification}
\description{
  \code{gausspr} is an implementation of Gaussian processes
  for classification and regression.
  
 }
\usage{

\S4method{gausspr}{formula}(x, data=NULL, ..., subset, na.action = na.omit)

\S4method{gausspr}{vector}(x,...)

\S4method{gausspr}{matrix}(x, y, type="classification", kernel="rbfdot", kpar=list(sigma = 0.1),
          var=1, tol=0.001, cross=0, fit=TRUE, ... , subset, na.action = na.omit)


}
%- maybe also 'usage' for other objects documented here.
\arguments{

\item{x}{a symbolic description of the model to be fit or a matrix or
  vector when a formula interface is not used. Note, that an when using
  the formula interface, that an intercept is always included, whether
  given in the formula or not.
  When not using a formula x is a matrix or vector
	  containg the variables in the model} 
	   
  \item{data}{an optional data frame containing the variables in the model.
          By default the variables are taken from the environment which
          `gausspr' is called from.}
	  
  \item{y}{a response vector with one label for each row/component of \code{x}. Can be either
    a factor (for classification tasks) or a numeric vector (for
    regression).}
  \item{type}{Type of problem. Either "classification" or "regression"}

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
    kpar parameter as well.}

  \item{var}{the initial noise variance}

  \item{tol}{tolerance of termination criterion (default: 0.001)}

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
\details{
 A Gaussian process is specified by a mean and a covariance function.
  The mean is a function of x (which is often the zero function), and
  the covariance
is a function \eqn{C(x,x')} which expresses the expected covariance between the
value of the function y at the points x and x'.
The actual function \eqn{y(x)} in any data modelling problem is assumed to be
a single sample from this Gaussian distribution.
Model parameter estimation in classification is done by a gradient descent algorithm.

}
\value{
An S4 object of class "gausspr" containing the fitted model along with
information.
 Accessor functions can be used to access the slots of the
  object which include :
  \item{alpha}{The resulting model parameters}
  \item{error}{Training error (if fit == TRUE)}
  
  }
  \references{
    Christopher K.I. Williams, Carl Edward Rasmussen\cr
    \emph{Gaussian Processes for Regression}\cr
    Advances in Neural Information Processing Systems, NIPS\cr
    \url{http://books.nips.cc/papers/files/nips08/0514.pdf}

  }
\author{Alexandros Karatzoglou \cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}}


\seealso{\code{\link{rvm}}, \code{\link{ksvm}} }

\examples{
# train model
data(iris)
test <- gausspr(Species~.,data=iris,var=2)
test
alpha(test)

# predict on the training set
predict(test,iris[,-5])


# create regression data
x <- seq(-20,20,0.1)
y <- sin(x)/x + rnorm(401,sd=0.03)

# regression with gaussian processes
foo <- gausspr(x, y)
foo

# predict and plot
ytest <- predict(foo, x)
plot(x, y, type ="l")
lines(x, ytest, col="red")
}
\keyword{classif}
\keyword{regression}
\keyword{nonlinear}
