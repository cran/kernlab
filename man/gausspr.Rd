\name{gausspr}
\alias{gausspr}
\alias{gausspr,formula-method}
\alias{gausspr,vector-method}
\alias{gausspr,matrix-method}
\alias{coef,gausspr-method}
\alias{show,gausspr-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Gaussian processes for regression and classification}
\description{
  \code{gausspr} is an implementation of Gaussian processes
  for classification and regression.
  
 }
\usage{


\S4method{gausspr}{formula}(x, data=NULL, ..., subset, na.action = na.omit, scaled = TRUE)

\S4method{gausspr}{vector}(x,...)

\S4method{gausspr}{matrix}(x, y, scaled = TRUE, type= NULL, kernel="rbfdot",
          kpar="automatic", var=1, variance.model = FALSE, tol=0.0005,
          cross=0, fit=TRUE, ... , subset, na.action = na.omit)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{a symbolic description of the model to be fit or a matrix or
  vector when a formula interface is not used. 
  When not using a formula x is a matrix or vector containing the variables in the model} 
	   
  \item{data}{an optional data frame containing the variables in the model.
          By default the variables are taken from the environment which
          `gausspr' is called from.}
	  
  \item{y}{a response vector with one label for each row/component of \code{x}. Can be either
    a factor (for classification tasks) or a numeric vector (for
    regression).}

  \item{type}{Type of problem. Either "classification" or "regression".
    Depending on whether \code{y} is a factor or not, the default
    setting for \code{type} is \code{classification} or \code{regression},
    respectively, but can be overwritten by setting an explicit value.\cr}

  \item{scaled}{A logical vector indicating the variables to be
    scaled. If \code{scaled} is of length 1, the value is recycled as
    many times as needed and all non-binary variables are scaled.
    Per default, data are scaled internally (both \code{x} and \code{y}
    variables) to zero mean and unit variance. The center and scale
    values are returned and used for later predictions.}

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

  \item{var}{the initial noise variance, (only for regression) (default
    : 0.001)}

  \item{variance.model}{build model for variance estimation (only for regression) (default
    : FALSE)}

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
    found. The default action is \code{na.omit}, which leads to
    rejection of cases with missing values on any required variable. An
    alternative is \code{na.fail}, which causes an error if \code{NA}
    cases are found. (NOTE: If given, this argument must be named.)}
	
  \item{\dots}{ additional parameters}

}
\details{
 A Gaussian process is specified by a mean and a covariance function.
  The mean is a function of \eqn{x} (which is often the zero function), and
  the covariance
is a function \eqn{C(x,x')} which expresses the expected covariance between the
value of the function \eqn{y} at the points \eqn{x} and \eqn{x'}.
The actual function \eqn{y(x)} in any data modeling problem is assumed to be
a single sample from this Gaussian distribution.
Laplace approximation is used for the parameter estimation in gaussian
processes for classification.\cr

The predict function can return class probabilities for 
classification problems by setting the \code{type} parameter to "probabilities".
For regression setting the \code{type} parameter to "variance" returns the estimated variance at each predicted point.
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
    C. K. I. Williams and D. Barber \cr
    Bayesian classification with Gaussian processes. \cr
    IEEE Transactions on Pattern Analysis and Machine Intelligence, 20(12):1342-1351, 1998\cr
    \url{http://www.dai.ed.ac.uk/homes/ckiw/postscript/pami_final.ps.gz}
  }
\author{Alexandros Karatzoglou \cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}}


\seealso{\code{\link{predict.gausspr}}, \code{\link{rvm}}, \code{\link{ksvm}}, \code{\link{gausspr-class}}, \code{\link{lssvm}} }

\examples{
# train model
data(iris)
test <- gausspr(Species~.,data=iris,var=2)
test
alpha(test)

# predict on the training set
predict(test,iris[,-5])
# class probabilities 
predict(test, iris[,-5], type="probabilities")

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


#predict and variance
x = c(-4, -3, -2, -1,  0, 0.5, 1, 2)
y = c(-2,  0,  -0.5,1,  2, 1, 0, -1)
plot(x,y)
foo2 <- gausspr(x, y, variance.model = TRUE)
xtest <- seq(-4,2,0.2)
lines(xtest, predict(foo2, xtest))
lines(xtest,
      predict(foo2, xtest)+2*sqrt(predict(foo2,xtest, type="variance")),
      col="red")
lines(xtest,
      predict(foo2, xtest)-2*sqrt(predict(foo2,xtest, type="variance")),
      col="red")

}
\keyword{classif}
\keyword{regression}
\keyword{nonlinear}
\keyword{methods}
