\name{kqr}
\alias{kqr}
\alias{kqr,formula-method}
\alias{kqr,vector-method}
\alias{kqr,matrix-method}
\alias{kqr,list-method}
\alias{kqr,kernelMatrix-method}
\alias{coef,kqr-method}
\alias{show,kqr-method}


\title{Kernel Quantile Regression.}
\description{The Kernel Quantile Regression algorithm \code{kqr} performs
  non-parametric Quantile Regression.}
\usage{
\S4method{kqr}{formula}(x, data=NULL, ..., subset, na.action = na.omit, scaled = TRUE)

\S4method{kqr}{vector}(x,...)

\S4method{kqr}{matrix}(x, y, scaled = TRUE, tau = 0.5, C = 0.1, kernel = "rbfdot",
    kpar = "automatic", reduced = FALSE, rank = dim(x)[1]/6,
    fit = TRUE, cross = 0, na.action = na.omit)

\S4method{kqr}{kernelMatrix}(x, y, tau = 0.5, C = 0.1, fit = TRUE, cross = 0)

\S4method{kqr}{list}(x, y, tau = 0.5, C = 0.1, kernel = "strigdot",
    kpar= list(length=4, C=0.5), fit = TRUE, cross = 0)
}

\arguments{
  \item{x}{e data or a symbolic description of the model to be fit.
    When not using a formula x can be a matrix or vector containing
    the training data or a kernel matrix of class \code{kernelMatrix}
    of the training data or a list of character vectors (for use
    with the string kernel). Note, that the intercept is always
    excluded, whether given in the formula or not.}
  
  \item{data}{an optional data frame containing the variables in the model.
    By default the variables are taken from the environment which
    \code{kqr} is called from.}
	  
  \item{y}{a numeric vector or a column matrix containing the response.}

  \item{scaled}{A logical vector indicating the variables to be
    scaled. If \code{scaled} is of length 1, the value is recycled as
    many times as needed and all non-binary variables are scaled.
    Per default, data are scaled internally (both \code{x} and \code{y}
    variables) to zero mean and unit variance. The center and scale
    values are returned and used for later predictions. (default: TRUE)}

  \item{tau}{the quantile to be estimated, this is generally a number
    strictly between 0 and 1. For 0.5 the median is calculated.
    (default: 0.5)}

  \item{C}{the cost regularization parameter. This parameter controls
    the smoothness of the fitted function, essentially higher
    values for C lead to less smooth functions.(default: 1)}
	
 \item{kernel}{the kernel function used in training and predicting.
   This parameter can be set to any function, of class kernel, which computes a dot product between two
   vector arguments. \code{kernlab} provides the most popular kernel functions
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
     \item \code{lenght, lambda, normalized} for the "stringdot" kernel
     where length is the length of the strings considered, lambda the
     decay factor and normalized a logical parameter determining if the
     kernel evaluations should be normalized.
   }
   
   Hyper-parameters for user defined kernels can be passed
   through the \code{kpar} parameter as well. In the case of a Radial
   Basis kernel function (Gaussian) kpar can also be set to the
   string "automatic" which uses the heuristics in  'sigest' to
   calculate a good 'sigma' value for the Gaussian RBF or
   Laplace kernel, from the data. (default = "automatic").
 }
 
 \item{reduced}{use an incomplete cholesky decomposition to calculate a
   decomposed form \eqn{Z} of the kernel Matrix \eqn{K} (where \eqn{K = ZZ'}) and
   perform the calculations with \eqn{Z}. This might be useful when
   using \code{kqr} with large datasets since normally an n times n
   kernel matrix would be computed. Setting \code{reduced} to \code{TRUE}
   makes use of \code{csi} to compute a decomposed form instead and
   thus only a \eqn{n \times m} matrix where \eqn{m < n} and \eqn{n} the sample size is
   stored in memory (default: FALSE)}
 
 \item{rank}{the rank m of the decomposed matrix calculated when using an
   incomplete cholesky decomposition. This parameter is only
   taken into account when \code{reduced} is \code{TRUE}(default :
   dim(x)[1]/6)}
	
 \item{fit}{indicates whether the fitted values should be computed and
   included in the model or not (default: 'TRUE')}
 
 \item{cross}{if a integer value k>0 is specified, a k-fold cross
   validation on the training data is performed to assess the
   quality of the model: the Pinball loss and the for quantile regression}
 
 \item{subset}{An index vector specifying the cases to be used in the
   training sample.  (NOTE: If given, this argument must be
   named.)}

 \item{na.action}{A function to specify the action to be taken if \code{NA}s are
   found. The default action is \code{na.omit}, which leads to
   rejection of cases with missing values on any required variable. An
   alternative is \code{na.fail}, which causes an error if \code{NA}
   cases are found. (NOTE: If given, this argument must be named.)}

 \item{\dots}{additional parameters.}
}
 
\details{In quantile regression a function is fitted to the data so that
  it satisfies the property that a portion \eqn{tau} of the data
  \eqn{y|n} is below the estimate. While the error bars of many
  regression problems can be viewed as such estimates quantile
  regression estimates this quantity directly. Kernel quantile regression
  is similar to nu-Support Vector Regression in that it minimizes a
  regularized loss function in RKHS. The difference between nu-SVR and
  kernel quantile regression is in the type of loss function used which
  in the case of quantile regression is the pinball loss (see reference
  for details.). Minimizing the regularized loss boils down to a
  quadratic problem which is solved using an interior point QP solver
  \code{ipop} implemented in \code{kernlab}.
  
}
\value{
  An S4 object of class \code{kqr} containing the fitted model along with
  information.Accessor functions can be used to access the slots of the
  object which include :
  \item{alpha}{The resulting model parameters which can be also accessed
    by \code{coef}.}
  \item{kernelf}{the kernel function used.}
  \item{error}{Training error (if fit == TRUE)}
  see \code{kqr-class} for more details.
}

  \references{Ichiro Takeuchi, Quoc V. Le, Timothy D. Sears, Alexander J. Smola\cr
    \emph{Nonparametric Quantile Estimation}\cr
    Journal of Machine Learning Research 7,2006,1231-1264  \cr
     \url{http://www.jmlr.org/papers/volume7/takeuchi06a/takeuchi06a.pdf}
   }
   
  \author{Alexandros Karatzoglou \cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}}
  
\seealso{\code{\link{predict.kqr}}, \code{\link{kqr-class}}, \code{\link{ipop}}, \code{\link{rvm}}, \code{\link{ksvm}}}

\examples{
# create data
x <- sort(runif(300))
y <- sin(pi*x) + rnorm(300,0,sd=exp(sin(2*pi*x)))

# first calculate the median
qrm <- kqr(x, y, tau = 0.5, C=0.15)

# predict and plot
plot(x, y)
ytest <- predict(qrm, x)
lines(x, ytest, col="blue")

# calculate 0.9 quantile
qrm <- kqr(x, y, tau = 0.9, kernel = "rbfdot",
           kpar= list(sigma=10), C=0.15)
ytest <- predict(qrm, x)
lines(x, ytest, col="red")

# calculate 0.1 quantile
qrm <- kqr(x, y, tau = 0.1,C=0.15)
ytest <- predict(qrm, x)
lines(x, ytest, col="green")

# print first 10 model coefficients
coef(qrm)[1:10]
}

\keyword{regression}
\keyword{nonlinear}
\keyword{methods}
