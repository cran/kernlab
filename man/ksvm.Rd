\name{ksvm}
\alias{ksvm}
\alias{ksvm,formula-method}
\alias{ksvm,vector-method}
\alias{ksvm,matrix-method}
\alias{ksvm,kernelMatrix-method}
\alias{ksvm,list-method}
\alias{show,ksvm-method}
\alias{coef,ksvm-method}
\title{Support Vector Machines}
\description{
  Support Vector Machines are an excellent tool for classification, 
  novelty detection, and regression. \code{ksvm} supports the 
  well known C-svc, nu-svc, (classification) one-class-svc (novelty)
  eps-svr, nu-svr (regression) formulations along with 
  native multi-class classification formulations and 
 the bound-constraint SVM formulations.\cr
  \code{ksvm} also supports class-probabilities output and 
  confidence intervals for regression.
}
\usage{
\S4method{ksvm}{formula}(x, data = NULL, ..., subset, na.action = na.omit, scaled = TRUE)

\S4method{ksvm}{vector}(x, ...)

\S4method{ksvm}{matrix}(x, y = NULL, scaled = TRUE, type = NULL,
     kernel ="rbfdot", kpar = "automatic",
     C = 1, nu = 0.2, epsilon = 0.1, prob.model = FALSE,
     class.weights = NULL, cross = 0, fit = TRUE, cache = 40,
     tol = 0.001, shrinking = TRUE, ..., 
     subset, na.action = na.omit)

\S4method{ksvm}{kernelMatrix}(x, y = NULL, type = NULL,
     C = 1, nu = 0.2, epsilon = 0.1, prob.model = FALSE,
     class.weights = NULL, cross = 0, fit = TRUE, cache = 40,
     tol = 0.001, shrinking = TRUE, ...)

\S4method{ksvm}{list}(x, y = NULL, type = NULL,
     kernel = "stringdot", kpar = list(length = 4, lambda = 0.5),
     C = 1, nu = 0.2, epsilon = 0.1, prob.model = FALSE,
     class.weights = NULL, cross = 0, fit = TRUE, cache = 40,
     tol = 0.001, shrinking = TRUE, ...,
     na.action = na.omit)

}

\arguments{
  \item{x}{a symbolic description of the model to be fit.  When not
    using a formula x can be a matrix or vector containing the training
    data 
    or a kernel matrix of class \code{kernelMatrix} of the training data
    or a list of character vectors (for use with the string
    kernel). Note, that the intercept is always excluded, whether
    given in the formula or not.}

  \item{data}{an optional data frame containing the training data, when using a formula.
          By default the data is taken from the environment which
          `ksvm' is called from.}
	  
  \item{y}{a response vector with one label for each row/component of \code{x}. Can be either
    a factor (for classification tasks) or a numeric vector (for
    regression).}

  \item{scaled}{A logical vector indicating the variables to be
    scaled. If \code{scaled} is of length 1, the value is recycled as
    many times as needed and all non-binary variables are scaled.
    Per default, data are scaled internally (both \code{x} and \code{y}
    variables) to zero mean and unit variance. The center and scale
    values are returned and used for later predictions.}

  \item{type}{\code{ksvm} can be used for classification
    , for regression, or for novelty detection.
    Depending on whether \code{y} is
    a factor or not, the default setting for \code{type} is \code{C-svc}
    or \code{eps-svr},
    respectively, but can be overwritten by setting an explicit value.\cr
    Valid options are:

    \itemize{
      \item \code{C-svc}   C classification

      \item \code{nu-svc}  nu classification

      \item \code{C-bsvc}  bound-constraint svm classification

      \item \code{spoc-svc}  Crammer, Singer native multi-class

      \item \code{kbb-svc}  Weston, Watkins native multi-class

      \item \code{one-svc}  novelty detection

      \item \code{eps-svr}  epsilon regression

      \item \code{nu-svr}   nu regression

      \item \code{eps-bsvr}  bound-constraint svm regression
    }
  }

  \item{kernel}{the kernel function used in training and predicting.
    This parameter can be set to any function, of class kernel, which
    computes the inner product in feature space between two
    vector arguments (see \code{\link{kernels}}). \cr
    kernlab provides the most popular kernel functions
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
  
  \item{C}{cost of constraints violation (default: 1) this is the
    `C'-constant of the regularization term in the Lagrange
    formulation.}

  \item{nu}{parameter needed for \code{nu-svc},
    \code{one-svc}, and \code{nu-svr}. The \code{nu}
    parameter sets the upper bound on the training error and the lower
    bound on the fraction of data points to become Support Vectors (default: 0.2).}

  \item{epsilon}{epsilon in the insensitive-loss function used for
    \code{eps-svr}, \code{nu-svr} and \code{eps-bsvm} (default: 0.1)}

  \item{prob.model}{if set to \code{TRUE} builds a model for calculating class
    probabilities or in case of regression, calculates the scaling
    parameter of the Laplacian distribution fitted on the residuals.
    Fitting is done  on output data created by performing a
    3-fold cross-validation on the training data. For details see
    references. (default: \code{FALSE})}

  \item{class.weights}{a named vector of weights for the different
    classes, used for asymmetric class sizes. Not all factor levels have
    to be supplied (default weight: 1). All components have to be named.}

  \item{cache}{cache memory in MB (default 40)}

  \item{tol}{tolerance of termination criterion (default: 0.001)}
 
  \item{shrinking}{option whether to use the shrinking-heuristics
    (default: \code{TRUE})}

  \item{cross}{if a integer value k>0 is specified, a k-fold cross
    validation on the training data is performed to assess the quality
    of the model: the accuracy rate for classification and the Mean
    Squared Error for regression}

  \item{fit}{indicates whether the fitted values should be computed
    and included in the model or not (default: \code{TRUE})}

  \item{\dots}{additional parameters for the low level fitting function}

  \item{subset}{An index vector specifying the cases to be used in the
          training sample.  (NOTE: If given, this argument must be
          named.)}

   \item{na.action}{A function to specify the action to be taken if \code{NA}s are
          found. The default action is \code{na.omit}, which leads to rejection of cases
          with missing values on any required variable. An alternative
	  is \code{na.fail}, which causes an error if \code{NA} cases
	  are found. (NOTE: If given, this argument must be named.)}
      }
    
    \value{
      An S4 object of class \code{"ksvm"} containing the fitted model,
  Accessor functions can be used to access the slots of the object (see
  examples) which include:
  \item{alpha}{The resulting support vectors, (alpha vector) (possibly scaled).}
  \item{alphaindex}{The index of the resulting support vectors in the data
    matrix. Note that this index refers to the pre-processed data (after
    the possible effect of \code{na.omit} and \code{subset})}
  \item{coef}{The corresponding coefficients times the training labels.}
  \item{b}{The negative intercept.}
  \item{nSV}{The number of Support Vectors}
 \item{obj}{The value of the objective function. In case of one-against-one classification this is a vector of values} 
  \item{error}{Training error}
  \item{cross}{Cross validation error, (when cross > 0)}
  \item{prob.model}{Contains the width of the Laplacian fitted on the
    residuals in case of regression, or the parameters of the sigmoid
    fitted on the decision values in case of classification.}
}

  
\details{
  \code{ksvm} uses John Platt's SMO algorithm for solving the SVM QP problem an
  most SVM formulations. On the \code{spoc-svc}, \code{kbb-svc}, \code{C-bsvc} and
  \code{eps-bsvr} formulations a chunking algorithm based on the TRON QP
  solver is used. \cr
  For multiclass-classification with \eqn{k} classes, \eqn{k > 2}, \code{ksvm} uses the
  `one-against-one'-approach, in which \eqn{k(k-1)/2} binary classifiers are
  trained; the appropriate class is found by a voting scheme,
  The \code{spoc-svc} and the \code{kbb-svc} formulations deal with the
  multiclass-classification problems by solving a single quadratic problem involving all the classes.\cr
  If the predictor variables include factors, the formula interface must be used to get a
  correct model matrix. \cr
  In classification when \code{prob.model} is \code{TRUE} a 3-fold cross validation is
  performed on the data and a sigmoid function is fitted on the
  resulting decision values \eqn{f}.
  The data can be passed to the \code{ksvm} function in a \code{matrix} or a
  \code{data.frame}, in addition \code{ksvm} also supports input in the form of a
  kernel matrix of class \code{kernelMatrix} or as a list of character
  vectors where a string kernel has to be used.\cr
  The \code{plot} function for binary classification \code{ksvm} objects
  displays a contour plot of the decision values with the corresponding
  support vectors highlighted.\cr
 The predict function can return class probabilities for 
 classification problems by setting the \code{type} parameter to
 "probabilities". \cr
 The problem of model selection is partially addressed by an empirical
 observation for the RBF kernels (Gaussian , Laplace) where the optimal values of the
 \eqn{sigma} width parameter are shown to lie in between the 0.1 and 0.9
 quantile of the \eqn{\|x- x'\|} statistics. When using an RBF kernel
 and setting \code{kpar} to "automatic", \code{ksvm} uses the \code{sigest} function
 to estimate the quantiles and uses the median of the values.
}
\note{Data is scaled internally by default, usually yielding better results.}
\references{
  \itemize{
    \item
      Chang Chih-Chung, Lin Chih-Jen\cr
      \emph{LIBSVM: a library for Support Vector Machines}\cr
      \url{http://www.csie.ntu.edu.tw/~cjlin/libsvm}
   
   \item
      Chih-Wei Hsu, Chih-Jen Lin\cr
       \emph{BSVM}
       \url{http://www.csie.ntu.edu.tw/~cjlin/bsvm/}

     \item
     J. Platt\cr
     \emph{Probabilistic outputs for support vector machines and comparison to regularized likelihood methods} \cr
     Advances in Large Margin Classifiers, A. Smola, P. Bartlett, B. Schoelkopf and D. Schuurmans, Eds. Cambridge, MA: MIT Press, 2000.\cr
     \url{http://citeseer.nj.nec.com/platt99probabilistic.html}
     
    \item
    H.-T. Lin, C.-J. Lin and R. C. Weng\cr
    \emph{A note on Platt's probabilistic outputs for support vector machines}\cr
    \url{http://www.csie.ntu.edu.tw/~cjlin/papers/plattprob.ps}
    
    \item
     C.-W. Hsu and C.-J. Lin \cr
     \emph{A comparison on methods for multi-class support vector machines}\cr
     IEEE Transactions on Neural Networks, 13(2002) 415-425.\cr
     \url{http://www.csie.ntu.edu.tw/~cjlin/papers/multisvm.ps.gz}

     \item
     K. Crammer, Y. Singer\cr
     \emph{On the learnability and design of output codes for multiclass prolems}\cr
     Computational Learning Theory, 35-46, 2000.\cr
     \url{http://www.cs.huji.ac.il/~kobics/publications/mlj01.ps.gz}

     \item
     J. Weston, C. Watkins\cr
     \emph{Multi-class support vector machines} \cr
     In M. Verleysen, Proceedings of ESANN99 Brussels, 1999\cr
     \url{http://citeseer.ist.psu.edu/8884.html}
  }
}
\author{
  Alexandros Karatzoglou (SMO optimizers in C++ by Chih-Chung Chang & Chih-Jen Lin)\cr
  \email{alexandros.karatzoglou@ci.tuwien.ac.at}
}
\seealso{\code{\link{predict.ksvm}}, \code{\link{ksvm-class}}, \code{\link{couple}} }
  
\keyword{methods}
\keyword{regression}
\keyword{nonlinear}
\keyword{classif}
\keyword{neural}

\examples{

## simple example using the spam data set
data(spam)

## create test and training set
index <- sample(1:dim(spam)[1])
spamtrain <- spam[index[1:floor(dim(spam)[1]/2)], ]
spamtest <- spam[index[((ceiling(dim(spam)[1]/2)) + 1):dim(spam)[1]], ]

## train a support vector machine
filter <- ksvm(type~.,data=spamtrain,kernel="rbfdot",
               kpar=list(sigma=0.05),C=5,cross=3)
filter

## predict mail type on the test set
mailtype <- predict(filter,spamtest[,-58])

## Check results
table(mailtype,spamtest[,58])


## Another example with the famous iris data
data(iris)

## Create a kernel function using the build in rbfdot function
rbf <- rbfdot(sigma=0.1)
rbf

## train a bound constraint support vector machine
irismodel <- ksvm(Species~.,data=iris,type="C-bsvc",
                  kernel=rbf,C=10,prob.model=TRUE)

irismodel

## get fitted values
fitted(irismodel)

## Test on the training set with probabilities as output
predict(irismodel, iris[,-5], type="probabilities")


## Demo of the plot function
x <- rbind(matrix(rnorm(120),,2),matrix(rnorm(120,mean=3),,2))
y <- matrix(c(rep(1,60),rep(-1,60)))

svp <- ksvm(x,y,type="C-svc")
plot(svp,data=x)


### Use kernelMatrix
K <- as.kernelMatrix(crossprod(t(x)))

svp2 <- ksvm(K, y, type="C-svc")

svp2

# test data
xtest <- rbind(matrix(rnorm(20),,2),matrix(rnorm(20,mean=3),,2))
# test kernel matrix i.e. inner/kernel product of test data with
# Support Vectors

Ktest <- as.kernelMatrix(crossprod(t(xtest),t(x[SVindex(svp2), ])))

predict(svp2, Ktest)


#### Use custom kernel 

k <- function(x,y) {(sum(x*y) +1)*exp(-0.001*sum((x-y)^2))}
class(k) <- "kernel"

data(promotergene)

## train svm using custom kernel
gene <- ksvm(Class~.,data=promotergene[c(1:20, 80:100),],kernel=k,
             C=5,cross=5)

gene


#### Use text with string kernels
data(reuters)
is(reuters)
tsv <- ksvm(reuters,rlabels,kernel="stringdot",
            kpar=list(length=5),cross=3,C=10)
tsv


## regression
# create data
x <- seq(-20,20,0.1)
y <- sin(x)/x + rnorm(401,sd=0.03)

# train support vector machine
regm <- ksvm(x,y,epsilon=0.01,kpar=list(sigma=16),cross=3)
plot(x,y,type="l")
lines(x,predict(regm,x),col="red")
}
