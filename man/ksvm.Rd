\name{ksvm}
\alias{ksvm}
\alias{ksvm,formula-method}
\alias{ksvm,vector-method}
\alias{ksvm,matrix-method}
\alias{show,ksvm-method}
\alias{coeff,ksvm-method}
\title{Support Vector Machines}
\description{
  Support Vector Machines are an excellent tool
  for classification as well as regression. \code{ksvm} supports the
  classical C-classification, nu-classification, one-class-classification
  e-regression, nu-regression along with the Crammer-Singer method
  for multi-class classification spoc-classification.\cr
  A probabilistic prediction function for classification is also included.
}
\usage{
\S4method{ksvm}{formula}(x, data = NULL, ..., subset, na.action = na.omit, scaled = TRUE)

\S4method{ksvm}{vector}{x, ...}

\S4method{ksvm}{matrix}(x, y = NULL, scaled = TRUE, type = NULL, kernel ="rbfdot", kpar = list(sigma = 0.1),
C = 1, nu = 0.2, epsilon = 0.1, prob.model = FALSE, class.weights = NULL, cachesize = 40, tol = 0.001,
shrinking = TRUE, cross = 0, fit = TRUE, ..., subset, na.action = na.omit)
}

\arguments{
  \item{x}{a symbolic description of the model to be fit. Note, that an
          intercept is always included, whether given in the formula or
          not. When not using a formula x is a matrix or vector
	  containg the variables in the model}
  \item{data}{an optional data frame containing the variables in the model.
          By default the variables are taken from the environment which
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
    a factor or not, the default setting for \code{type} is \code{C-classification} or \code{eps-regression}, respectively, but can be overwritten by setting an explicit value.\cr
    Valid options are:
    \itemize{
      \item \code{C-classification}
      \item \code{nu-classification}
      \item \code{spoc-classification} (Crammer Singer multi-class)
      \item \code{one-classification} (novelty detection)
      \item \code{eps-regression}
      \item \code{nu-regression}
    }
  }
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
    kpar parameter as well. In the case of a Radial Basis kernel function (Gaussian)
  kpar can also be set to the string "automatic" which uses the heuristics in 
  \code{sigest} to calculate a good \code{sigma} value for the data.}
  \item{C}{cost of constraints violation (default: 1)---it is the
    `C'-constant of the regularization term in the Lagrange formulation.}
  \item{nu}{parameter needed for \code{nu-classification},
    \code{one-classification}, and \code{nu-regression}. The \code{nu}
    parameter sets the upper bound on the training error and the lower
    bound on the fraction of data points to become Support Vectors (default: 0.2).}
  \item{epsilon}{epsilon in the insensitive-loss function used for
    \code{eps-regression} and \code{nu-regression} (default: 0.1)}
  \item{prob.model}{if set to \code{TRUE} a model for calculating class
    probabilities is fitted on output data created by performing a
    3-fold cross-validation on the training data. For details see
    references. (default: \code{FALSE})}
  \item{class.weights}{a named vector of weights for the different
    classes, used for asymmetric class sizes. Not all factor levels have
    to be supplied (default weight: 1). All components have to be named.}
  \item{cachesize}{cache memory in MB (default 40)}
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
  \item{coefs}{The corresponding coefficients times the training labels.}
  \item{b}{The negative intercept.}
  \item{nSV}{The number of Support Vectors}
  \item{error}{Training error}
  \item{cross}{Cross validation error, (when cross > 0)}
}

  
\details{For multiclass-classification with k levels, k>2, \code{ksvm} uses the
  `one-against-one'-approach, in which k(k-1)/2 binary classifiers are
  trained; the appropriate class is found by a voting scheme.\cr
  If the predictor variables include factors, the formula interface must be used to get a
  correct model matrix.
  The predict function can return probabilistic output (probability
  matrix) in the case of
  classification by setting the \code{type} parameter to "probabilities".
  
}
\note{Data is scaled internally, usually yielding better results.}
\references{
  \itemize{
    \item
      Chang, Chih-Chung and Lin, Chih-Jen:\cr
      \emph{LIBSVM: a library for Support Vector Machines}\cr
      \url{http://www.csie.ntu.edu.tw/~cjlin/libsvm}

     \item
      Exact formulations of models, algorithms, etc. can be found in the
      document:\cr
      Chang, Chih-Chung and Lin, Chih-Jen:\cr
      \emph{LIBSVM: a library for Support Vector Machines}\cr
      \url{http://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.ps.gz}
   
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
  
  }
}
\author{
  Alexandros Karatzoglou (SMO optimizers in C/C++ by Chih-Chung Chang & Chih-Jen Lin)\cr
  \email{alexandros.karatzoglou@ci.tuwien.ac.at}
}
\seealso{\code{\link{predict.ksvm}}, \code{\link{couple}} }
  
\keyword{methods}
\keyword{regression}
\keyword{nonlinear}
\keyword{classif}

\examples{

## simple example using the spam data set
data(spam)


## create test and training set
index <- sample(1:dim(spam)[1])
spamtrain <- spam[index[1:floor(2 * dim(spam)[1]/3)], ]
spamtest <- spam[index[((2 * ceiling(dim(spam)[1]/3)) + 1):dim(spam)[1]], ]

## train a support vector machine
filter <- ksvm(type~.,data=spamtrain,kernel="rbfdot",kpar=list(sigma=0.05),C=5,cross=3)
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

## train a support vector machine
irismodel <- ksvm(Species~.,data=iris,kernel=rbf,C=10,prob.model=TRUE)

irismodel

## get fitted values
fit(irismodel)

## Test on the training set with probabilities as output
predict(irismodel, iris[,-5], type="probabilities")


## regression
# create data
x <- seq(-20,20,0.1)
y <- sin(x)/x + rnorm(401,sd=0.03)

# train support vector machine
regm <- ksvm(x,y,epsilon=0.01,kpar=list(sigma=16),cross=3)
plot(x,y,type="l")
lines(x,predict(regm,x),col="red")
}
