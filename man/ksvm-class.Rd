\name{ksvm-class}
\docType{class}
\alias{ksvm-class}
\alias{SVindex}
\alias{alphaindex}
\alias{prob.model}
\alias{scaling}
\alias{prior}
\alias{show}
\alias{param}
\alias{b}
\alias{obj}
\alias{nSV}
\alias{coef,vm-method}
\alias{SVindex,ksvm-method}
\alias{alpha,ksvm-method}
\alias{alphaindex,ksvm-method}
\alias{cross,ksvm-method}
\alias{error,ksvm-method}
\alias{param,ksvm-method}
\alias{fitted,ksvm-method}
\alias{prior,ksvm-method}
\alias{prob.model,ksvm-method}
\alias{kernelf,ksvm-method}
\alias{kpar,ksvm-method}
\alias{lev,ksvm-method}
\alias{kcall,ksvm-method}
\alias{scaling,ksvm-method}
\alias{type,ksvm-method}
\alias{xmatrix,ksvm-method}
\alias{ymatrix,ksvm-method}
\alias{b,ksvm-method}
\alias{obj,ksvm-method}
\alias{nSV,ksvm-method}


\title{Class "ksvm" }
\description{An S4 class containing the output (model) of the
  \code{ksvm} Support Vector Machines function }
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ksvm", ...)}
  or by calls to the \code{ksvm} function.

}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}  containing
      the support vector machine type
      ("C-svc", "nu-svc", "C-bsvc", "spoc-svc",
      "one-svc", "eps-svr", "nu-svr", "eps-bsvr")}
    \item{\code{param}:}{Object of class \code{"list"} containing the
      Support Vector Machine parameters (C, nu, epsilon)}
    \item{\code{kernelf}:}{Object of class \code{"function"} containing
      the kernel function}
    \item{\code{kpar}:}{Object of class \code{"list"} containing the
      kernel function parameters (hyperparameters)}
    \item{\code{kcall}:}{Object of class \code{"ANY"} containing the      \code{ksvm} function call}
    \item{\code{scaling}:}{Object of class \code{"ANY"} containing the
      scaling information performed on the data}
    \item{\code{terms}:}{Object of class \code{"ANY"} containing the
      terms representation of the symbolic model used (when using a formula)}
    \item{\code{xmatrix}:}{Object of class \code{"input"} (\code{"list"}
      for multiclass problems 
      or \code{"matrix"} for binary classification and regression
      problems) containing the support vectors calculated from
      the data matrix used during computations (possibly scaled and
      without NA). In the case of multi-class classification each list
      entry contains the support vectors from each binary classification
    problem from the one-against-one method.}
    \item{\code{ymatrix}:}{Object of class \code{"output"}
      the response \code{"matrix"} or \code{"factor"} or \code{"vector"} or
      \code{"logical"}}
    \item{\code{fitted}:}{Object of class \code{"output"} with the fitted values,
      predictions using the training set.}
    \item{\code{lev}:}{Object of class \code{"vector"} with the levels of the
      response (in the case of classification)}
    \item{\code{prob.model}:}{Object of class \code{"list"} with the
      class prob. model}
    \item{\code{prior}:}{Object of class \code{"list"} with the
      prior of the training set}
    \item{\code{nclass}:}{Object of class \code{"numeric"}  containing
      the number of classes (in the case of classification)}
    \item{\code{alpha}:}{Object of class \code{"listI"} containing the
      resulting alpha vector (\code{"list"} or \code{"matrix"} in case of multiclass classification) (support vectors)}
    \item{\code{coef}:}{Object of class \code{"ANY"} containing the
      resulting coefficients}
    \item{\code{alphaindex}:}{Object of class \code{"list"} containing}
    \item{\code{b}:}{Object of class \code{"numeric"} containing the
      resulting offset }
    \item{\code{SVindex}:}{Object of class \code{"vector"} containing
      the indexes of the support vectors}
    \item{\code{nSV}:}{Object of class \code{"numeric"} containing the
      number of support vectors }
      \item{\code{obj}:}{Object of class \code{vector} containing the value of the objective function. When using 
      one-against-one in multiclass classification this is a vector.}
    \item{\code{error}:}{Object of class \code{"numeric"} containing the
    training error}
    \item{\code{cross}:}{Object of class \code{"numeric"} containing the
      cross-validation error }
    \item{\code{n.action}:}{Object of class \code{"ANY"} containing the
      action performed for NA }
  }
}
\section{Methods}{
  \describe{
    \item{SVindex}{\code{signature(object = "ksvm")}: return the indexes
    of support vectors}
    \item{alpha}{\code{signature(object = "ksvm")}: returns the complete
5    alpha vector (wit zero values)}
    \item{alphaindex}{\code{signature(object = "ksvm")}: returns the
      indexes of non-zero alphas (support vectors)}
    \item{cross}{\code{signature(object = "ksvm")}: returns the
      cross-validation error }
    \item{error}{\code{signature(object = "ksvm")}: returns the training
      error }
       \item{obj}{\code{signature(object = "ksvm")}: returns the value of the objective function}
    \item{fitted}{\code{signature(object = "vm")}: returns the fitted
      values (predict on training set) }
    \item{kernelf}{\code{signature(object = "ksvm")}: returns the kernel
    function}
    \item{kpar}{\code{signature(object = "ksvm")}: returns the kernel
      parameters (hyperparameters)}
    \item{lev}{\code{signature(object = "ksvm")}: returns the levels in
      case of classification  }
    \item{prob.model}{\code{signature(object="ksvm")}: returns class
      prob. model values}
    \item{param}{\code{signature(object="ksvm")}: returns 
      the parameters of the SVM in a list (C, epsilon, nu etc.)}
    \item{prior}{\code{signature(object="ksvm")}: returns 
      the prior of the training set}
    \item{kcall}{\code{signature(object="ksvm")}: returns the
    \code{ksvm} function call}
    \item{scaling}{\code{signature(object = "ksvm")}: returns the
      scaling values }
    \item{show}{\code{signature(object = "ksvm")}: prints the object information}
    \item{type}{\code{signature(object = "ksvm")}: returns the problem type}
    \item{xmatrix}{\code{signature(object = "ksvm")}: returns the data
      matrix used}
    \item{ymatrix}{\code{signature(object = "ksvm")}: returns the
      response vector}
  }
}

\author{Alexandros Karatzoglou \cr \email{alexandros.karatzolgou@ci.tuwien.ac.at}}


\seealso{
  \code{\link{ksvm}}, 
  \code{\link{rvm-class}},
  \code{\link{gausspr-class}}
}
\examples{
## simple example using the promotergene data set
data(promotergene)

## train a support vector machine
gene <- ksvm(Class~.,data=promotergene,kernel="rbfdot",
             kpar=list(sigma=0.015),C=50,cross=4)
gene

# the kernel  function
kernelf(gene)
# the alpha values
alpha(gene)
# the coefficients
coef(gene)
# the fitted values
fitted(gene)
# the cross validation error
cross(gene)


}
\keyword{classes}
