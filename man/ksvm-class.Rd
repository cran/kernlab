\name{ksvm-class}
\docType{class}
\alias{ksvm-class}
\alias{SVindex}
\alias{cross}
\alias{alpha}
\alias{alphaindex}
\alias{coeff}
\alias{cross}
\alias{error}
\alias{fit}
\alias{type}
\alias{kernelf}
\alias{xmatrix}
\alias{ymatrix}
\alias{scaling}
\alias{lev}
\alias{kcall}
\alias{show}
\alias{SVindex,ksvm-method}
\alias{alpha,ksvm-method}
\alias{alphaindex,ksvm-method}
\alias{cross,ksvm-method}
\alias{error,ksvm-method}
\alias{fit,ksvm-method}
\alias{kernelf,ksvm-method}
\alias{kpar,ksvm-method}
\alias{lev,ksvm-method}
\alias{kcall,ksvm-method}
\alias{scaling,ksvm-method}
\alias{type,ksvm-method}
\alias{xmatrix,ksvm-method}
\alias{ymatrix,ksvm-method}

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
      the problem support vector machine problem type
      ("C-classification", "nu-classification", "spoc-classification",
      "one-classification", "eps-regression", "nu-regression")}
    \item{\code{param}:}{Object of class \code{"list"} containing the
      Support Vector Machine parameters (C, nu, epsilon)}
    \item{\code{kernelf}:}{Object of class \code{"function"} containing
      the kernel function}
    \item{\code{kpar}:}{Object of class \code{"list"} containing the
      kernel function parameters (hyperparameters)}
    \item{\code{kcall}:}{Object of class \code{"ANY"} containing the      \code{ksvm} function call}
    \item{\code{scaling}:}{Object of class \code{"ANY"} containing the
      scaling information performed on the data}
    \item{\code{kterms}:}{Object of class \code{"ANY"} containing the
      terms representation of the symbolic model used (when using a formula)}
    \item{\code{xmatrix}:}{Object of class \code{"matrix"} the data
      matrix used during computations (possibly scaled and whithout NA)}
    \item{\code{ymatrix}:}{Object of class \code{"ANY"} the response matrix/vector }
    \item{\code{fit}:}{Object of class \code{"ANY"} with the fitted values,
      predictions using the training set.}
    \item{\code{lev}:}{Object of class \code{"vector"} with the levels of the
      response (in the case of classifiaction)
    }
    \item{\code{nclass}:}{Object of class \code{"numeric"}  containing
      the number of classes (in the case of classification)}
    \item{\code{alpha}:}{Object of class \code{"ANY"} containing the
      resulting alpha vector (list or matrix in case of multiclass classification) (support vectors)}
    \item{\code{coeff}:}{Object of class \code{"ANY"} containing the
      resulting coefficients}
    \item{\code{alphaindex}:}{Object of class \code{"list"} containing}
    \item{\code{b}:}{Object of class \code{"numeric"} containing the
      resulting offset }
    \item{\code{SVindex}:}{Object of class \code{"vector"} containing
      the indexes of the support vectors}
    \item{\code{nSV}:}{Object of class \code{"numeric"} containing the
      number of suppport vector machines }
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
    alpha vector (wit zero values)}
    \item{alphaindex}{\code{signature(object = "ksvm")}: returns the
      indexes of non-zero alphas (support vectors}
    \item{cross}{\code{signature(object = "ksvm")}: returns the
      cross-validation error }
    \item{error}{\code{signature(object = "ksvm")}: returns the training
      error }
    \item{fit}{\code{signature(object = "ksvm")}: returns the fitted
      values (predict on training set) }
    \item{kernelf}{\code{signature(object = "ksvm")}: returns the kernel
    function}
    \item{kpar}{\code{signature(object = "ksvm")}: returns the kernel
      parameters (hyperparameters)}
    \item{lev}{\code{signature(object = "ksvm")}: returns the levels in
      case of classification  }
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
## simple example using the spam data set
data(spam)

## create test and training set
spamtrain <- spam[1:(2 * dim(spam)[1]/3), ]
spamtest <- spam[((2 * dim(spam)[1]/3) + 1):length(spam), ]

## train a support vector machine
filter <- ksvm(type~.,data=spamtrain,kernel="rbfdot",kpar=list(sigma=0.05),C=5,cross=3)
filter

# the kernel  function
kernelf(filter)
# the alpha values
alpha(filter)
# the coefficients
coeff(filter)
# the fitted values
fit(filter)
# the cross validation error
cross(filter)


}
\keyword{classes}
