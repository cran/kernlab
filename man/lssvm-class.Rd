\name{lssvm-class}
\docType{class}
\alias{lssvm-class}
\alias{alpha,lssvm-method}
\alias{b,lssvm-method}
\alias{cross,lssvm-method}
\alias{error,lssvm-method}
\alias{kcall,lssvm-method}
\alias{kernelf,lssvm-method}
\alias{kpar,lssvm-method}
\alias{param,lssvm-method}
\alias{lev,lssvm-method}
\alias{type,lssvm-method}
\alias{alphaindex,lssvm-method}
\alias{xmatrix,lssvm-method}
\alias{ymatrix,lssvm-method}
\alias{scaling,lssvm-method}
\alias{nSV,lssvm-method}

\title{Class "lssvm"}
\description{The Gaussian Processes object }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("lssvm", ...)}.
    or by calling the \code{lssvm} function 
}
\section{Slots}{
  \describe{
    \item{\code{kernelf}:}{Object of class \code{"kfunction"} contains
      the kernel function used}
    \item{\code{kpar}:}{Object of class \code{"list"} contains the
      kernel parameter used }
     \item{\code{param}:}{Object of class \code{"list"} contains the
      regularization parameter used.}
    \item{\code{kcall}:}{Object of class \code{"call"} contains the used
      function call }
    \item{\code{type}:}{Object of class \code{"character"} contains
      type of problem }
 \item{\code{coef}:}{Object of class \code{"ANY"} contains
      the model parameter }
    \item{\code{terms}:}{Object of class \code{"ANY"} contains the
      terms representation of the symbolic model used (when using a formula)}
    \item{\code{xmatrix}:}{Object of class \code{"matrix"} containing
      the data matrix used }
    \item{\code{ymatrix}:}{Object of class \code{"output"} containing the
      response matrix}
    \item{\code{fitted}:}{Object of class \code{"output"} containing the
      fitted values  }
\item{\code{b}:}{Object of class \code{"numeric"} containing the
      offset }
    \item{\code{lev}:}{Object of class \code{"vector"} containing the
      levels of the response (in case of classification) }
        \item{\code{scaling}:}{Object of class \code{"ANY"} containing the
      scaling information performed on the data}
    \item{\code{nclass}:}{Object of class \code{"numeric"} containing
      the number of classes (in case of classification) }
    \item{\code{alpha}:}{Object of class \code{"listI"} containing the
      computes alpha values }
    \item{\code{alphaindex}}{Object of class \code{"list"} containing
      the indexes for the alphas in various classes (in multi-class problems).}
    \item{\code{error}:}{Object of class \code{"numeric"} containing the
    training error}
    \item{\code{cross}:}{Object of class \code{"numeric"} containing the
    cross validation error}
    \item{\code{n.action}:}{Object of class \code{"ANY"} containing the
      action performed in NA }
    \item{\code{nSV}:}{Object of class \code{"numeric"} containing the
      number of model parameters }
  }
}
\section{Methods}{
  \describe{ 
    \item{alpha}{\code{signature(object = "lssvm")}: returns the alpha
    vector}
  \item{cross}{\code{signature(object = "lssvm")}: returns the cross
    validation error }
    \item{error}{\code{signature(object = "lssvm")}: returns the
      training error }
    \item{fitted}{\code{signature(object = "vm")}: returns the fitted values }
    \item{kcall}{\code{signature(object = "lssvm")}: returns the call performed}
    \item{kernelf}{\code{signature(object = "lssvm")}: returns the
      kernel function used}
    \item{kpar}{\code{signature(object = "lssvm")}: returns the kernel
      parameter used}
     \item{param}{\code{signature(object = "lssvm")}: returns the regularization
    parameter used}
    \item{lev}{\code{signature(object = "lssvm")}: returns the
      response levels (in classification) }
    \item{type}{\code{signature(object = "lssvm")}: returns the type
      of problem}
     \item{scaling}{\code{signature(object = "ksvm")}: returns the
      scaling values }
    \item{xmatrix}{\code{signature(object = "lssvm")}: returns the
      data matrix used}
    \item{ymatrix}{\code{signature(object = "lssvm")}: returns the
      response matrix used}
  }
}

\author{Alexandros Karatzoglou\cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}}



\seealso{
  \code{\link{lssvm}}, 
  \code{\link{ksvm-class}}
}
\examples{

# train model
data(iris)
test <- lssvm(Species~.,data=iris,var=2)
test
alpha(test)
error(test)
lev(test)
}
\keyword{classes}
