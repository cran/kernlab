\name{gausspr-class}
\docType{class}
\alias{gausspr-class}
\alias{alpha,gausspr-method}
\alias{cross,gausspr-method}
\alias{error,gausspr-method}
\alias{fit,gausspr-method}
\alias{kcall,gausspr-method}
\alias{kernelf,gausspr-method}
\alias{kpar,gausspr-method}
\alias{lev,gausspr-method}
\alias{type,gausspr-method}
\alias{alphaindex,gausspr-method}
\alias{xmatrix,gausspr-method}
\alias{ymatrix,gausspr-method}

\title{Class "gausspr"}
\description{The Gaussian Processes object }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gausspr", ...)}.
    or by calling the \code{gausspr} function 
}
\section{Slots}{
  \describe{
    \item{\code{tol}:}{Object of class \code{"numeric"} contains
      tolerance of termination criteria}
    \item{\code{kernelf}:}{Object of class \code{"function"} contains
      the kernel function used}
    \item{\code{kpar}:}{Object of class \code{"list"} contains the
      kernel parameter used }
    \item{\code{kcall}:}{Object of class \code{"ANY"} contains the used
      function call }
    \item{\code{type}:}{Object of class \code{"character"} contains
      type of problem }
    \item{\code{kterms}:}{Object of class \code{"ANY"} contains the
      terms representation of the symbolic model used (when using a formula)}}
    \item{\code{xmatrix}:}{Object of class \code{"matrix"} containing
      the data matrix used }
    \item{\code{ymatrix}:}{Object of class \code{"ANY"} containing the
      response matrix}
    \item{\code{fit}:}{Object of class \code{"ANY"} containing the
      fitted values  }
    \item{\code{lev}:}{Object of class \code{"vector"} containing the
      levels of the response (in case of classification) }
    \item{\code{nclass}:}{Object of class \code{"numeric"} containing
      the number of classes (in case of classification) }
    \item{\code{alpha}:}{Object of class \code{"ANY"} containing the
      computes alpha values }
    \item{\code{alphaindex}}{Object of class \code{"list"} containing
      the indexes for the alphas in various classes (in multi-class problems).}
    \item{\code{nvar}:}{Object of class \code{"numeric"} containing the
      computed variance}
    \item{\code{error}:}{Object of class \code{"numeric"} containing the
    training error}
    \item{\code{cross}:}{Object of class \code{"numeric"} containing the
    cross validation error}
    \item{\code{n.action}:}{Object of class \code{"ANY"} containing the
      action performed in NA }
  }
}
\section{Methods}{
  \describe{ 
    \item{alpha}{\code{signature(object = "gausspr")}: returns the alpha
    vector}
  \item{cross}{\code{signature(object = "gausspr")}: returns the cross
    validation error }
    \item{error}{\code{signature(object = "gausspr")}: returns the
      training error }
    \item{fit}{\code{signature(object = "gausspr")}: returns the fitted values }
    \item{kcall}{\code{signature(object = "gausspr")}: returns the call performed}
    \item{kernelf}{\code{signature(object = "gausspr")}: returns the
      kernel function used}
    \item{kpar}{\code{signature(object = "gausspr")}: returns the kernel
    parameter used}
    \item{lev}{\code{signature(object = "gausspr")}: returns the
      response levels (in classification) }
    \item{type}{\code{signature(object = "gausspr")}: returns the type
      of problem}
    \item{xmatrix}{\code{signature(object = "gausspr")}: returns the
      data matrix used}
    \item{ymatrix}{\code{signature(object = "gausspr")}: returns the
      response matrix used}
  }
}

\author{Alexandros Karatzoglou\cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}}



\seealso{
  \code{\link{gausspr}}, 
  \code{\link{ksvm-class}}
}
\examples{

# train model
data(iris)
test <- gausspr(Species~.,data=iris,var=2)
test
alpha(test)
error(test)
lev(test)
}
\keyword{classes}
