\name{rvm-class}
\docType{class}
\alias{rvm-class}
\alias{RVindex}
\alias{mlike}
\alias{nvar}
\alias{RVindex,rvm-method}
\alias{alpha,rvm-method}
\alias{cross,rvm-method}
\alias{error,rvm-method}
\alias{kcall,rvm-method}
\alias{kernelf,rvm-method}
\alias{kpar,rvm-method}
\alias{lev,rvm-method}
\alias{mlike,rvm-method}
\alias{nvar,rvm-method}
\alias{type,rvm-method}
\alias{xmatrix,rvm-method}
\alias{ymatrix,rvm-method}

\title{Class "rvm"}
\description{Relevance Vector Machine Class}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rvm", ...)}.
or by calling the \code{rvm} function.
}
\section{Slots}{
  \describe{

    \item{\code{tol}:}{Object of class \code{"numeric"} contains
      tolerance of termination criteria used.}
    \item{\code{kernelf}:}{Object of class \code{"kfunction"} contains
      the kernel function used }
    \item{\code{kpar}:}{Object of class \code{"list"} contains the
      hyperparameter used}
    \item{\code{kcall}:}{Object of class \code{"call"} contains the
      function call}
    \item{\code{type}:}{Object of class \code{"character"} contains type
    of problem}
    \item{\code{terms}:}{Object of class \code{"ANY"}  containing the
      terms representation of the symbolic model used (when using a
      formula interface)}
    \item{\code{xmatrix}:}{Object of class \code{"matrix"} contains the data
      matrix used during computation}
    \item{\code{ymatrix}:}{Object of class \code{"output"} contains the
      response matrix}
    \item{\code{fitted}:}{Object of class \code{"output"} with the fitted
      values, (predict on training set).}
    \item{\code{lev}:}{Object of class \code{"vector"} contains the
      levels of the response (in classification)}
    \item{\code{nclass}:}{Object of class \code{"numeric"} contains the
      number of classes (in classification)}
    \item{\code{alpha}:}{Object of class \code{"listI"} containing the the
    resulting alpha vector}
     \item{\code{coef}:}{Object of class \code{"ANY"} containing the the
    resulting model parameters}
    \item{\code{nvar}:}{Object of class \code{"numeric"} containing the
      calculated variance (in case of regression)}
    \item{\code{mlike}:}{Object of class \code{"numeric"} containing the
    computed maximum likelihood}
    \item{\code{RVindex}:}{Object of class \code{"vector"} containing
      the indexes of the resulting relevance vectors }
    \item{\code{nRV}:}{Object of class \code{"numeric"} containing the
      number of relevance vectors}
    \item{\code{cross}:}{Object of class \code{"numeric"} containing the
      resulting cross validation error }
    \item{\code{error}:}{Object of class \code{"numeric"} containing the
    training error}
    \item{\code{n.action}:}{Object of class \code{"ANY"} containing the
      action performed on NA}

  }
}
\section{Methods}{
  \describe{
    \item{RVindex}{\code{signature(object = "rvm")}: returns the index
      of the relevance vectors }
    \item{alpha}{\code{signature(object = "rvm")}: returns the resulting
    alpha vector}
    \item{cross}{\code{signature(object = "rvm")}: returns the resulting
    cross validation error}
    \item{error}{\code{signature(object = "rvm")}: returns the training
      error  }
    \item{fitted}{\code{signature(object = "vm")}: returns the fitted values }
    \item{kcall}{\code{signature(object = "rvm")}: returns the function call }
    \item{kernelf}{\code{signature(object = "rvm")}: returns the used
      kernel function }
    \item{kpar}{\code{signature(object = "rvm")}: returns the parameters
    of the kernel function}
    \item{lev}{\code{signature(object = "rvm")}: returns the levels of
      the response (in classification)}
    \item{mlike}{\code{signature(object = "rvm")}: returns the estimated
      maximum likelihood}
    \item{nvar}{\code{signature(object = "rvm")}: returns the calculated
    variance (in regression)}
    \item{type}{\code{signature(object = "rvm")}: returns the type of problem}
    \item{xmatrix}{\code{signature(object = "rvm")}: returns the data
      matrix used during computation}
    \item{ymatrix}{\code{signature(object = "rvm")}: returns the used response }
  }
}

\author{Alexandros Karatzoglou\cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}}


\seealso{
  \code{\link{rvm}}, 
   \code{\link{ksvm-class}}
}

\examples{

# create data
x <- seq(-20,20,0.1)
y <- sin(x)/x + rnorm(401,sd=0.05)

# train relevance vector machine
foo <- rvm(x, y)
foo

alpha(foo)
RVindex(foo)
fitted(foo)
kernelf(foo)
nvar(foo)

## show slots
slotNames(foo)

}
\keyword{classes}
