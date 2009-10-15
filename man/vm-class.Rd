\name{vm-class}
\docType{class}

\alias{vm-class}
\alias{cross}
\alias{alpha}
\alias{error}
\alias{type}
\alias{kernelf}
\alias{xmatrix}
\alias{ymatrix}
\alias{lev}
\alias{kcall}

\alias{alpha,vm-method}
\alias{cross,vm-method}
\alias{error,vm-method}
\alias{fitted,vm-method}
\alias{kernelf,vm-method}
\alias{kpar,vm-method}
\alias{lev,vm-method}
\alias{kcall,vm-method}
\alias{type,vm-method}
\alias{xmatrix,vm-method}
\alias{ymatrix,vm-method}

\title{Class "vm" }
\description{An S4 VIRTUAL class used as a base for the various vector
  machine classes in \pkg{kernlab}}

\section{Objects from the Class}{
  Objects from the class cannot be created directly but only contained
  in other classes.
  }

  \section{Slots}{
  \describe{

    \item{\code{alpha}:}{Object of class \code{"listI"} containing the
      resulting alpha vector (list in case of multiclass classification) (support vectors)}

    \item{\code{type}:}{Object of class \code{"character"}  containing
      the vector machine type e.g.,
      ("C-svc", "nu-svc", "C-bsvc", "spoc-svc",
      "one-svc", "eps-svr", "nu-svr", "eps-bsvr")}

    \item{\code{kernelf}:}{Object of class \code{"function"} containing
      the kernel function}

    \item{\code{kpar}:}{Object of class \code{"list"} containing the
      kernel function parameters (hyperparameters)}

    \item{\code{kcall}:}{Object of class \code{"call"} containing the function call}

    \item{\code{terms}:}{Object of class \code{"ANY"} containing the
      terms representation of the symbolic model used (when using a formula)}

    \item{\code{xmatrix}:}{Object of class \code{"input"} the data
      matrix used during computations (support vectors) (possibly scaled and without NA)}

    \item{\code{ymatrix}:}{Object of class \code{"output"} the response matrix/vector }

    \item{\code{fitted}:}{Object of class \code{"output"} with the fitted values,
      predictions using the training set.}

    \item{\code{lev}:}{Object of class \code{"vector"} with the levels of the
      response (in the case of classification)}

     \item{\code{nclass}:}{Object of class \code{"numeric"}  containing
      the number of classes (in the case of classification)}
    
    \item{\code{error}:}{Object of class \code{"vector"} containing the
    training error}
  
   \item{\code{cross}:}{Object of class \code{"vector"} containing the
      cross-validation error }

    \item{\code{n.action}:}{Object of class \code{"ANY"} containing the
      action performed for NA }
  }
}
\section{Methods}{
  \describe{
    
    \item{alpha}{\code{signature(object = "vm")}: returns the complete
    alpha vector (wit zero values)}

    \item{cross}{\code{signature(object = "vm")}: returns the
      cross-validation error }

    \item{error}{\code{signature(object = "vm")}: returns the training
      error }

    \item{fitted}{\code{signature(object = "vm")}: returns the fitted
      values (predict on training set) }

    \item{kernelf}{\code{signature(object = "vm")}: returns the kernel
    function}

    \item{kpar}{\code{signature(object = "vm")}: returns the kernel
      parameters (hyperparameters)}

    \item{lev}{\code{signature(object = "vm")}: returns the levels in
      case of classification  }

    \item{kcall}{\code{signature(object="vm")}: returns the function call}

    \item{type}{\code{signature(object = "vm")}: returns the problem type}

    \item{xmatrix}{\code{signature(object = "vm")}: returns the data
      matrix used(support vectors)}

    \item{ymatrix}{\code{signature(object = "vm")}: returns the
      response vector}
  }
}

\author{Alexandros Karatzoglou \cr \email{alexandros.karatzolgou@ci.tuwien.ac.at}}


\seealso{
  \code{\link{ksvm-class}}, 
  \code{\link{rvm-class}},
  \code{\link{gausspr-class}}
}

\keyword{classes}
