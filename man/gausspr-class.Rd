\name{gausspr-class}
\docType{class}
\alias{gausspr-class}
\alias{alpha,gausspr-method}
\alias{cross,gausspr-method}
\alias{error,gausspr-method}
\alias{kcall,gausspr-method}
\alias{kernelf,gausspr-method}
\alias{kpar,gausspr-method}
\alias{lev,gausspr-method}
\alias{type,gausspr-method}
\alias{alphaindex,gausspr-method}
\alias{xmatrix,gausspr-method}
\alias{ymatrix,gausspr-method}
\alias{scaling,gausspr-method}

\title{Class "gausspr"}
\description{The Gaussian Processes object class}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("gausspr", ...)}.
    or by calling the \code{gausspr} function 
}
\section{Slots}{
  \describe{
    \item{\code{tol}:}{Object of class \code{"numeric"} contains
      tolerance of termination criteria}
    \item{\code{kernelf}:}{Object of class \code{"kfunction"} contains
      the kernel function used}
    \item{\code{kpar}:}{Object of class \code{"list"} contains the
      kernel parameter used }
    \item{\code{kcall}:}{Object of class \code{"list"} contains the used
      function call }
    \item{\code{type}:}{Object of class \code{"character"} contains
      type of problem }
    \item{\code{terms}:}{Object of class \code{"ANY"} contains the
      terms representation of the symbolic model used (when using a formula)}
    \item{\code{xmatrix}:}{Object of class \code{"input"} containing
      the data matrix used }
    \item{\code{ymatrix}:}{Object of class \code{"output"} containing the
      response matrix}
    \item{\code{fitted}:}{Object of class \code{"output"} containing the
      fitted values  }
    \item{\code{lev}:}{Object of class \code{"vector"} containing the
      levels of the response (in case of classification) }
    \item{\code{nclass}:}{Object of class \code{"numeric"} containing
      the number of classes (in case of classification) }
    \item{\code{alpha}:}{Object of class \code{"listI"} containing the
      computes alpha values }
    \item{\code{alphaindex}}{Object of class \code{"list"} containing
      the indexes for the alphas in various classes (in multi-class
      problems).}
      \item{\code{sol}}{Object of class \code{"matrix"} containing the solution to the Gaussian Process formulation, it is used to compute the variance in regression problems.} 
     \item{\code{scaling}}{Object of class \code{"ANY"} containing
      the scaling coefficients of the data (when case \code{scaled = TRUE} is used).}
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
    \item{fitted}{\code{signature(object = "vm")}: returns the fitted values }
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
    \item{scaling}{\code{signature(object = "gausspr")}: returns the
      scaling coefficients of the data (when \code{scaled = TRUE} is used)}
  
  }
}

\author{Alexandros Karatzoglou\cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}}



\seealso{
  \code{\link{gausspr}}, 
  \code{\link{ksvm-class}},
  \code{\link{vm-class}}
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
