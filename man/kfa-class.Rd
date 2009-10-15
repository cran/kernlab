\name{kfa-class}
\docType{class}
\alias{kfa-class}
\alias{alpha,kfa-method}
\alias{alphaindex,kfa-method}
\alias{kcall,kfa-method}
\alias{kernelf,kfa-method}
\alias{predict,kfa-method}
\alias{xmatrix,kfa-method}

\title{Class "kfa"}
\description{The class of the object returned by the Kernel Feature
  Analysis \code{kfa} function}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("kfa", ...)} or by
calling the \code{kfa} method. The objects contain the features along with the
alpha values. 
}
\section{Slots}{
  \describe{
    \item{\code{alpha}:}{Object of class \code{"matrix"} containing the
      alpha values }
    \item{\code{alphaindex}:}{Object of class \code{"vector"} containing
    the indexes of the selected feature}
    \item{\code{kernelf}:}{Object of class \code{"kfunction"} containing
      the kernel function used}
    \item{\code{xmatrix}:}{Object of class \code{"matrix"} containing
      the selected features}
    \item{\code{kcall}:}{Object of class \code{"call"} containing the
      \code{kfa} function call}
    \item{\code{terms}:}{Object of class \code{"ANY"} containing the
      formula terms}
  }
}
\section{Methods}{
  \describe{
    \item{alpha}{\code{signature(object = "kfa")}: returns the alpha values }
    \item{alphaindex}{\code{signature(object = "kfa")}: returns the
      index of the selected features}
    \item{kcall}{\code{signature(object = "kfa")}: returns the function call }
    \item{kernelf}{\code{signature(object = "kfa")}: returns the kernel
      function used }
    \item{predict}{\code{signature(object = "kfa")}: used to embed more
      data points to the feature base}
    \item{xmatrix}{\code{signature(object = "kfa")}: returns the
      selected features. }
  }
}

\author{Alexandros Karatzoglou\cr
\email{alexandros.karatzoglou@ci.tuwien.ac.at}}




\seealso{\code{\link{kfa}}, \code{\link{kpca-class}} }

\examples{
data(promotergene)
f <- kfa(~.,data=promotergene)
}
\keyword{classes}
