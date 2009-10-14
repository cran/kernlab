\name{prc-class}
\docType{class}
\alias{prc-class}

\alias{eig}
\alias{pcv}

\alias{eig,prc-method}
\alias{kcall,prc-method}
\alias{kernelf,prc-method}
\alias{pcv,prc-method}
\alias{xmatrix,prc-method}

\title{Class "prc"}
\description{Principal Components Class}
\section{Objects of class "prc"}{Objects from the class cannot be created directly but only contained
  in other classes.}

\section{Slots}{
  \describe{
    \item{\code{pcv}:}{Object of class \code{"matrix"} containing the
      principal component vectors }

    \item{\code{eig}:}{Object of class \code{"vector"} containing the
      corresponding eigenvalues}

    \item{\code{kernelf}:}{Object of class \code{"kfunction"} containing
      the kernel function used}

    \item{\code{kpar}:}{Object of class \code{"list"} containing the
      kernel parameters used }

    \item{\code{xmatrix}:}{Object of class \code{"input"} containing
      the data matrix used }

    \item{\code{kcall}:}{Object of class \code{"ANY"} containing the
      function call }

    \item{\code{n.action}:}{Object of class \code{"ANY"} containing the
      action performed on NA }
  }
}
\section{Methods}{
  \describe{

    \item{eig}{\code{signature(object = "prc")}: returns the eigenvalues }

    \item{kcall}{\code{signature(object = "prc")}: returns the
      performed call}

    \item{kernelf}{\code{signature(object = "prc")}: returns the used
      kernel function}

    \item{pcv}{\code{signature(object = "prc")}: returns the principal
      component vectors }

    \item{predict}{\code{signature(object = "prc")}: embeds new data }

    \item{xmatrix}{\code{signature(object = "prc")}: returns the used
      data matrix }
  }
}

\author{Alexandros Karatzoglou\cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}}

\seealso{
  \code{\link{kpca-class}},\code{\link{kha-class}}, \code{\link{kfa-class}} 
}

\keyword{classes}
