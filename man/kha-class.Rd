\name{kha-class}
\docType{class}
\alias{kha-class}
\alias{eig,kha-method}
\alias{kcall,kha-method}
\alias{kernelf,kha-method}
\alias{pcv,kha-method}
\alias{xmatrix,kha-method}
\alias{eskm,kha-method}

\title{Class "kha"}
\description{ The Kernel Hebbian Algorithm class}
\section{Objects objects of class "kha"}{
Objects can be created by calls of the form \code{new("kha", ...)}.
    or by calling the \code{kha} function.
}
\section{Slots}{
  \describe{
    \item{\code{pcv}:}{Object of class \code{"matrix"} containing the
      principal component vectors }
    \item{\code{eig}:}{Object of class \code{"vector"} containing the
      corresponding normalization values}
    \item{\code{eskm}:}{Object of class \code{"vector"} containing the
      kernel sum}
    \item{\code{kernelf}:}{Object of class \code{"kfunction"} containing
      the kernel function used}
    \item{\code{kpar}:}{Object of class \code{"list"} containing the
      kernel parameters used }
    \item{\code{xmatrix}:}{Object of class \code{"matrix"} containing
      the data matrix used }
    \item{\code{kcall}:}{Object of class \code{"ANY"} containing the
      function call }
    \item{\code{n.action}:}{Object of class \code{"ANY"} containing the
      action performed on NA }
  }
}
\section{Methods}{
  \describe{

    \item{eig}{\code{signature(object = "kha")}: returns the
      normalization values }

    \item{kcall}{\code{signature(object = "kha")}: returns the
      performed call}
    \item{kernelf}{\code{signature(object = "kha")}: returns the used
      kernel function}
    \item{pcv}{\code{signature(object = "kha")}: returns the principal
      component vectors }
    \item{eskm}{\code{signature(object = "kha")}: returns the kernel sum}
    \item{predict}{\code{signature(object = "kha")}: embeds new data }
    \item{xmatrix}{\code{signature(object = "kha")}: returns the used
      data matrix }
  }
}

\author{Alexandros Karatzoglou\cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}}

\seealso{ \code{\link{kha}},
  \code{\link{ksvm-class}}, 
   \code{\link{kcca-class}} 
}

\examples{
# another example using the iris
data(iris)
test <- sample(1:50,20)

kpc <- kha(~.,data=iris[-test,-5], kernel="rbfdot",
           kpar=list(sigma=0.2),features=2, eta=0.001, maxiter=65)

#print the principal component vectors
pcv(kpc)
kernelf(kpc)
eig(kpc)
}
\keyword{classes}
