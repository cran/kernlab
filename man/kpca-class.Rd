\name{kpca-class}
\docType{class}
\alias{kpca-class}
\alias{rotated}
\alias{eig,kpca-method}
\alias{kcall,kpca-method}
\alias{kernelf,kpca-method}
\alias{pcv,kpca-method}
\alias{rotated,kpca-method}
\alias{xmatrix,kpca-method}

\title{Class "kpca"}
\description{ The Kernel Principal Components Analysis class}
\section{Objects of class "kpca"}{
Objects can be created by calls of the form \code{new("kpca", ...)}.
    or by calling the \code{kpca} function.
}
\section{Slots}{
  \describe{
    \item{\code{pcv}:}{Object of class \code{"matrix"} containing the
      principal component vectors }
    \item{\code{eig}:}{Object of class \code{"vector"} containing the
      corresponding eigenvalues}
    \item{\code{rotated}:}{Object of class \code{"matrix"} containing the
      projection of the data on the principal components}
    \item{\code{kernelf}:}{Object of class \code{"function"} containing
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

    \item{eig}{\code{signature(object = "kpca")}: returns the eigenvalues }

    \item{kcall}{\code{signature(object = "kpca")}: returns the
      performed call}
    \item{kernelf}{\code{signature(object = "kpca")}: returns the used
      kernel function}
    \item{pcv}{\code{signature(object = "kpca")}: returns the principal
      component vectors }
    \item{predict}{\code{signature(object = "kpca")}: embeds new data }
    \item{rotated}{\code{signature(object = "kpca")}: returns the
      projected data}
    \item{xmatrix}{\code{signature(object = "kpca")}: returns the used
      data matrix }
  }
}

\author{Alexandros Karatzoglou\cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}}

\seealso{
  \code{\link{ksvm-class}}, 
   \code{\link{kcca-class}} 
}
\examples{
# another example using the iris
data(iris)
test <- sample(1:50,20)

kpc <- kpca(~.,data=iris[-test,-5],kernel="rbfdot",
            kpar=list(sigma=0.2),features=2)

#print the principal component vectors
pcv(kpc)
rotated(kpc)
kernelf(kpc)
eig(kpc)
}
\keyword{classes}
