\name{inc.chol-class}
\docType{class}
\alias{inc.chol-class}
\alias{diag.residues}
\alias{maxresiduals}
\alias{pivots}
\alias{diag.residues,inc.chol-method}
\alias{maxresiduals,inc.chol-method}
\alias{pivots,inc.chol-method}

\title{Class "inc.chol" }
\description{  The reduced Cholesky decomposition object}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("inc.chol", ...)}.
  or by calling the  \code{chol.reduce} function.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"} contains the decomposed matrix}
    \item{\code{pivots}:}{Object of class \code{"vector"} contains the
      pivots performed}
    \item{\code{diag.residues}:}{Object of class \code{"vector"}
      contains the diagonial residues}
    \item{\code{maxresiduals}:}{Object of class \code{"vector"} contains
    the maximum residues}
  }
}
\section{Extends}{
Class \code{"matrix"}, directly.
}
\section{Methods}{
  \describe{
    \item{diag.residues}{\code{signature(object = "inc.chol")}: returns
      the diagonial residues }
    \item{maxresiduals}{\code{signature(object = "inc.chol")}: returns
      the maximum residues }
    \item{pivots}{\code{signature(object = "inc.chol")}: returns the
      pivots performed }
  }
}

\author{Alexandros Karatzoglou\cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}}


\seealso{
 \code{\link{chol.reduce}}
  }
\examples{
data(iris)
datamatrix <- as.matrix(iris[,-5])
# initialize kernel function
rbf <- rbfdot(sigma=0.1)
rbf
Z <- chol.reduce(datamatrix,kernel=rbf)
dim(Z)
pivots(Z)
diag.residues(Z)
maxresiduals(Z)
}
\keyword{classes}
