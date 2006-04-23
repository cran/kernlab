\name{inchol-class}
\docType{class}
\alias{inchol-class}
\alias{diagresidues}
\alias{maxresiduals}
\alias{pivots}
\alias{diagresidues,inchol-method}
\alias{maxresiduals,inchol-method}
\alias{pivots,inchol-method}

\title{Class "inchol" }
\description{The reduced Cholesky decomposition object}

\section{Objects from the Class}{Objects can be created by calls of the form \code{new("inchol", ...)}.
  or by calling the  \code{inchol} function.}

\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"} contains
      the decomposed matrix}

    \item{\code{pivots}:}{Object of class \code{"vector"} contains
      the pivots performed}

    \item{\code{diagresidues}:}{Object of class \code{"vector"} contains
      the diagonial residues}

    \item{\code{maxresiduals}:}{Object of class \code{"vector"} contains
    the maximum residues}
  }
}
\section{Extends}{
Class \code{"matrix"}, directly.
}
\section{Methods}{
  \describe{

    \item{diagresidues}{\code{signature(object = "inchol")}: returns
      the diagonial residues}

    \item{maxresiduals}{\code{signature(object = "inchol")}: returns
      the maximum residues}

    \item{pivots}{\code{signature(object = "inchol")}: returns
      the pivots performed}
  }
}

\author{Alexandros Karatzoglou\cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}}


\seealso{ \code{\link{inchol}}, \code{\link{csi-class}}, \code{\link{csi}}}

  \examples{
data(iris)
datamatrix <- as.matrix(iris[,-5])
# initialize kernel function
rbf <- rbfdot(sigma=0.1)
rbf
Z <- inchol(datamatrix,kernel=rbf)
dim(Z)
pivots(Z)
diagresidues(Z)
maxresiduals(Z)
}
\keyword{classes}
