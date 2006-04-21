\name{csi-class}
\docType{class}
\alias{csi-class}
\alias{Q}
\alias{R}
\alias{predgain}
\alias{truegain}
\alias{diagresidues,csi-method}
\alias{maxresiduals,csi-method}
\alias{pivots,csi-method}
\alias{predgain,csi-method}
\alias{truegain,csi-method}
\alias{Q,csi-method}
\alias{R,csi-method}

\title{Class "csi"}

\description{The reduced Cholesky decomposition object}

\section{Objects from the Class}{Objects can be created by calls of the form \code{new("csi", ...)}.
  or by calling the  \code{csi} function.}

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

    \item{predgain}{Object of class \code{"vector"} contains
      the predicted gain before adding each column}

    \item{truegain}{Object of class \code{"vector"} contains
      the actual gain after adding each column}

    \item{Q}{Object of class \code{"matrix"} contains
      Q from the QR decomposition of the kernel matrix}

    \item{R}{Object of class \code{"matrix"} contains
      R from the QR decomposition of the kernel matrix}
    
  }
}

\section{Extends}{
Class \code{"matrix"}, directly.
}
\section{Methods}{
  \describe{
    
    \item{diagresidues}{\code{signature(object = "csi")}: returns
      the diagonial residues}

    \item{maxresiduals}{\code{signature(object = "csi")}: returns
      the maximum residues}

    \item{pivots}{\code{signature(object = "csi")}: returns
      the pivots performed}

    \item{predgain}{\code{signature(object = "csi")}: returns
      the predicted gain before adding each column}

    \item{truegain}{\code{signature(object = "csi")}: returns
      the actual gain after adding each column}

    \item{Q}{\code{signature(object = "csi")}: returns
      Q from the QR decomposition of the kernel matrix}

    \item{R}{\code{signature(object = "csi")}: returns
      R from the QR decomposition of the kernel matrix}
  }
}

\author{Alexandros Karatzoglou\cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}}

\seealso{ \code{\link{csi}}, \code{\link{inchol-class}}}

\examples{
data(iris)

## create multidimensional y matrix
yind <- t(matrix(1:3,3,150))
ymat <- matrix(0, 150, 3)
ymat[yind==as.integer(iris[,5])] <- 1

datamatrix <- as.matrix(iris[,-5])
# initialize kernel function
rbf <- rbfdot(sigma=0.1)
rbf
Z <- csi(datamatrix,ymat, kernel=rbf, rank = 30)
dim(Z)
pivots(Z)
# calculate kernel matrix
K <- crossprod(t(Z))
# difference between approximated and real kernel matrix
(K - kernelMatrix(kernel=rbf, datamatrix))[6,]

}
\keyword{classes}
