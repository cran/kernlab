\name{specc-class}
\docType{class}
\alias{specc-class}
\alias{centers}
\alias{size}
\alias{cluster}
\alias{centers,specc-method}
\alias{cluster,specc-method}
\alias{size,specc-method}
\alias{kernelf,specc-method}


\title{Class "specc"}
\description{ The Spectral Clustering Class}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("specc", ...)}.
    or by calling the function \code{specc}. 
}
\section{Slots}{
  \describe{
    \item{\code{cluster}:}{Object of class \code{"vector"} containing
      the cluster to which each point is allocated}
    \item{\code{centers}:}{Object of class \code{"matrix"} containing
      the cluser centers}
    \item{\code{size}:}{Object of class \code{"vector"} containing the
      number of points in each  cluster}
    \item{\code{kernelf}}{Object of class \code{kernel} containing the
      used kernel function.}
  }
}
\section{Methods}{
  \describe{
    \item{centers}{\code{signature(object = "specc")}: returns the
      cluster centers}
    \item{cluster}{\code{signature(object = "specc")}: returns the
      cluster to which each point is allocated}
    \item{size}{\code{signature(object = "specc")}: returns the number
      of points in each cluster }

  }
}

\author{Alexandros Karatzoglou\cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}}



\seealso{
  \code{\link{specc}}, 
   \code{\link{kpca-class}} 
}
\examples{
## Cluster the spirals data set.
data(spirals)

sc <- specc(spirals, centers=2, kernel="rbfdot", kpar=list(sigma=400))
cluster(sc)
centers(sc)
size(sc)
}
\keyword{classes}
