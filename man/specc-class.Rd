\name{specc-class}
\docType{class}
\alias{specc-class}
\alias{centers}
\alias{size}
\alias{withinss}
\alias{centers,specc-method}
\alias{withinss,specc-method}
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
    \item{\code{.Data}:}{Object of class \code{"vector"} containing the cluster assignments}
    \item{\code{centers}:}{Object of class \code{"matrix"} containing
      the cluster centers}
    \item{\code{size}:}{Object of class \code{"vector"} containing the
      number of points in each  cluster}
    \item{\code{withinss}:}{Object of class \code{"vector"} containing
      the within-cluster sum of squares for each cluster}
    \item{\code{kernelf}}{Object of class \code{kernel} containing the
      used kernel function.}
  }
}
\section{Methods}{
  \describe{
    \item{centers}{\code{signature(object = "specc")}: returns the
      cluster centers}
    \item{withinss}{\code{signature(object = "specc")}: returns the
      within-cluster sum of squares for each cluster}
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

sc <- specc(spirals, centers=2)

centers(sc)
size(sc)
}
\keyword{classes}
