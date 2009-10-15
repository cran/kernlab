\name{kmmd-class}
\docType{class}
\alias{kmmd-class}
\alias{kernelf,kmmd-method}
\alias{H0,kmmd-method}
\alias{AsympH0,kmmd-method}
\alias{Radbound,kmmd-method}
\alias{Asymbound,kmmd-method}
\alias{mmdstats,kmmd-method}

\title{Class "kqr"}
\description{The Kernel Maximum Mean Discrepancy object class}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("kmmd", ...)}.
    or by calling the \code{kmmd} function 
}
\section{Slots}{
  \describe{
    \item{\code{kernelf}:}{Object of class \code{"kfunction"} contains
      the kernel function used}
    \item{\code{xmatrix}:}{Object of class \code{"kernelMatrix"} containing the
      data used }
    \item{H0}{Object of class \code{"logical"} contains value of : is H0 rejected (logical)}
    \item{\code{AsympH0}}{Object of class \code{"logical"} contains
      value : is H0 rejected according to the asymptotic bound (logical)}
      \item{\code{mmdstats}}{Object of class \code{"vector"} contains the test statistics (vector of two)}
    \item{\code{Radbound}}{Object of class \code{"numeric"} contains the Rademacher bound}
    \item{\code{Asymbound}}{Object of class \code{"numeric"} contains the asymptotic bound}
    
  }
}
\section{Methods}{
  \describe{ 
      \item{kernelf}{\code{signature(object = "kmmd")}: returns the
      kernel function used}
    \item{H0}{\code{signature(object = "kmmd")}: returns the value of H0
      being rejected}
     \item{AsympH0}{\code{signature(object = "kmmd")}: returns the value of H0
      being rejected according to the asymptotic bound}
    \item{mmdstats}{\code{signature(object = "kmmd")}: returns the values
      of the mmd statistics}
    \item{Radbound}{\code{signature(object = "kmmd")}: returns the
      value of the Rademacher bound}
    \item{Asymbound}{\code{signature(object = "kmmd")}: returns the
      value of the asymptotic bound}
      }
}

\author{Alexandros Karatzoglou\cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}}

\seealso{
  \code{\link{kmmd}}, 
}
\examples{
# create data
x <- matrix(runif(300),100)
y <- matrix(runif(300)+1,100)


mmdo <- kmmd(x, y)

H0(mmdo)

}
\keyword{classes}
