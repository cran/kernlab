\name{kcca-class}
\docType{class}
\alias{kcca-class}
\alias{kcor}
\alias{xcoef}
\alias{ycoef}
%%\alias{yvar}
%%\alias{xvar}
\alias{kcor,kcca-method}
\alias{xcoef,kcca-method}
\alias{xvar,kcca-method}
\alias{ycoef,kcca-method}
\alias{yvar,kcca-method}

\title{Class "kcca"}
\description{The "kcca" class  }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("kcca", ...)}.
    or by the calling the \code{kcca} function. 
}
\section{Slots}{
  \describe{
    \item{\code{kcor}:}{Object of class \code{"vector"} describing the correlations}
    \item{\code{xcoef}:}{Object of class \code{"matrix"} estimated coefficients for the \code{x} variables}
    \item{\code{ycoef}:}{Object of class \code{"matrix"} estimated coefficients for the \code{y} variables }
  %%  \item{\code{xvar}:}{Object of class \code{"matrix"} holds the
  %%    canonical variates for \code{x}}
  %%  \item{\code{yvar}:}{Object of class \code{"matrix"} holds the
  %%    canonical variates for \code{y}}
  }
}
\section{Methods}{
  \describe{
    \item{kcor}{\code{signature(object = "kcca")}: returns the correlations}
    \item{xcoef}{\code{signature(object = "kcca")}: returns the estimated coefficients for the \code{x} variables}
    \item{ycoef}{\code{signature(object = "kcca")}: returns the estimated coefficients for the \code{y} variables }
  %%  \item{xvar}{\code{signature(object = "kcca")}: returns the canonical
  %%  variates for \code{x}}
   %% \item{yvar}{\code{signature(object = "kcca")}: returns the canonical
  %%  variates for \code{y}}
  }
}

\author{Alexandros Karatzoglou \cr
\email{alexandros.karatzoglou@ci.tuwien.ac.at}}


\seealso{
  \code{\link{kcca}},
  \code{\link{kpca-class}} 
}
\examples{

## dummy data
x <- matrix(rnorm(30),15)
y <- matrix(rnorm(30),15)

kcca(x,y,ncomps=2)

}
\keyword{classes}
