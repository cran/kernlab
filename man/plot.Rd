\name{plot}
\alias{plot.ksvm}
\alias{plot,ksvm,missing-method}
\alias{plot,ksvm-method}
\title{plot method for support vector object}


\description{Plot a binary classification support vector machine object.
The \code{plot} function returns a contour plot of the decision values. }


\usage{
\S4method{plot}{ksvm}(object, data=NULL, grid = 50, slice = list())
}

\arguments{

  \item{object}{a \code{ksvm} classification object created by the
    \code{ksvm} function}
  \item{data}{a data frame or matrix containing data to be plotted}
  \item{grid}{granularity for the contour plot.}
  \item{slice}{a list of named numeric values for the dimensions held
    constant (only needed if more than two variables are
    used). Dimensions not specified are fixed at 0. }

}

\seealso{\code{\link{ksvm}}}

\author{Alexandros Karatzoglou\cr
  \email{alexandros.karatzoglou@ci.tuwien.ac.at}}
   
\keyword{methods}
\keyword{regression}
\keyword{classif}


\examples{
## Demo of the plot function
x <- rbind(matrix(rnorm(120),,2),matrix(rnorm(120,mean=3),,2))
y <- matrix(c(rep(1,60),rep(-1,60)))

svp <- ksvm(x,y,type="C-svc")
plot(svp,data=x)

}

