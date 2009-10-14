\name{ipop-class}
\docType{class}
\alias{ipop-class}
\alias{primal,ipop-method}
\alias{dual,ipop-method}
\alias{how,ipop-method}
\alias{primal}
\alias{dual}
\alias{how}

\title{Class "ipop"}
\description{The quadratic problem solver class}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ipop", ...)}.
   or by calling the \code{ipop} function.
}
\section{Slots}{
  \describe{
    \item{\code{primal}:}{Object of class \code{"vector"} the primal
      solution of the problem}
    \item{\code{dual}:}{Object of class \code{"numeric"} the dual of the
    problem}
    \item{\code{how}:}{Object of class \code{"character"} convergence information}
  }
}
\section{Methods}{
\describe{
  \item{primal}{Object of class \code{ipop}}{Return the primal of the problem}
  \item{dual}{Object of class \code{ipop}}{Return the dual of the problem}
  \item{how}{Object of class \code{ipop}}{Return information on convergence}
  }
}
  \author{Alexandros Karatzoglou\cr
    \email{alexandros.karatzoglou@ci.tuwien.ac.at}}

\seealso{
  \code{\link{ipop}}

}
\examples{
## solve the Support Vector Machine optimization problem
data(spam)

## sample a scaled part (300 points) of the spam data set
m <- 300
set <- sample(1:dim(spam)[1],m)
x <- scale(as.matrix(spam[,-58]))[set,]
y <- as.integer(spam[set,58])
y[y==2] <- -1

##set C parameter and kernel
C <- 5
rbf <- rbfdot(sigma = 0.1)

## create H matrix etc.
H <- kernelPol(rbf,x,,y)
c <- matrix(rep(-1,m))
A <- t(y)
b <- 0
l <- matrix(rep(0,m))
u <- matrix(rep(C,m))
r <- 0

sv <- ipop(c,H,A,b,l,u,r)
primal(sv)
dual(sv)
how(sv)

}
\keyword{classes}
