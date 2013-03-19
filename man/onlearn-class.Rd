\name{onlearn-class}
\docType{class}
\alias{onlearn-class}
\alias{alpha,onlearn-method}
\alias{b,onlearn-method}
\alias{buffer,onlearn-method}
\alias{fit,onlearn-method}
\alias{kernelf,onlearn-method}
\alias{kpar,onlearn-method}
\alias{predict,onlearn-method}
\alias{rho,onlearn-method}
\alias{rho}
\alias{show,onlearn-method}  
\alias{type,onlearn-method}
\alias{xmatrix,onlearn-method}
\alias{buffer}

\title{Class "onlearn"}
\description{ The class of objects used by the Kernel-based Online
  learning algorithms} 
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("onlearn", ...)}.
or by calls to the function \code{inlearn}.
}
\section{Slots}{
  \describe{
    \item{\code{kernelf}:}{Object of class \code{"function"} containing
      the used kernel function}
    \item{\code{buffer}:}{Object of class \code{"numeric"} containing
      the size of the buffer}
    \item{\code{kpar}:}{Object of class \code{"list"} containing the
      hyperparameters of the kernel function.}
    \item{\code{xmatrix}:}{Object of class \code{"matrix"} containing
      the data points (similar to support vectors) }
    \item{\code{fit}:}{Object of class \code{"numeric"} containing the
      decision function value of the last data point}
    \item{\code{onstart}:}{Object of class \code{"numeric"} used for indexing }
    \item{\code{onstop}:}{Object of class \code{"numeric"} used for indexing}
    \item{\code{alpha}:}{Object of class \code{"ANY"} containing the
      model parameters}
    \item{\code{rho}:}{Object of class \code{"numeric"} containing model
    parameter}
    \item{\code{b}:}{Object of class \code{"numeric"} containing the offset}
    \item{\code{pattern}:}{Object of class \code{"factor"} used for
      dealing with factors}
    \item{\code{type}:}{Object of class \code{"character"} containing
      the problem type (classification, regression, or novelty }
  }
}
\section{Methods}{
  \describe{
    \item{alpha}{\code{signature(object = "onlearn")}: returns the model
    parameters}
    \item{b}{\code{signature(object = "onlearn")}: returns the offset }
    \item{buffer}{\code{signature(object = "onlearn")}: returns the
      buffer size}
    \item{fit}{\code{signature(object = "onlearn")}: returns the last
      decision function value}
    \item{kernelf}{\code{signature(object = "onlearn")}: return the
      kernel function used}
    \item{kpar}{\code{signature(object = "onlearn")}: returns the
      hyper-parameters used}
    \item{onlearn}{\code{signature(obj = "onlearn")}: the learning function}
    \item{predict}{\code{signature(object = "onlearn")}: the predict function}
    \item{rho}{\code{signature(object = "onlearn")}: returns model parameter}
    \item{show}{\code{signature(object = "onlearn")}: show function}
    \item{type}{\code{signature(object = "onlearn")}: returns the type
      of problem}
    \item{xmatrix}{\code{signature(object = "onlearn")}: returns the
      stored data points}
  }
}

\author{Alexandros Karatzoglou\cr
\email{alexandros.karatzoglou@ci.tuwien.ac.at}}


\seealso{
  \code{\link{onlearn}}, \code{\link{inlearn}}
}
\examples{

## create toy data set
x <- rbind(matrix(rnorm(100),,2),matrix(rnorm(100)+3,,2))
y <- matrix(c(rep(1,50),rep(-1,50)),,1)

## initialize onlearn object
on <- inlearn(2,kernel="rbfdot",kpar=list(sigma=0.2),
              type="classification")

## learn one data point at the time
for(i in sample(1:100,100))
on <- onlearn(on,x[i,],y[i],nu=0.03,lambda=0.1)

sign(predict(on,x))

}
\keyword{classes}
