\name{kqr-class}
\docType{class}
\alias{kqr-class}
\alias{alpha,kqr-method}
\alias{cross,kqr-method}
\alias{error,kqr-method}
\alias{kcall,kqr-method}
\alias{kernelf,kqr-method}
\alias{kpar,kqr-method}
\alias{param,kqr-method}
\alias{alphaindex,kqr-method}
\alias{b,kqr-method}
\alias{xmatrix,kqr-method}
\alias{ymatrix,kqr-method}
\alias{scaling,kqr-method}

\title{Class "kqr"}
\description{The Kernel Quantile Regression object class}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("kqr", ...)}.
    or by calling the \code{kqr} function 
}
\section{Slots}{
  \describe{
    \item{\code{kernelf}:}{Object of class \code{"kfunction"} contains
      the kernel function used}
    \item{\code{kpar}:}{Object of class \code{"list"} contains the
      kernel parameter used }
    \item{\code{coef}:}{Object of class \code{"ANY"} containing the model parameters}
     \item{\code{param}:}{Object of class \code{"list"} contains the
      cost parameter C and tau parameter used }
    \item{\code{kcall}:}{Object of class \code{"list"} contains the used
      function call }
    \item{\code{terms}:}{Object of class \code{"ANY"} contains the
      terms representation of the symbolic model used (when using a formula)}
    \item{\code{xmatrix}:}{Object of class \code{"input"} containing
      the data matrix used }
    \item{\code{ymatrix}:}{Object of class \code{"output"} containing the
      response matrix}
    \item{\code{fitted}:}{Object of class \code{"output"} containing the
      fitted values  }
    \item{\code{alpha}:}{Object of class \code{"listI"} containing the
      computes alpha values }
    \item{\code{b}:}{Object of class \code{"numeric"} containing the
      offset of the model.}
     \item{\code{scaling}}{Object of class \code{"ANY"} containing
      the scaling coefficients of the data (when case \code{scaled = TRUE} is used).}
    \item{\code{error}:}{Object of class \code{"numeric"} containing the
    training error}
    \item{\code{cross}:}{Object of class \code{"numeric"} containing the
    cross validation error}
    \item{\code{n.action}:}{Object of class \code{"ANY"} containing the
      action performed in NA }
    \item{\code{nclass}:}{Inherited from class \code{vm}, not used in kqr}
    \item{\code{lev}:}{Inherited from class \code{vm}, not used in kqr}
    \item{\code{type}:}{Inherited from class \code{vm}, not used in kqr}
  }
}
\section{Methods}{
  \describe{ 
    \item{coef}{\code{signature(object = "kqr")}: returns the
      coefficients (alpha) of the model}
    \item{alpha}{\code{signature(object = "kqr")}: returns the alpha
      vector (identical to \code{coef})}
    \item{b}{\code{signature(object = "kqr")}: returns the offset beta
      of the model.}
  \item{cross}{\code{signature(object = "kqr")}: returns the cross
    validation error }
    \item{error}{\code{signature(object = "kqr")}: returns the
      training error }
    \item{fitted}{\code{signature(object = "vm")}: returns the fitted values }
    \item{kcall}{\code{signature(object = "kqr")}: returns the call performed}
    \item{kernelf}{\code{signature(object = "kqr")}: returns the
      kernel function used}
    \item{kpar}{\code{signature(object = "kqr")}: returns the kernel
    parameter used}
    \item{param}{\code{signature(object = "kqr")}: returns the
    cost regularization parameter C and tau used}
        \item{xmatrix}{\code{signature(object = "kqr")}: returns the
      data matrix used}
    \item{ymatrix}{\code{signature(object = "kqr")}: returns the
      response matrix used}
    \item{scaling}{\code{signature(object = "kqr")}: returns the
      scaling coefficients of the data (when \code{scaled = TRUE} is used)}
  
  }
}

\author{Alexandros Karatzoglou\cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}}

\seealso{
  \code{\link{kqr}}, 
  \code{\link{vm-class}},
  \code{\link{ksvm-class}}
}
\examples{


# create data
x <- sort(runif(300))
y <- sin(pi*x) + rnorm(300,0,sd=exp(sin(2*pi*x)))

# first calculate the median
qrm <- kqr(x, y, tau = 0.5, C=0.15)

# predict and plot
plot(x, y)
ytest <- predict(qrm, x)
lines(x, ytest, col="blue")

# calculate 0.9 quantile
qrm <- kqr(x, y, tau = 0.9, kernel = "rbfdot",
           kpar = list(sigma = 10), C = 0.15)
ytest <- predict(qrm, x)
lines(x, ytest, col="red")

# print model coefficients and other information
coef(qrm)
b(qrm)
error(qrm)
kernelf(qrm)
}
\keyword{classes}
