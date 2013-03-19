\name{predict.kqr}
\alias{predict.kqr}
\alias{predict,kqr-method}
\title{Predict method for kernel Quantile Regression object}


\description{Prediction of test data for kernel quantile regression}


\usage{
\S4method{predict}{kqr}(object, newdata)
}

\arguments{

  \item{object}{an S4 object of class \code{kqr} created by the
    \code{kqr} function}
  \item{newdata}{a data frame, matrix, or kernelMatrix containing new data}
}

\value{The value of the quantile given by the computed \code{kqr}
  model in a vector of length equal to the the rows of \code{newdata}.
   }

\author{Alexandros Karatzoglou\cr
  \email{alexandros.karatzoglou@ci.tuwien.ac.at}}
   
\keyword{methods}
\keyword{regression}


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
           kpar= list(sigma=10), C=0.15)
ytest <- predict(qrm, x)
lines(x, ytest, col="red")
}

