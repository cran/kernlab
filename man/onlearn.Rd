\name{onlearn}
\alias{onlearn}
\alias{onlearn,onlearn-method}

\title{Kernel Online Learning algorithms}
\description{
Online Kernel-based Learning algorithms for classification, novelty
detection, and regression. 
}
\usage{
\S4method{onlearn}{onlearn}(obj, x, y = NULL, nu = 0.2, lambda = 1e-04)
}

\arguments{
  \item{obj}{\code{obj} an object of class \code{onlearn} created by the
    initialization function \code{inlearn} containing the kernel to be
    used during learning and the parameters of the
    learned model}
  \item{x}{vector or matrix containing the data. Factors have
    to be numerically coded. If \code{x} is a matrix the code is
    run internally one sample at the time.}
  \item{y}{the class label in case of classification. Only binary
    classification is supported and class labels have to be -1 or +1.
  }
  \item{nu}{the parameter similarly to the \code{nu} parameter in SVM
    bounds the training error.}
  \item{lambda}{the learning rate}
}
\details{
  The online algorithms are based on a simple stochastic gradient descent
  method in feature space.
  The state of the algorithm is stored in an object of class
  \code{onlearn} and has to be passed to the function at each iteration.
}
\value{
  The function returns an \code{S4} object of class \code{onlearn}
  containing the model parameters and the last fitted value which can be
  retrieved by the accessor method \code{fit}. The value returned in the
  classification and novelty detection problem is the decision function
  value phi.
  The accessor methods \code{alpha} returns the model parameters. 
}
\references{ Kivinen J. Smola A.J. Williamson R.C. \cr
  \emph{Online Learning with Kernels}\cr
  IEEE Transactions on Signal Processing vol. 52, Issue 8, 2004\cr
  \url{http://mlg.anu.edu.au/~smola/papers/KivSmoWil03.pdf}}

\author{Alexandros Karatzoglou\cr
  \email{alexandros.karatzoglou@ci.tuwien.ac.at}}


\seealso{\code{\link{inlearn}}}
\examples{

## create toy data set
x <- rbind(matrix(rnorm(100),,2),matrix(rnorm(100)+3,,2))
y <- matrix(c(rep(1,50),rep(-1,50)),,1)

## initialize onlearn object
on <- inlearn(2,kernel="rbfdot",kpar=list(sigma=0.2),
              type="classification")

ind <- sample(1:100,100)
## learn one data point at the time
for(i in ind)
on <- onlearn(on,x[i,],y[i],nu=0.03,lambda=0.1)

## or learn all the data 
on <- onlearn(on,x[ind,],y[ind],nu=0.03,lambda=0.1)

sign(predict(on,x))
}

\keyword{classif}
\keyword{neural}
\keyword{regression}
\keyword{ts}
