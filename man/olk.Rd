\name{olk}
\alias{olk}
\alias{olk,olk-method}

\title{Kernel Online Learning algorithms}
\description{
Online Kernel-based Learning algorithms for classification, novelty
detection, and regression.
}
\usage{
\S4method{olk}{olk}(obj, x, y = NULL, C = 0.5, r = 1e-4, epsilon = 0.005)
}

\arguments{
  \item{obj}{\code{obj} an object of class \code{olk} created by the
    initialization function \code{initOLK} containing the kernel to be
    used during learning and the parameters of the
    learned model}
  \item{x}{vector or matrix containing the data. Factors have
    to be numerically coded. If \code{x} is a matrix the code is
    run internally one sample at the time.}
  \item{y}{the class label in case of classification and must be -1 or
    +1. For regression the value of is any real number.
  }
  \item{nu}{the parameter similarly to the \code{nu} parameter in SVM
    bounds the training error.}
  \item{lambda}{the learning rate}
}
\details{
  The online algorithms are based on a simple stochastic gradient descent
  method in feature space.
  The state of the algorithm is stored in an object of class
  \code{olk} and has to be passed to the function at each iteration.
}
\value{
  The function returns an \code{S4} object of class \code{olk}
  containing the model parameters and the last fitted value which can be
  retrieved by the accessor method \code{fit}. The value returned in the
  classification and novelty detection problem is the decision function
  value ft.
  The accessor methods \code{alpha} returns the model parameters. 
}
\references{ Li, Guoqi, Changyun Wen, Zheng Guo Li, Aimin Zhang, Feng
  Yang, and Kezhi Mao \cr
  \emph{Model-based online learning with kernels} \cr
  Neural Networks and Learning Systems, IEEE Transactions on 24, no. 3
  (2013): 356-369. \cr

\author{Stephen Tridgell\cr
  \email{stephen.tridgell@sydney.edu.au}}


\seealso{\code{\link{initOLK}}}
\examples{

## create toy data set
x <- rbind(matrix(rnorm(100),,2),matrix(rnorm(100)+3,,2))
y <- matrix(c(rep(1,50),rep(-1,50)),,1)

## initialize olk object
on <- initOLK(2,kernel="rbfdot",kpar=list(sigma=0.2),
              type="classification")

ind <- sample(1:100,100)
## learn one data point at the time
for(i in ind)
on <- olk(on,x[i,],y[i],nu=0.03,lambda=0.1)

## or learn all the data 
on <- olk(on,x[ind,],y[ind],nu=0.03,lambda=0.1)

sign(predict(on,x))
}

\keyword{classif}
\keyword{neural}
\keyword{regression}
\keyword{ts}
