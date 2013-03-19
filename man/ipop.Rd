\name{ipop}
\alias{ipop}
\alias{ipop,ANY,matrix-method}

\title{Quadratic Programming Solver}
\description{
  ipop solves the quadratic programming problem :\cr
  \eqn{\min(c'*x + 1/2 * x' * H * x)}\cr
  subject to: \cr
  \eqn{b <= A * x <= b + r}\cr
  \eqn{l <= x <= u}
}
\usage{
ipop(c, H, A, b, l, u, r, sigf = 7, maxiter = 40, margin = 0.05,
     bound = 10, verb = 0)
}

\arguments{
  \item{c}{Vector or one column matrix appearing in the quadratic function}
  \item{H}{square matrix appearing in the quadratic function, or the
    decomposed form \eqn{Z} of the \eqn{H} matrix where \eqn{Z} is a
    \eqn{n x m} matrix with \eqn{n > m} and \eqn{ZZ' = H}.}
  \item{A}{Matrix defining the constrains under which we minimize the
    quadratic function}
  \item{b}{Vector or one column matrix defining the constrains}
  \item{l}{Lower bound vector or one column matrix}
  \item{u}{Upper bound vector or one column matrix}
  \item{r}{Vector or one column matrix defining constrains}
  \item{sigf}{Precision (default: 7 significant figures)}
  \item{maxiter}{Maximum number of iterations}
  \item{margin}{how close we get to the constrains}
  \item{bound}{Clipping bound for the variables}
  \item{verb}{Display convergence information during runtime}
}
\details{
  ipop uses an interior point method to solve the quadratic programming
  problem. \cr
  The \eqn{H} matrix can also be provided in the decomposed form \eqn{Z}
  where \eqn{ZZ' = H} in that case the Sherman Morrison Woodbury formula
  is used internally.
}
\value{
  An S4 object with the following slots
  \item{primal}{Vector containing the primal solution of the quadratic problem}
  \item{dual}{The dual solution of the problem}
  \item{how}{Character string describing the type of convergence}

  all slots can be accessed through accessor functions (see example)
}
\references{
  R. J. Vanderbei\cr
  \emph{LOQO: An interior point code for quadratic programming}\cr
  Optimization Methods and Software 11, 451-484, 1999 \cr
  \url{http://www.sor.princeton.edu/~rvdb/ps/loqo3.ps.gz}
}
\author{Alexandros Karatzoglou (based on Matlab code by Alex Smola) \cr
\email{alexandros.karatzoglou@ci.tuwien.ac.at}}


\seealso{\code{solve.QP}, \code{\link{inchol}}, \code{\link{csi}}}
\examples{
## solve the Support Vector Machine optimization problem
data(spam)

## sample a scaled part (500 points) of the spam data set
m <- 500
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
sv
dual(sv)

}
\keyword{optimize}

