\name{couple}
\alias{couple}

\title{Probabilities Coupling function}
\description{
  \code{couple} is used to link class-probability estimates produced by
   pairwise coupling  in multi-class classification problems.
}
\usage{
couple(probin, coupler = "minpair")
}

\arguments{
  \item{probin}{ The pairwise coupled class-probability estimates}
  \item{coupler}{The type of coupler to use. Currently \code{minpar} and
    \code{pkpd} and \code{vote} are supported (see reference for more
    details).
  If \code{vote} is selected the returned value is a primitive estimate
  passed on given votes.}
  
}
\details{
  As binary classification problems are much easier to solve many
  techniques exist to decompose multi-class classification problems into
  many binary classification problems (voting, error codes,
  etc.). Pairwise coupling (one against one) constructs a rule for
  discriminating between every pair of classes and then selecting the
  class
  with the most winning two-class decisions. 
   By using Platt's probabilities output for SVM one can get a class
   probability for each of the \eqn{k(k-1)/2} models created in the pairwise 
   classification. The couple method implements various techniques to combine
   these probabilities.    
}
\value{
  A matrix with the resulting probability estimates. 
}
\references{
  Ting-Fan Wu, Chih-Jen Lin, ruby C. Weng\cr
  \emph{Probability Estimates for Multi-class Classification by Pairwise
    Coupling}\cr
  Neural Information Processing Symposium 2003 \cr
  \url{http://books.nips.cc/papers/files/nips16/NIPS2003_0538.pdf}
  }
\author{Alexandros Karatzoglou \cr \email{alexandros.karatzoglou@ci.tuwien.ac.at} }


 

\seealso{ \code{\link{predict.ksvm}}, \code{\link{ksvm}}}
\examples{
## create artificial pairwise probabilities
pairs <- matrix(c(0.82,0.12,0.76,0.1,0.9,0.05),2)

couple(pairs)

couple(pairs, coupler="pkpd")

couple(pairs, coupler ="vote")
}
\keyword{classif}

