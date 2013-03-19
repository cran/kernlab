\name{promotergene}
\alias{promotergene}
\docType{data}
\title{E. coli promoter gene sequences (DNA)}
\description{
  Promoters have a region where a protein (RNA polymerase) must make contact
   and the helical DNA sequence must have a valid conformation so that
   the two pieces of the contact region spatially align.
   The data contains DNA sequences of promoters and non-promoters.
 }
\usage{data(promotergene)}
\format{
  A data frame with 106 observations and 58 variables.
  The first variable \code{Class} is a factor with levels \code{+} for a promoter gene  
  and \code{-} for a non-promoter gene. 
  The remaining 57 variables \code{V2 to V58} are factors describing the sequence. 
  The DNA bases are coded as follows: \code{a} adenine \code{c} cytosine \code{g} 
  guanine \code{t} thymine
}

\source{
  UCI Machine Learning data repository \cr
\url{ftp://ftp.ics.uci.edu/pub/machine-learning-databases/molecular-biology/promoter-gene-sequences}
}
\references{
  Towell, G., Shavlik, J. and Noordewier, M. \cr
  \emph{Refinement of Approximate Domain Theories by Knowledge-Based
  Artificial Neural Networks.} \cr
  In Proceedings of the Eighth National Conference on Artificial Intelligence (AAAI-90)
}
  

\examples{
data(promotergene)

## Create classification model using Gaussian Processes

prom <- gausspr(Class~.,data=promotergene,kernel="rbfdot",
                kpar=list(sigma=0.02),cross=4)
prom

## Create model using Support Vector Machines

promsv <- ksvm(Class~.,data=promotergene,kernel="laplacedot",
               kpar="automatic",C=60,cross=4)
promsv
}
\keyword{datasets}
