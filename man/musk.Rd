\name{musk}
\alias{musk}
\docType{data}
\title{Musk data set}
\description{
This dataset describes a set of 92 molecules of which 47 are judged
   by human experts to be musks and the remaining 45 molecules are
   judged to be non-musks.
}
\usage{data(musk)}
\format{
  A data frame with 476 observations on the following 167 variables.

  Variables 1-162  are "distance features" along rays.  The distances are
measured in hundredths of Angstroms.  The distances may be negative or
positive, since they are actually measured relative to an origin placed
along each ray.  The origin was defined by a "consensus musk" surface
that is no longer used.  Hence, any experiments with the data should
treat these feature values as lying on an arbitrary continuous scale.  In
particular, the algorithm should not make any use of the zero point or
the sign of each feature value. 

Variable 163 is the distance of the oxygen atom in the molecule to a
designated point in 3-space. This is also called OXY-DIS.

Variable 164  is the X-displacement from the designated point.

Variable 165 is the Y-displacement from the designated point.

Variable 166 is the  Z-displacement from the designated point. 

Class: 0 for non-musk, and 1 for musk
}
   

\source{
 UCI Machine Learning data repository \cr
}

\examples{
data(musk)

muskm <- ksvm(Class~.,data=musk,kernel="rbfdot",C=1000)

muskm

}
\keyword{datasets}
