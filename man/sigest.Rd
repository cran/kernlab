\name{sigest}
\alias{sigest}
\alias{sigest,formula-method}
\alias{sigest,matrix-method}

\title{Hyperparameter estimation for the Gaussian Radial Basis kernel}
\description{
 Given a range of values for the "sigma" inverse width parameter in the Gaussian Radial Basis kernel
 for use with Support Vector Machines. The estimation is based on the
 data to be used.  
}
\usage{
\S4method{sigest}{formula}(x, data=NULL, frac = 0.5, na.action = na.omit, scaled = TRUE)
\S4method{sigest}{matrix}(x, frac = 0.5, scaled = TRUE, na.action = na.omit)
}

\arguments{
  \item{x}{a symbolic description of the model upon the estimation is
    based. When not using a formula x is a matrix or vector
    containing the data}
 \item{data}{an optional data frame containing the variables in the model.
          By default the variables are taken from the environment which
          `ksvm' is called from.}

\item{frac}{Fraction of data to use for estimation. By default a quarter
  of the data is used to estimate the range of the sigma hyperparameter.}
	
	\item{scaled}{A logical vector indicating the variables to be
    scaled. If \code{scaled} is of length 1, the value is recycled as
    many times as needed and all non-binary variables are scaled.
    Per default, data are scaled internally to zero mean and unit
    variance
    (since this the default action in \code{ksvm} as well). The center and scale
    values are returned and used for later predictions. }
  \item{na.action}{A function to specify the action to be taken if \code{NA}s are
          found. The default action is \code{na.omit}, which leads to rejection of cases
          with missing values on any required variable. An alternative
	  is \code{na.fail}, which causes an error if \code{NA} cases
	  are found. (NOTE: If given, this argument must be named.)}
 
}



\details{
\code{sigest} estimates the range of values for the sigma parameter
which would return good results when used with a Support Vector
Machine (\code{ksvm}). The estimation is based upon the 0.1 and 0.9 quantile 
of \eqn{\|x -x'\|^2}. Basically any value in between those two bounds will
produce good results.
}
\value{
  Returns a vector of length 3 defining the range (0.1 quantile, median
  and 0.9 quantile) of
  the sigma hyperparameter. 
  }
\references{ B. Caputo, K. Sim, F. Furesjo, A. Smola, \cr
\emph{Appearance-based object recognition using SVMs: which kernel should I use?}\cr
Proc of NIPS workshop on Statitsical methods for computational experiments in visual processing and computer vision, Whistler, 2002.
}
\author{Alexandros Karatzoglou \cr 
  \email{alexandros.karatzoglou@ci.tuwien.ac.at}}



\seealso{\code{\link{ksvm}}}
\examples{

## estimate good sigma values for promotergene
data(promotergene)
srange <- sigest(Class~.,data = promotergene)
srange

s <- srange[2]
s
## create test and training set
ind <- sample(1:dim(promotergene)[1],20)
genetrain <- promotergene[-ind, ]
genetest <- promotergene[ind, ]

## train a support vector machine
gene <- ksvm(Class~.,data=genetrain,kernel="rbfdot",
             kpar=list(sigma = s),C=50,cross=3)
gene

## predict gene type on the test set
promoter <- predict(gene,genetest[,-1])

## Check results
table(promoter,genetest[,1])
}
\keyword{classif}
\keyword{regression}
