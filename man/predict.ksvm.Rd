\name{predict.ksvm}
\alias{predict.ksvm}
\alias{predict,ksvm-method}
\title{predict method for support vector object}


\description{Prediction of test data using support vector machines}


\usage{
\S4method{predict}{ksvm}(object, newdata, type = "response", coupler = "minpair")
}

\arguments{

  \item{object}{an S4 object of class \code{ksvm} created by the
    \code{ksvm} function}
  \item{newdata}{a data frame or matrix containing new data}
  \item{type}{one of \code{response}, \code{probabilities} ,\code{votes}
    indicating the type of output: predicted values, matrix of class
    probabilities, or matrix of vote counts.}
  \item{coupler}{Coupling method used in the multiclass case, can be one
    of \code{minpair} or \code{pkpd} (see reference for more details).}

}

\value{
   If \code{type(object)} is \code{C-classification},
     \code{nu-classification} or \code{spoc-classification}
     the vector returned depends on the argument \code{type}:
     
     \item{response}{predicted classes (the classes with majority vote).}

     \item{probabilities}{matrix of class probabilities (one column for each class and
       one row for each input).}

     \item{votes}{matrix of vote counts (one column for each class and one row
       for each new input)}
     
     If \code{type(object)} is \code{eps-regression},or
     \code{nu-regression} a vector of predicted values is returned.
     If \code{type(object)} is \code{one-classification} a vector of
     logical values is returned.
   }


   \references{
 \itemize{
     \item
       T.F. Wu, C.J. Lin, R.C. Weng. \cr
       \emph{Probability estimates for Multi-class Classification by
	 Pairwise Coupling}\cr
        \url{http://www.csie.ntu.edu.tw/~cjlin/papers/svmprob/svmprob.pdf}

      \item
	H.T. Lin, C.J. Lin, R.C. Weng\cr
	\emph{A note on Platt's probabilistic outputs for support vector
	  machines}\cr
	\url{http://www.csie.ntu.edu.tw/~cjlin/papers/plattprob.ps}
    }
}
\author{Alexandros Karatzoglou\cr
  \email{alexandros.karatzoglou@ci.tuwien.ac.at}}
   
\keyword{methods}
\keyword{regression}
\keyword{classif}


\examples{

## example using the spam data set
data(spam)
## create test and training set
traindex <- sample(1:dim(spam)[1],2*dim(spam)[1]/3)
spamtrain <- spam[traindex, ]
spamtest <- spam[-traindex, ]

## train a support vector machine
filter <- ksvm(type~.,data=spamtrain,kernel="rbfdot",kpar=list(sigma=0.05),C=5,cross=3)
filter

## predict mail type probabilities on the test set
mailtype <- predict(filter,spamtest,type="probabilities")
mailtype[c(1:10,1000:1010),]

}

