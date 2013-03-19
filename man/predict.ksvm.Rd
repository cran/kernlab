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
  \item{type}{one of \code{response}, \code{probabilities}
    ,\code{votes}, \code{decision}
    indicating the type of output: predicted values, matrix of class
    probabilities, matrix of vote counts, or matrix of decision values.}
  \item{coupler}{Coupling method used in the multiclass case, can be one
    of \code{minpair} or \code{pkpd} (see reference for more details).}

}

\value{
   If \code{type(object)} is \code{C-svc},
     \code{nu-svc}, \code{C-bsvm} or \code{spoc-svc}
     the vector returned depends on the argument \code{type}:
     
     \item{response}{predicted classes (the classes with majority vote).}

     \item{probabilities}{matrix of class probabilities (one column for each class and
       one row for each input).}

     \item{votes}{matrix of vote counts (one column for each class and one row
       for each new input)}
     
     If \code{type(object)} is \code{eps-svr}, \code{eps-bsvr} or
     \code{nu-svr} a vector of predicted values is returned.
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

## example using the promotergene data set
data(promotergene)

## create test and training set
ind <- sample(1:dim(promotergene)[1],20)
genetrain <- promotergene[-ind, ]
genetest <- promotergene[ind, ]

## train a support vector machine
gene <- ksvm(Class~.,data=genetrain,kernel="rbfdot",
             kpar=list(sigma=0.015),C=70,cross=4,prob.model=TRUE)
gene

## predict gene type probabilities on the test set
genetype <- predict(gene,genetest,type="probabilities")
genetype
}

