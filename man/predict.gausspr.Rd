\name{predict.gausspr}
\alias{predict.gausspr}
\alias{predict,gausspr-method}
\title{predict method for Gaussian Processes object}


\description{Prediction of test data using Gaussian Processes}


\usage{
\S4method{predict}{gausspr}(object, newdata, type = "response", coupler = "minpair")
}

\arguments{

  \item{object}{an S4 object of class \code{gausspr} created by the
    \code{gausspr} function}
  \item{newdata}{a data frame or matrix containing new data}
  \item{type}{one of \code{response}, \code{probabilities} 
    indicating the type of output: predicted values or matrix of class
    probabilities}
  \item{coupler}{Coupling method used in the multiclass case, can be one
    of \code{minpair} or \code{pkpd} (see reference for more details).}

}

\value{
      \item{response}{predicted classes (the classes with majority vote)
       or the response value in regression.}

     \item{probabilities}{matrix of class probabilities (one column for each class and
       one row for each input).}
   }


   \references{
 \itemize{

  \item
    C. K. I. Williams and D. Barber \cr
    Bayesian classification with Gaussian processes. \cr
    IEEE Transactions on Pattern Analysis and Machine Intelligence, 20(12):1342-1351, 1998\cr
    \url{http://www.dai.ed.ac.uk/homes/ckiw/postscript/pami_final.ps.gz}
  
   \item
       T.F. Wu, C.J. Lin, R.C. Weng. \cr
       \emph{Probability estimates for Multi-class Classification by
	 Pairwise Coupling}\cr
        \url{http://www.csie.ntu.edu.tw/~cjlin/papers/svmprob/svmprob.pdf}
 
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
gene <- gausspr(Class~.,data=genetrain,kernel="rbfdot",
                kpar=list(sigma=0.015))
gene

## predict gene type probabilities on the test set
genetype <- predict(gene,genetest,type="probabilities")
genetype
}

