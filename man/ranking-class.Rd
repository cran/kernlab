\name{ranking-class}
\docType{class}
\alias{ranking-class}
\alias{edgegraph}
\alias{convergence}
\alias{convergence,ranking-method}
\alias{edgegraph,ranking-method}
\alias{show,ranking-method}

\title{Class "ranking"}
\description{Object of the class \code{"ranking"} are created from the
  \code{ranking} function and extend the class \code{matrix}}
\section{Objects from the Class}{

Objects can be created by calls of the form \code{new("ranking", ...)}.

}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"matrix"} containing the
      data ranking and scores}
    \item{\code{convergence}:}{Object of class \code{"matrix"}
      containing the convergence matrix}
    \item{\code{edgegraph}:}{Object of class \code{"matrix"} containing
      the edgegraph}
  }
}
\section{Extends}{
Class \code{"matrix"}, directly.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "ranking")}: displays the
      ranking score matrix}
  }
}

\author{Alexandros Karatzoglou \cr
\email{alexandros.karatzoglou@ci.tuwien.ac.at}
}

\seealso{
  \code{\link{ranking}}
}
\examples{
data(spirals)

## create data set to be ranked
ran<-spirals[rowSums(abs(spirals)<0.55)==2,]

## rank points according to "relevance" to point 54 (up left)
ranked<-ranking(ran,54,kernel="rbfdot",
                kpar=list(sigma=100),edgegraph=TRUE)

ranked
edgegraph(ranked)[1:10,1:10]
}
\keyword{classes}
