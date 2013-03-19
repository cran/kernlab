\name{ranking}
\alias{ranking}
\alias{ranking,matrix-method}
\alias{ranking,list-method}
\alias{ranking,kernelMatrix-method}

\title{Ranking}
\description{
 A universal ranking algorithm which assigns importance/ranking to data points
 given a query.
}
\usage{
\S4method{ranking}{matrix}(x, y,
        kernel ="rbfdot", kpar = list(sigma = 1),
        scale = FALSE, alpha = 0.99, iterations = 600,
        edgegraph = FALSE, convergence = FALSE ,...)

\S4method{ranking}{kernelMatrix}(x, y,
        alpha = 0.99, iterations = 600, convergence = FALSE,...)

\S4method{ranking}{list}(x, y,
        kernel = "stringdot", kpar = list(length = 4, lambda = 0.5),
        alpha = 0.99, iterations = 600, convergence = FALSE, ...)

}

\arguments{
  \item{x}{a matrix containing the data to be ranked, or the kernel
    matrix of data to be ranked or a list of character vectors}
  \item{y}{The index of the query point in the data matrix or a vector
    of length equal to the rows of the data matrix having a one at the
    index of the query points index and zero at all the other points.}

 
    \item{kernel}{the kernel function used in training and predicting.
    This parameter can be set to any function, of class kernel, which computes a dot product between two
    vector arguments. kernlab provides the most popular kernel functions
    which can be used by setting the kernel parameter to the following
    strings:
    \itemize{
      \item \code{rbfdot} Radial Basis kernel function "Gaussian"
      \item \code{polydot} Polynomial kernel function
      \item \code{vanilladot} Linear kernel function
      \item \code{tanhdot} Hyperbolic tangent kernel function
      \item \code{laplacedot} Laplacian kernel function
      \item \code{besseldot} Bessel kernel function
      \item \code{anovadot} ANOVA RBF kernel function
      \item \code{splinedot} Spline kernel 
    }
    The kernel parameter can also be set to a user defined function of
    class kernel by passing the function name as an argument.
  }

  \item{kpar}{the list of hyper-parameters (kernel parameters).
    This is a list which contains the parameters to be used with the
    kernel function. For valid parameters for existing kernels are :
    \itemize{
      \item \code{sigma} inverse kernel width for the Radial Basis
      kernel function "rbfdot" and the Laplacian kernel "laplacedot".
      \item \code{degree, scale, offset} for the Polynomial kernel "polydot"
      \item \code{scale, offset} for the Hyperbolic tangent kernel
      function "tanhdot"
      \item \code{sigma, order, degree} for the Bessel kernel "besseldot".
      \item \code{sigma, degree} for the ANOVA kernel "anovadot".
    }
    Hyper-parameters for user defined kernels can be passed through the
    kpar parameter as well.}
  \item{scale}{If TRUE the data matrix columns are scaled to zero mean
    and unit variance.}
  \item{alpha}{ The \code{alpha} parameter takes values between 0 and 1
    and is used to control the authoritative scores received from the
    unlabeled points. For 0 no global structure is found the algorithm
    ranks the points similarly to the original distance metric.}
  \item{iterations}{Maximum number of iterations}
  \item{edgegraph}{Construct edgegraph (only supported with the RBF
    kernel)}
  \item{convergence}{Include convergence matrix in results}
  \item{\dots}{Additional arguments}
}
\details{
A simple universal ranking algorithm which exploits the intrinsic
global geometric structure of the data. In many real world
applications this should be superior to a local method in which the data
are simply ranked by pairwise Euclidean distances.
Firstly a weighted network is defined on the data and an authoritative
score is assigned to each query. The query points act as source nodes
that continually pump their authoritative scores to the remaining points
via the weighted network and the remaining points further spread the
scores they received to their neighbors. This spreading process is
repeated until convergence and the points are ranked according to their
score at the end of the iterations.
}
\value{
  An S4 object of class \code{ranking} which extends the \code{matrix}
  class.
  The first column of the returned matrix contains the original index of
  the points in the data matrix the second column contains the final
  score received by each point and the third column the ranking of the point.
  The object contains the following slots :
  \item{edgegraph}{Containing the edgegraph of the data points. }
  \item{convergence}{Containing the convergence matrix}
}
\references{
  D. Zhou, J. Weston, A. Gretton, O. Bousquet, B. Schoelkopf \cr
  \emph{Ranking on Data Manifolds}\cr
  Advances in Neural Information Processing Systems 16.\cr
  MIT Press Cambridge Mass. 2004 \cr
  \url{http://www.kyb.mpg.de/publications/pdfs/pdf2334.pdf}
}

\author{Alexandros Karatzoglou \cr
\email{alexandros.karatzoglou@ci.tuwien.ac.at}}

\seealso{ \code{\link{ranking-class}}, \code{\link{specc}} }
\examples{
data(spirals)

## create data from spirals
ran <- spirals[rowSums(abs(spirals) < 0.55) == 2,]

## rank points according to similarity to the most upper left point  
ranked <- ranking(ran, 54, kernel = "rbfdot",
                  kpar = list(sigma = 100), edgegraph = TRUE)
ranked[54, 2] <- max(ranked[-54, 2])
c<-1:86
op <- par(mfrow = c(1, 2),pty="s")
plot(ran)
plot(ran, cex=c[ranked[,3]]/40)

}
\keyword{cluster}
\keyword{classif}
