\name{specc}
\alias{specc}
\alias{specc,matrix-method}
\alias{specc,formula-method}

\title{Spectral Clustering}
\description{
A spectral clustering  algorithm. This algorithm clusters points using
eigenvectors of kernel matrixes derived from the data.  
}
\usage{
\S4method{specc}{formula}(x, data = NULL, na.action = na.omit, ...)

\S4method{specc}{matrix}(x, centers, kernel = "rbfdot", kpar =
list(sigma = 0.1), iterations = 200, na.action = na.omit, ...)
}

\arguments{
  \item{x}{the matrix of data to be clustered or a symbolic
    description of the model to be fit.}
  \item{data}{an optional data frame containing the variables in the model.
    By default the variables are taken from the environment which
    `specc' is called from.}
  \item{centers}{Either the number of clusters or a set of initial cluster
    centers. If the first, a random set of rows in the eigenvectors
    matrix are chosen as the initial centers.}
	
	\item{kernel}{the kernel function used in training and predicting.
	  This parameter can be set to any function, of class kernel, which computes a dot product between two
	  vector arguments. kernlab provides the most popular kernel functions
	  which can be used by setting the kernel parameter to the following
	  strings:
	  \itemize{
      \item \code{rbfdot} (Radial Basis kernel function)
      \item \code{polydot} (Polynomial kernel function)
      \item \code{vanilladot} (Linear kernel function)
      \item \code{tanhdot} (Hyperbolic tangent kernel function)
    }
    The kernel parameter can also be set to a user defined function of
    class kernel by passing the function name as an argument.
  }

  \item{kpar}{the list of hyper-parameters (kernel parameters).
    This is a list which contains the parameters to be used with the
    kernel function. For valid parameters for existing kernels are :
    \itemize{
      \item \code{sigma} (inverse kernel width for the Radial Basis kernel function "rbfdot")
      \item \code{degree, scale, offset} (for the Polynomial kernel "polydot")
      \item \code{scale, offset} (for the Hyperbolic tangent kernel
      function "tanhdot")
    }
    Hyper-parameters for user defined kernels can be passed through the
    kpar parameter as well.}

  \item{iterations}{The maximum number of iterations allowed. }

  \item{na.action}{ The action to perform on NA}

  \item{\dots}{additional parameters}
    
}
\details{
 In Spectral Clustering one uses the top \code{k} (number of clusters) eigenvectors of a matrix derived
 from the distance between points. Very good results are obtained by
 using a standard clustering technique
 to cluster the resulting eigenvector matrixes. 
}
\value{
 An S4 object of class \code{specc}
  \item{cluster}{A vector of integers indicating the cluster to which
   each point is allocated.}
  \item{centers}{A matrix of cluster centers.}
 \item{size}{The number of point in each cluster}
}
\references{
  Andrew Y. Ng, Michael I. Jordan, Yair Weiss\cr
  \emph{On Spectral Clustering: Analysis and an Algorithm}\cr
  Neural Information Processing Symposium 2001\cr
  \url{http://www.nips.cc/NIPS2001/papers/psgz/AA35.ps.gz}

}
\author{ Alexandros Karatzoglou \cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}
}


\seealso{ \code{\link{kpca}}, \code{\link{kcca}} }
\examples{
## Cluster the spirals data set.
data(spirals)

sc <- specc(spirals, centers=2, kernel="rbfdot", kpar=list(sigma=400))

centers(sc)
size(sc)

plot(spirals, col=cluster(sc))

}
\keyword{cluster}

