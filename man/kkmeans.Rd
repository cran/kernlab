\name{kkmeans}
\alias{kkmeans}
\alias{kkmeans,matrix-method}
\alias{kkmeans,formula-method}
\title{Kernel k-means}
\description{
  A weigthed kernel version of the famous k-means algorithm.
}
\usage{
\S4method{kkmeans}{formula}(x, data = NULL, na.action = na.omit, ...)

\S4method{kkmeans}{matrix}(x, centers, kernel = "rbfdot", kpar = list(sigma = 0.1),alg="kkmeans", p=1, na.action = na.omit, ...)
}

\arguments{
  \item{x}{the matrix of data to be clustered or a symbolic
    description of the model to be fit.}
  \item{data}{an optional data frame containing the variables in the model.
    By default the variables are taken from the environment which
    `kkmeans' is called from.}
  \item{centers}{Either the number of clusters or a set of initial cluster
    centers. If the first, a random set of rows in the eigenvectors
    matrix are chosen as the initial centers.}


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

  \item{kpar}

  {a character string or the list of hyper-parameters (kernel parameters).
    The default character string \code{"automatic"} uses a heuristic the determine a
    suitable value for the width parameter of the RBF kernel.

    A list can also be used containing the parameters to be used with the
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
  \item{alg}{the algorithm to use. Options currently include
    \code{kkmeans} and \code{kerninghan}. }
  
  \item{p}{a parameter used to keep the affinity matrix positive semidefinite}

  \item{na.action}{ The action to perform on NA}

  \item{\dots}{additional parameters}
    
}
\details{
  The algorithm is implemented using the triangle inequality to avoid
  unnecessary and computational expensive distance calculations.
  This leads to significant speedup particularly on large data sets with
  a high number of clusters.
  With a particular choice of weights this algorithm becomes
  equivalent to Kernighan-Lin, and the norm-cut graph partitioning
  algorithms.
}
\value{
 An S4 object of class \code{specc} wich extends the class \code{vector}
 containing integers indicating the cluster to which
 each point is allocated. The following slots contain useful information
 
  \item{centers}{A matrix of cluster centers.}
  \item{size}{The number of point in each cluster}
  \item{withinss}{The within-cluster sum of squares for each cluster}
  \item{kernelf}{The kernel function used}
}
\references{
  Inderjit Dhillon, Yuqiang Guan, Brian Kulis\cr
  A Unified view of Kernel k-means, Spectral Clustering and Graph
  Partitioning\cr
  UTCS Technical Report\cr
  \url{http://www.cs.utexas.edu/users/kulis/pubs/spectral_techreport.pdf}
  
}
\author{ Alexandros Karatzoglou \cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}
}

\seealso{\code{\link{specc}}, \code{\link{kpca}}, \code{\link{kcca}} }
\examples{
## Cluster the iris data set.
data(iris)

sc <- kkmeans(as.matrix(iris[,-5]), centers=3)

sc
centers(sc)
size(sc)
withinss(sc)


}
\keyword{cluster}

