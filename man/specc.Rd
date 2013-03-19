\name{specc}
\alias{specc}
\alias{specc,matrix-method}
\alias{specc,formula-method}
\alias{specc,list-method}
\alias{specc,kernelMatrix-method}
\alias{show,specc-method}
\title{Spectral Clustering}
\description{
A spectral clustering algorithm. Clustering is performed by
embedding the data into the subspace of the eigenvectors 
of an affinity matrix.
}
\usage{
\S4method{specc}{formula}(x, data = NULL, na.action = na.omit, ...)

\S4method{specc}{matrix}(x, centers,
      kernel = "rbfdot", kpar = "automatic", 
      nystrom.red = FALSE, nystrom.sample = dim(x)[1]/6,
      iterations = 200, mod.sample = 0.75, na.action = na.omit, ...)

\S4method{specc}{kernelMatrix}(x, centers, nystrom.red = FALSE, iterations = 200, ...)

\S4method{specc}{list}(x, centers,
      kernel = "stringdot", kpar = list(length=4, lambda=0.5),
      nystrom.red = FALSE, nystrom.sample = length(x)/6,
      iterations = 200, mod.sample = 0.75, na.action = na.omit, ...)
}

\arguments{
 \item{x}{the matrix of data to be clustered, or a symbolic
    description of the model to be fit, or a kernel Matrix of class
    \code{kernelMatrix}, or a list of character vectors.}

  \item{data}{an optional data frame containing the variables in the model.
    By default the variables are taken from the environment which
    `specc' is called from.}

\item{centers}{Either the number of clusters or a set of initial cluster
    centers. If the first, a random set of rows in the eigenvectors
    matrix are chosen as the initial centers.}

\item{kernel}{the kernel function used in computing the affinity matrix. 
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
      \item \code{stringdot} String kernel
    }
    The kernel parameter can also be set to a user defined function of
    class kernel by passing the function name as an argument.
  }

  \item{kpar}{a character string or the list of hyper-parameters (kernel parameters).
    The default character string \code{"automatic"} uses a heuristic to determine a
    suitable value for the width parameter of the RBF kernel.
    The second option \code{"local"} (local scaling) uses a more advanced heuristic
     and sets a width parameter for every point in the data set. This is
    particularly useful when the data incorporates multiple scales.
    A list can also be used containing the parameters to be used with the
    kernel function. Valid parameters for existing kernels are :
    \itemize{
      \item \code{sigma} inverse kernel width for the Radial Basis
      kernel function "rbfdot" and the Laplacian kernel "laplacedot".
      \item \code{degree, scale, offset} for the Polynomial kernel "polydot"
      \item \code{scale, offset} for the Hyperbolic tangent kernel
      function "tanhdot"
      \item \code{sigma, order, degree} for the Bessel kernel "besseldot". 
      \item \code{sigma, degree} for the ANOVA kernel "anovadot".
      \item \code{length, lambda, normalized} for the "stringdot" kernel
       where length is the length of the strings considered, lambda the
       decay factor and normalized a logical parameter determining if the
       kernel evaluations should be normalized.
    }
    
    Hyper-parameters for user defined kernels can be passed through the
    kpar parameter as well.}

  \item{nystrom.red}{use nystrom method to calculate eigenvectors. When
    \code{TRUE} a sample of the dataset is used to calculate the
    eigenvalues, thus only a \eqn{n x m} matrix where \eqn{n} the sample size
    is stored in memory (default: \code{FALSE}}

  \item{nystrom.sample}{number of data points to use for estimating the
    eigenvalues when using the nystrom method. (default : dim(x)[1]/6)}

  \item{mod.sample}{proportion of data to use when estimating sigma (default: 0.75)}	

  \item{iterations}{the maximum number of iterations allowed. }

  \item{na.action}{the action to perform on NA}

  \item{\dots}{additional parameters}
    
}
\details{
  Spectral clustering  works by embedding the data points of the
  partitioning problem into the
  subspace of the \eqn{k} largest eigenvectors of a normalized affinity/kernel matrix.
Using a simple clustering method like \code{kmeans} on the embedded points usually
leads to good performance. It can be shown that spectral clustering methods boil down to 
 graph partitioning.\cr
The data can be passed to the \code{specc} function in a \code{matrix} or a
\code{data.frame}, in addition \code{specc} also supports input in the form of a
kernel matrix of class \code{kernelMatrix} or as a list of character
vectors where a string kernel has to be used.}
\value{
 An S4 object of class \code{specc} which extends the class \code{vector}
 containing integers indicating the cluster to which
 each point is allocated. The following slots contain useful information
 
  \item{centers}{A matrix of cluster centers.}
  \item{size}{The number of point in each cluster}
  \item{withinss}{The within-cluster sum of squares for each cluster}
  \item{kernelf}{The kernel function used}
}
\references{
  Andrew Y. Ng, Michael I. Jordan, Yair Weiss\cr
  \emph{On Spectral Clustering: Analysis and an Algorithm}\cr
  Neural Information Processing Symposium 2001\cr
  \url{http://www.nips.cc/NIPS2001/papers/psgz/AA35.ps.gz}

}
\author{Alexandros Karatzoglou \cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}
}


\seealso{\code{\link{kkmeans}}, \code{\link{kpca}}, \code{\link{kcca}} }
\examples{
## Cluster the spirals data set.
data(spirals)

sc <- specc(spirals, centers=2)

sc
centers(sc)
size(sc)
withinss(sc)

plot(spirals, col=sc)

}
\keyword{cluster}

