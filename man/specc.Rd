\name{specc}
\alias{specc}
\alias{specc,matrix-method}
\alias{specc,formula-method}
\alias{show,specc-method}
\title{Spectral Clustering}
\description{
A spectral clustering  algorithm. This algorithm clusters points using
eigenvectors of kernel matrixes derived from the data.  
}
\usage{
\S4method{specc}{formula}(x, data = NULL, na.action = na.omit, ...)

\S4method{specc}{matrix}(x, centers, kernel = "rbfdot", kpar = list(sigma = 0.1), 
       iterations = 200, mod.sample = 0.6, na.action = na.omit, ...)
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
	
 \item{mod.sample} {Proportion of data to use when estimating sigma default 0.6}	

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
 An S4 object of class \code{specc} wich extends the class \code{vector}
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
\author{ Alexandros Karatzoglou \cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}
}


\seealso{ \code{\link{kpca}}, \code{\link{kcca}} }
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

