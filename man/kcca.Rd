\name{kcca}
\alias{kcca}
\alias{kcca,matrix-method}
\title{Kernel Canonical Correlation Analysis}
\description{
Computes the canonical correlation analysis in a feature space.
}
\usage{
\S4method{kcca}{matrix}(x, y, kernel="rbfdot", kpar=list(sigma=0.1), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a matrix containing data index by row}
  \item{y}{a matrix containing data index by row}
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

  \item{\dots}{adittional parameters for the \code{kpca} function}
}
\details{
The kernel version of canonical correlation analysis. 
}
\value{
  An S4 object containg the following slots:
  \item{kcor}{Correlation coefficients in feature space}
  \item{xcoef}{estimated coefficients for the \code{x} variables in the
    feature space}
  \item{ycoef}{estimated coefficients for the \code{y} variables in the
    feature space}
  \item{xvar}{The canonical variates for \code{x}}
  \item{yvar}{The canonical variates for \code{y}}
    
}
\references{ Malte Kuss, Thore Graepel \cr 
	\emph{The Geometry Of Kernel Canonical Correlation Analysis}\cr
	\url{http://www.kyb.tuebingen.mpg.de/publications/pdfs/pdf2233.pdf}}
\author{ Alexandros Karatzoglou \cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}}


\seealso{\code{\link{cancor}}  \code{\link{kpca}}}
\examples{

}
\keyword{multivariate}

