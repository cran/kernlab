\name{kmmd}
\alias{kmmd}
\alias{kmmd,matrix-method}
\alias{kmmd,list-method}
\alias{kmmd,kernelMatrix-method}
\alias{show,kmmd-method}
\alias{H0}
\alias{Asymbound}
\alias{Radbound}
\alias{mmdstats}
\alias{AsympH0}

\title{Kernel Maximum Mean Discrepancy.}
\description{The Kernel Maximum Mean Discrepancy \code{kmmd} performs
 a non-parametric distribution test.}
\usage{

\S4method{kmmd}{matrix}(x, y, kernel="rbfdot",kpar="automatic", alpha = 0.05,
     asymptotic = FALSE, replace = TRUE, ntimes = 150, frac = 1, ...)

\S4method{kmmd}{kernelMatrix}(x, y, Kxy, alpha = 0.05,
     asymptotic = FALSE, replace = TRUE, ntimes = 100, frac = 1, ...)

\S4method{kmmd}{list}(x, y, kernel="stringdot", 
     kpar = list(type = "spectrum", length = 4), alpha = 0.05,
     asymptotic = FALSE, replace = TRUE, ntimes = 150, frac = 1, ...)

}

\arguments{
  \item{x}{data values, in a \code{matrix},
    \code{list}, or \code{kernelMatrix}}
  
  \item{y}{data values, in a \code{matrix},
    \code{list}, or \code{kernelMatrix}}
	  
  \item{Kxy}{\code{kernlMatrix} between \eqn{x} and \eqn{y} values (only for the
    kernelMatrix interface)}

 \item{kernel}{the kernel function used in training and predicting.
   This parameter can be set to any function, of class kernel, which computes a dot product between two
   vector arguments. \code{kernlab} provides the most popular kernel functions
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
 
 \item{kpar}{the list of hyper-parameters (kernel parameters).
   This is a list which contains the parameters to be used with the
   kernel function. Valid parameters for existing kernels are :
   \itemize{
     \item \code{sigma} inverse kernel width for the Radial Basis
     kernel function "rbfdot" and the Laplacian kernel "laplacedot".
     \item \code{degree, scale, offset} for the Polynomial kernel "polydot"
     \item \code{scale, offset} for the Hyperbolic tangent kernel
     function "tanhdot"
     \item \code{sigma, order, degree} for the Bessel kernel "besseldot". 
     \item \code{sigma, degree} for the ANOVA kernel "anovadot".
     \item \code{lenght, lambda, normalized} for the "stringdot" kernel
     where length is the length of the strings considered, lambda the
     decay factor and normalized a logical parameter determining if the
     kernel evaluations should be normalized.
   }
   
   Hyper-parameters for user defined kernels can be passed
   through the \code{kpar} parameter as well. In the case of a Radial
   Basis kernel function (Gaussian) kpar can also be set to the
   string "automatic" which uses the heuristics in  'sigest' to
   calculate a good 'sigma' value for the Gaussian RBF or
   Laplace kernel, from the data. (default = "automatic").
 }
 
 \item{alpha}{the confidence level of the test (default: 0.05)}
 
 \item{asymptotic}{calculate the bounds asymptotically (suitable for
   smaller datasets) (default: FALSE)}
	
 \item{replace}{use replace when sampling for computing the asymptotic
   bounds (default : TRUE)}
 
 \item{ntimes}{number of times repeating the sampling procedure (default
   : 150)}
 
 \item{frac}{fraction of points to sample (frac : 1) }

 \item{\dots}{additional parameters.}
}
 
\details{\code{kmmd} calculates the kernel maximum mean discrepancy for
  samples from two distributions and conducts a test as to whether the samples are
  from different distributions with level \code{alpha}.
  
}
\value{
  An S4 object of class \code{kmmd} containing the
  results of whether the H0 hypothesis is rejected or not. H0 being
  that the samples \eqn{x} and \eqn{y} come from the same distribution.
 The object contains the following slots :
 \item{\code{H0}}{is H0 rejected (logical)}
 \item{\code{AsympH0}}{is H0 rejected according to the asymptotic bound (logical)}
 \item{\code{kernelf}}{the kernel function used.}
 \item{\code{mmdstats}}{the test statistics (vector of two)}
 \item{\code{Radbound}}{the Rademacher bound}
 \item{\code{Asymbound}}{the asymptotic bound}
 
  see \code{kmmd-class} for more details.
}

  \references{Gretton, A., K. Borgwardt, M. Rasch, B. Schoelkopf and A. Smola\cr
    \emph{A Kernel Method for the Two-Sample-Problem}\cr
    Neural Information Processing Systems 2006, Vancouver   \cr
     \url{http://www.kyb.mpg.de/publications/attachments/mmd_final_4193[1].pdf}
   }
   
  \author{Alexandros Karatzoglou \cr \email{alexandros.karatzoglou@ci.tuwien.ac.at}}
  
\seealso{\code{ksvm}}

\examples{
# create data
x <- matrix(runif(300),100)
y <- matrix(runif(300)+1,100)


mmdo <- kmmd(x, y)

mmdo
}

\keyword{htest}
\keyword{nonlinear}
\keyword{nonparametric}
