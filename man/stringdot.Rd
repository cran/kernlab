\name{stringdot}
\alias{stringdot}
\title{String Kernel Functions}
\description{
  String kernels.
}
\usage{
stringdot(length = 4, lambda = 1.1, type = "spectrum", normalized = TRUE)
}

\arguments{

  \item{length}{The length of the substrings considered}

  \item{lambda}{The decay factor}

  \item{type}{Type of string kernel, currently the following kernels are
    supported : \cr
    
    \code{spectrum} the kernel considers only matching substring of
    exactly length \eqn{n} (also know as string kernel). Each such matching
    substring is given a constant weight. The length parameter in this
    kernel has to be \eqn{length > 1}.\cr

    \code{boundrange}
    this kernel (also known as boundrange) considers only matching substrings of length less than or equal to a
    given number N. This type of string kernel requires a length
    parameter \eqn{length > 1}\cr
        
    \code{constant}
    The kernel considers all matching substrings and assigns constant weight (e.g. 1) to each
    of them. This \code{constant} kernel does not require any additional
    parameter.\cr


    \code{exponential}
    Exponential Decay kernel where the substring weight decays as the
    matching substring gets longer. The kernel requires a decay factor \eqn{
      \lambda > 1}\cr

    \code{string} essentially identical to the spectrum kernel, only
    computed using a more conventional way.\cr

    \code{fullstring} essentially identical to the boundrange kernel
    only computed in a more conventional way. \cr
 }
  \item{normalized}{normalize string kernel values, (default: \code{TRUE})}
}
\details{
  The kernel generating functions are used to initialize a kernel function
  which calculates the dot (inner) product between two feature vectors in a
  Hilbert Space. These functions or their function generating names
  can be passed as a \code{kernel} argument on almost all
  functions in \pkg{kernlab}(e.g., \code{ksvm}, \code{kpca}  etc.).

  The string kernels calculate similarities between two strings
  (e.g. texts or sequences) by matching the common substring
  in the strings.  Different types of string kernel exists and are
  mainly distinguished by how the matching is performed i.e. some string
  kernels count the exact matchings of \eqn{n} characters (spectrum
  kernel) between the strings, others allow gaps (mismatch kernel) etc.
 
   
  }
\value{
 Returns an S4 object of class \code{stringkernel} which extents the
 \code{function} class. The resulting function implements the given
 kernel calculating the inner (dot) product between two character vectors.
 \item{kpar}{a list containing the kernel parameters (hyperparameters)
   used.}
 The kernel parameters can be accessed by the \code{kpar} function.
 }

\author{Alexandros Karatzoglou\cr
  \email{alexandros.karatzoglou@ci.tuwien.ac.at}}

\note{ The \code{spectrum} and \code{boundrange} kernel are faster and
  more efficient implementations of the \code{string} and
  \code{fullstring} kernels
  which will be still included in \code{kernlab} for the next two versions.
  
}




\seealso{ \code{\link{dots} }, \code{\link{kernelMatrix} }, \code{\link{kernelMult}}, \code{\link{kernelPol}}}
\examples{

sk <- stringdot(type="string", length=5)

sk



}
\keyword{symbolmath}

