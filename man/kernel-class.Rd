\name{kernel-class}
\docType{class}
\alias{rbfkernel-class}
\alias{polykernel-class}
\alias{vanillakernel-class}
\alias{tanhkernel-class}
\alias{kernel-class}
\alias{kpar,kernel-method}
\title{Class "kernel" "rbfkernel" "polykernel", "tanhkernel", "vanillakernel"}
\description{  The build in \code{kernel} classes in \code{kernelab}}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rbfkernel", ...)},
\code{new{"polykernel"}}, \code{new{"tanhkernel"}},
\code{new{"vanillakernel"}}
or by calling the \code{rbfdot}, \code{polydot}, \code{tanhdot},
\code{vanilladot} functions.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"function"} containing
      the kernel function }
    \item{\code{kpar}:}{Object of class \code{"list"} containing the
      kernel parameters }
  }
}
\section{Extends}{
Class \code{"kernel"}, directly.
Class \code{"function"}, by class \code{"kernel"}.
}
\section{Methods}{
  \describe{
    \item{kernelMatrix}{\code{signature(kernel = "rbfkernel", x =
	"matrix")}: computes the kernel matrix }
    \item{kernelMult}{\code{signature(kernel = "rbfkernel", x =
	"matrix")}: computes the quadratic kernel expression}
    \item{kernelPol}{\code{signature(kernel = "rbfkernel", x =
	"matrix")}: computes the kernel expansion }
  }
}

\author{Alexandros Karatzoglou\cr \email{alexandros.karatzoglou@ci.tuwien.ac.at} }



\seealso{
  \code{\link{dots}}
  
}
\examples{

rbfkernel <- rbfdot(sigma = 0.1)
rbfkernel
is(rbfkernel)
kpar(rbfkernel)

}
\keyword{classes}
