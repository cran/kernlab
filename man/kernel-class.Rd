\name{kernel-class}
\docType{class}
\alias{rbfkernel-class}
\alias{polykernel-class}
\alias{vanillakernel-class}
\alias{tanhkernel-class}
\alias{anovakernel-class}
\alias{besselkernel-class}
\alias{laplacekernel-class}
\alias{splinekernel-class}
\alias{stringkernel-class}
\alias{fourierkernel-class}
\alias{kfunction-class}

\alias{kernel-class}
\alias{kpar,kernel-method}
\title{Class "kernel" "rbfkernel" "polykernel", "tanhkernel", "vanillakernel"}
\description{  The built-in kernel classes in \pkg{kernlab}}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("rbfkernel")},
\code{new{"polykernel"}}, \code{new{"tanhkernel"}},
\code{new{"vanillakernel"}}, \code{new{"anovakernel"}},
\code{new{"besselkernel"}}, \code{new{"laplacekernel"}},
\code{new{"splinekernel"}}, \code{new{"stringkernel"}}

or by calling the \code{rbfdot}, \code{polydot}, \code{tanhdot},
\code{vanilladot}, \code{anovadot}, \code{besseldot}, \code{laplacedot},
\code{splinedot}, \code{stringdot} functions etc..
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
	"matrix")}: computes the kernel matrix}
    \item{kernelMult}{\code{signature(kernel = "rbfkernel", x =
	"matrix")}: computes the quadratic kernel expression}
    \item{kernelPol}{\code{signature(kernel = "rbfkernel", x =
	"matrix")}: computes the kernel expansion}
    \item{kernelFast}{\code{signature(kernel = "rbfkernel", x =
	"matrix"),,a}: computes parts or the full kernel matrix, mainly
      used in kernel algorithms where columns of the kernel matrix are
      computed per invocation }
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
