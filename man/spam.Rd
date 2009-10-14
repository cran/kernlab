\name{spam}
\alias{spam}
\title{Spam E-mail Database}
\description{A data set collected at Hewlett-Packard Labs, that classifies 4601
e-mails as spam or non-spam. In addition to this class label there are 57
variables indicating the frequency of certain words and characters in the
e-mail.}
\usage{data(spam)}
\format{A data frame with 4601 observations and 58 variables.

The first 48 variables contain the frequency of the variable name
(e.g., business) in the e-mail. If the variable name starts with num (e.g.,
num650) the it indicates the frequency of the corresponding number (e.g., 650).
The variables 49-54 indicate the frequency of the characters `;', `(', `[', `!',
`\$', and `\#'. The variables 55-57 contain the average, longest 
and total run-length of capital letters. Variable 58 indicates the type of the
mail and is either \code{"nonspam"} or \code{"spam"}, i.e. unsolicited
commercial e-mail.}

\details{
The data set contains 2788 e-mails classified as \code{"nonspam"} and 1813
classified as \code{"spam"}.

The ``spam'' concept is diverse: advertisements for products/web
sites, make money fast schemes, chain letters, pornography...
This collection of spam e-mails came from the collectors' postmaster and
individuals who had filed spam.  The collection of non-spam
e-mails came from filed work and personal e-mails, and hence
the word 'george' and the area code '650' are indicators of
non-spam.  These are useful when constructing a personalized
spam filter.  One would either have to blind such non-spam
indicators or get a very wide collection of non-spam to
generate a general purpose spam filter.
}
\source{
\itemize{
\item Creators: Mark Hopkins, Erik Reeber, George Forman, Jaap Suermondt at
Hewlett-Packard Labs, 1501 Page Mill Rd., Palo Alto, CA 94304
\item Donor: George Forman (gforman at nospam hpl.hp.com)  650-857-7835
}
These data have been taken from the UCI Repository Of Machine Learning
Databases at \url{http://www.ics.uci.edu/~mlearn/MLRepository.html}}
\references{
T. Hastie, R. Tibshirani, J.H. Friedman. \emph{The Elements of Statistical
Learning.} Springer, 2001.
}

\keyword{datasets}
