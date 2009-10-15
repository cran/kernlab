\name{income}
\alias{income}
\title{Income Data}
\description{
Customer Income Data from a marketing survey.
}
\usage{data(income)}

\format{
  A data frame with 14 categorical variables (8993 observations).

  Explanation of the variable names:

  \tabular{rllll}{
    \tab 1 \tab \code{INCOME} \tab annual income of household \tab \cr
    \tab   \tab               \tab (Personal income if single) \tab ordinal\cr
    \tab 2 \tab \code{SEX} \tab sex \tab nominal\cr
    \tab 3 \tab \code{MARITAL.STATUS} \tab marital status \tab nominal\cr
    \tab 4 \tab \code{AGE} \tab age \tab ordinal\cr
    \tab 5 \tab \code{EDUCATION} \tab educational grade \tab ordinal\cr
    \tab 6 \tab \code{OCCUPATION} \tab type of work \tab nominal \cr
    \tab 7 \tab \code{AREA} \tab how long the interviewed person has lived\tab
    \cr
    \tab   \tab             \tab in the San Francisco/Oakland/San Jose area \tab ordinal\cr
    \tab 8 \tab \code{DUAL.INCOMES} \tab dual incomes (if married) \tab nominal\cr
    \tab 9 \tab \code{HOUSEHOLD.SIZE} \tab persons living in the
    household \tab ordinal\cr
    \tab 10 \tab \code{UNDER18} \tab persons in household under 18 \tab ordinal\cr
    \tab 11 \tab \code{HOUSEHOLDER} \tab householder status \tab nominal\cr
    \tab 12 \tab \code{HOME.TYPE} \tab type of home \tab nominal\cr
    \tab 13 \tab \code{ETHNIC.CLASS} \tab ethnic classification \tab nominal\cr
    \tab 14 \tab \code{LANGUAGE} \tab language most often spoken at
    home \tab nominal\cr
  }
}
\details{
A total of N=9409 questionnaires containing 502 questions were 
filled out by shopping mall customers in the San Francisco Bay area. 
The dataset is an extract from this survey. It consists of 
14 demographic attributes. The dataset is a mixture of nominal and
ordinal variables with a lot of missing data.
The goal is to predict the Anual Income of Household from the other 13
demographics attributes.
}
\source{
  Impact Resources, Inc., Columbus, OH (1987).
}
\keyword{datasets}
