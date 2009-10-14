\name{ticdata}
\alias{ticdata}
\title{The Insurance Company Data}
\description{
This data set used in the CoIL 2000 Challenge contains information on customers of an insurance
company. The data consists of 86 variables and includes product usage data and socio-demographic
data derived from zip area codes.
The data was collected to answer the following question: Can you
predict who would be interested in buying a caravan insurance policy and give an explanation why ?
}
\usage{data(ticdata)}
\format{

ticdata: Dataset to train and validate prediction models and build a description (9822
customer records). Each record consists of 86 attributes, containing
sociodemographic data (attribute 1-43) and product ownership (attributes 44-86).
The sociodemographic data is derived from zip codes. All customers
living in areas with the same zip code have the same sociodemographic
attributes. Attribute 86, \code{CARAVAN:Number of mobile home policies}, is the target variable.

Data Format

  
\tabular{rlll}{
\tab 1 \tab \code{STYPE} \tab Customer Subtype\cr
\tab 2 \tab \code{MAANTHUI} \tab Number of houses 1 - 10\cr
\tab 3 \tab \code{MGEMOMV} \tab Avg size household 1 - 6\cr
\tab 4 \tab \code{MGEMLEEF} \tab Average age\cr
\tab 5 \tab \code{MOSHOOFD} \tab Customer main type\cr
\tab 6 \tab \code{MGODRK} \tab Roman catholic   \cr
\tab 7 \tab \code{MGODPR} \tab Protestant ...   \cr
\tab 8 \tab \code{MGODOV} \tab Other religion    \cr
\tab 9 \tab \code{MGODGE} \tab No religion      \cr
\tab 10 \tab \code{MRELGE} \tab Married        \cr
\tab 11 \tab \code{MRELSA} \tab Living together  \cr
\tab 12 \tab \code{MRELOV} \tab Other relation   \cr
\tab 13 \tab \code{MFALLEEN} \tab Singles         \cr
\tab 14 \tab \code{MFGEKIND} \tab Household without children    \cr
\tab 15 \tab \code{MFWEKIND} \tab Household with children \cr
\tab 16 \tab \code{MOPLHOOG} \tab High level education    \cr
\tab 17 \tab \code{MOPLMIDD} \tab Medium level education  \cr
\tab 18 \tab \code{MOPLLAAG} \tab Lower level education  \cr
\tab 19 \tab \code{MBERHOOG} \tab High status      \cr
\tab 20 \tab \code{MBERZELF} \tab Entrepreneur       \cr
\tab 21 \tab \code{MBERBOER} \tab Farmer            \cr
\tab 22 \tab \code{MBERMIDD} \tab Middle management  \cr
\tab 23 \tab \code{MBERARBG} \tab Skilled labourers  \cr
\tab 24 \tab \code{MBERARBO} \tab Unskilled labourers \cr
\tab 25 \tab \code{MSKA} \tab Social class A          \cr
\tab 26 \tab \code{MSKB1} \tab Social class B1        \cr
\tab 27 \tab \code{MSKB2} \tab Social class B2        \cr
\tab 28 \tab \code{MSKC} \tab Social class C           \cr
\tab 29 \tab \code{MSKD} \tab Social class D           \cr
\tab 30 \tab \code{MHHUUR} \tab Rented house        \cr
\tab 31 \tab \code{MHKOOP} \tab Home owners         \cr
\tab 32 \tab \code{MAUT1} \tab 1 car                \cr
\tab 33 \tab \code{MAUT2} \tab 2 cars                \cr
\tab 34 \tab \code{MAUT0} \tab No car   \cr
\tab 35 \tab \code{MZFONDS} \tab National Health Service    \cr
\tab 36 \tab \code{MZPART} \tab Private health insurance     \cr
\tab 37 \tab \code{MINKM30} \tab Income  >30.000      \cr
\tab 38 \tab \code{MINK3045} \tab Income 30-45.000    \cr
\tab 39 \tab \code{MINK4575} \tab Income 45-75.000    \cr
\tab 40 \tab \code{MINK7512} \tab Income 75-122.000   \cr
\tab 41 \tab \code{MINK123M} \tab Income <123.000    \cr
\tab 42 \tab \code{MINKGEM} \tab Average income       \cr
\tab 43 \tab \code{MKOOPKLA} \tab Purchasing power class  \cr
\tab 44 \tab \code{PWAPART} \tab Contribution private third party insurance  \cr
\tab 45 \tab \code{PWABEDR} \tab Contribution third party insurance (firms)  \cr
\tab 46 \tab \code{PWALAND} \tab Contribution third party insurance (agriculture) \cr
\tab 47 \tab \code{PPERSAUT} \tab Contribution car policies          \cr
\tab 48 \tab \code{PBESAUT} \tab Contribution delivery van policies  \cr
\tab 49 \tab \code{PMOTSCO} \tab Contribution motorcycle/scooter policies   \cr
\tab 50 \tab \code{PVRAAUT} \tab Contribution lorry policies             \cr
\tab 51 \tab \code{PAANHANG} \tab Contribution trailer policies         \cr
\tab 52 \tab \code{PTRACTOR} \tab Contribution tractor policies        \cr
\tab 53 \tab \code{PWERKT} \tab Contribution agricultural machines policies   \cr
\tab 54 \tab \code{PBROM} \tab Contribution moped policies       \cr
\tab 55 \tab \code{PLEVEN} \tab Contribution life insurances     \cr
\tab 56 \tab \code{PPERSONG} \tab Contribution private accident insurance policies \cr
\tab 57 \tab \code{PGEZONG} \tab Contribution family accidents insurance policies \cr
\tab 58 \tab \code{PWAOREG} \tab Contribution disability insurance policies \cr
\tab 59 \tab \code{PBRAND} \tab Contribution fire policies           \cr
\tab 60 \tab \code{PZEILPL} \tab Contribution surfboard policies    \cr
\tab 61 \tab \code{PPLEZIER} \tab Contribution boat policies      \cr
\tab 62 \tab \code{PFIETS} \tab Contribution bicycle policies       \cr
\tab 63 \tab \code{PINBOED} \tab Contribution property insurance policies   \cr
\tab 64 \tab \code{PBYSTAND} \tab Contribution social security insurance policies  \cr
\tab 65 \tab \code{AWAPART} \tab Number of private third party insurance 1 - 12   \cr
\tab 66 \tab \code{AWABEDR} \tab Number of third party insurance (firms) ...   \cr
\tab 67 \tab \code{AWALAND} \tab Number of third party insurance (agriculture)  \cr
\tab 68 \tab \code{APERSAUT} \tab Number of car policies          \cr
\tab 69 \tab \code{ABESAUT} \tab Number of delivery van policies   \cr
\tab 70 \tab \code{AMOTSCO} \tab Number of motorcycle/scooter policies  \cr
\tab 71 \tab \code{AVRAAUT} \tab Number of lorry policies    \cr
\tab 72 \tab \code{AAANHANG} \tab Number of trailer policies  \cr
\tab 73 \tab \code{ATRACTOR} \tab Number of tractor policies  \cr
\tab 74 \tab \code{AWERKT} \tab Number of agricultural machines policies     \cr
\tab 75 \tab \code{ABROM} \tab Number of moped policies          \cr
\tab 76 \tab \code{ALEVEN} \tab Number of life insurances        \cr
\tab 77 \tab \code{APERSONG} \tab Number of private accident insurance policies   \cr
\tab 78 \tab \code{AGEZONG} \tab Number of family accidents insurance policies   \cr
\tab 79 \tab \code{AWAOREG} \tab Number of disability insurance policies         \cr
\tab 80 \tab \code{ABRAND} \tab Number of fire policies            \cr
\tab 81 \tab \code{AZEILPL} \tab Number of surfboard policies     \cr
\tab 82 \tab \code{APLEZIER} \tab Number of boat policies        \cr
\tab 83 \tab \code{AFIETS} \tab Number of bicycle policies         \cr
\tab 84 \tab \code{AINBOED} \tab Number of property insurance policies   \cr
\tab 85 \tab \code{ABYSTAND} \tab Number of social security insurance policies  \cr
\tab 86 \tab \code{CARAVAN} \tab Number of mobile home policies 0 - 1          \cr
}
Note: All the variables starting with M are zipcode variables. They give information on the
distribution of that variable, e.g., Rented house, in the zipcode area of
the customer.
}



\details{
  
  Information about the insurance company customers consists of 86
  variables and includes
  product usage data and socio-demographic data derived from zip area codes. The
  data was supplied by the Dutch data mining company Sentient Machine
  Research and is based on a real world business problem. The training
  set contains over 5000 descriptions of customers, including the
  information of whether or not they have a caravan insurance policy.
  The test set contains 4000 customers. The test and data set are merged in the
  ticdata set.
  More information about the data set and the CoIL 2000 Challenge along
  with publications based on the data set can be found at \url{http://www.liacs.nl/~putten/library/cc2000/}.
}
\source{
  \itemize{
    \item UCI KDD Archive:\url{http://kdd.ics.uci.edu} 
      \item Donor:  Sentient Machine Research \cr
    Peter van der Putten \cr
    Sentient Machine Research \cr
    Baarsjesweg 224 \cr
    1058 AA Amsterdam \cr
    The Netherlands \cr
    +31 20 6186927 \cr
    pvdputten@hotmail.com, putten@liacs.nl
  }
}
  \references{Peter van der Putten, Michel de Ruiter, Maarten van
    Someren  \emph{CoIL Challenge 2000 Tasks and Results: Predicting and
      Explaining Caravan Policy Ownership}\cr
    \url{http://www.liacs.nl/~putten/library/cc2000/}}
  

 



\keyword{datasets}
