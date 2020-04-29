\name{rlba}

\alias{rlba}
\alias{rlba2}

\title{Generate data from an LBA model.}

\description{Generate data from an LBA model.}

\usage{
rlba(n, b, A, vs, s, t0, st0 = 0, truncdrifts = TRUE)

rlba2(nn, bnd, sp, drift, sddr, nond, st0 = 0)
}
\arguments{
  \item{n}{The number of random variates required.}
  \item{b}{The boundary parameter.}
  \item{A}{The start point parameter.}
  \item{vs}{The drift rates for each of the accumulators.}
  \item{s}{The standard deviation of the drift rates.}
  \item{t0}{The non-decision time.}
  \item{st0}{Variability of the non-decision time, uniformly distributed as t0+U(0,st0).}
  \item{truncdrifts}{Logical; whether drifts that are generated should be
  truncated at zero.  Negative RTs may result when set to FALSE.}
  
  \item{nn}{The number of random variates required.}
  \item{bnd}{The boundary parameter.}
  \item{sp}{The start point parameter.}
  \item{drift}{The drift rates for each of the accumulators.}
  \item{sddr}{The standard deviation of the drift rates.}
  \item{nond}{The non-decision time.}
  
}

\details{
The \code{rlba} function is the one used by Brown and Heathcote; the \code{rlba2} 
function is the fully vectorized version of the same function; both provide random variates
under the lba distribution; \code{rlba2} allows all the arguments to be vectors of length 'n' 
as required to enable simulation in the \code{glba}-package where potentially all parameters
depend on covariates and are hence different for each datapoint.
}

\value{A data.frame with two variables, 'rt' and 'resp', for the response
time and the response, respectively.}

\references{
Function adapted from Brown and Heathcote, code extracted in 2015 from 
this website: http://www.newcl.org/publications/2008/lbasoftware.htm
}

\author{
Ingmar Visser, function adapted from Brown and Heathcote, code extracted in 2015 
from this website: http://www.newcl.org/publications/2008/lbasoftware.htm
}

\examples{
# to be added later
}
\keyword{statmodel}
