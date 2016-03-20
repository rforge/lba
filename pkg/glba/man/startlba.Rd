\name{startlba}
\alias{startlba}

\title{
Specify and fit lba models. 
}
\description{

Specification of lba models with explanatory variables on the parameters. Models are fitted
using full information maximum likelihood. 

}
\usage{
startlba(rt, resp, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rt}{A vector of response times.}
  \item{resp}{A vector of (binary) responses.}
  \item{...}{Not currently used.}
}

\details{

The \code{startlba} function returns a named vector with starting values
for the (intercepts) of the parameter models.  Currently defined as: the
standard deviation of the RTs is used as the standard deviation of the
drift rate; the drift rate for the correct responses is initialized as the
mean accuracy; the non-decision parameter is set to half the minimum of the
response times; the start point parameter is set to drift*(qrt2-nond) where
drift is the drift rate, qrt2 the second quartile of the RT distribution,
and nond is the non-decision time; the bound is set to drift*(qrt4-nond)
where qrt4 is the fourth value of quartiles of the RT distribution.

}

\value{
The \code{startlba} function returns a named vector with starting values
for the (intercepts) of the parameter models.
}

\author{
Ingmar Visser
}

\examples{
data(bh08)
# remove extreme RTs
bh08 <- bh08[bh08$rt>.180&bh08$rt<2,]
startlba(bh08$rt,bh08$acc)
}

\keyword{statmodel}
