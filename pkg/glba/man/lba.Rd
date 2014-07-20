\name{lba}
\alias{lba}
\title{
Specify and fit lba models. 
}
\description{

Specification of lba models with explanatory variables on the parameters. Models are fitted
using full information maximum likelihood. 

}
\usage{
lba(rt, response, data, weights, 
sddr = ~1, 
sp = ~1, 
bound = ~1, 
nond = ~1, 
drift = ~1, 
scaling = c("sum", "ratio", "constant-error"), 
loglink = c(FALSE, FALSE, FALSE, FALSE), 
startpars, 
fixed = NULL,
method = "Nelder-Mead",
hessian = FALSE,
lower = -Inf,
upper = Inf) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rt}{The name of the response time variable in data.frame 'data'.}
  
  \item{response}{The name of the response time variable in data.frame 'data'.}
  
  \item{data}{Data.frame containing the variables.}
  
  \item{weights}{Optional case weights.}

  \item{sddr}{The model for the standard deviation of the drift rate.}
  \item{sp}{The model for the start point parameter.}
  
  
  \item{bound}{The model for the boundary parameter.}
  
  \item{nond}{The model for the non-decision time.}
  
  \item{drift}{The model for the drift rate parameter.}
  
  \item{scaling}{Options for making the model identifiable; currently only
  the "sum" option is implemented (and hence the argument is ignored).}
  
  \item{loglink}{Wheter parameters (sddr, sp, nond, and bound) should be
  fitted using a log link function (to avoid boundary estimation issues).}
  
  \item{startpars}{Vector of starting values for the parameters.}
  
  \item{fixed}{Logical; when startpars is given 'fixed' can be used to
  specify whether parameters should remain fixed during estimation;
  sometimes required for more stable estimation.}
  
  \item{method}{Optimization method; passed to optim}
  
  \item{hessian}{Logical; indicating whether the estimated hessian at the minimum; passed to optim.}
  
  \item{lower, upper}{lower and upper boundaries for the parameters; only applicable with optim'ization method='L-BFGS-B'}  
  
}

\details{
This is the main function to specify lba models with expanatory variables.
}

\value{
A list with the following components:
	\item{pars }{The values of fitted parameters.}
}

\references{
LBA models as in Brown & Heathcote 2008.
}

\author{
Ingmar Visser
}

\examples{
data(bh08)
# remove extreme RTs
bh08 <- bh08[bh08$rt>.180&bh08$rt<2,]

# set some starting values
pst <- c(0.242, .321, .392, .233, 0.984, 0.791, 0.618)
# these are based on the quantile parameter estimates from B&H 2008:
# parsBH <- c(0.289, 0.321, 0.233, 0.392, 0.984, 0.791, 0.619)
# names(parsBH) <- c("s","A","ter","b","v1","v2","v3")

# use random starting value or those provided above
set.seed(2)
# fit the model
m1 <- lba(rt=rt,response=acc,drift=~diff-1,data=bh08,startpars=runif(7)) 

\dontrun{
# get the parameters
pp <- m1$pars
# simulate some data from the model for the three conditions
dt1 <- rlba(1000, b=pp[3]+pp[2], A=pp[2], vs=c(pp[5],1-pp[5]), s=pp[1], t0=pp[4], truncdrifts=TRUE)
dt2 <- rlba(1000, b=pp[3]+pp[2], A=pp[2], vs=c(pp[6],1-pp[6]), s=pp[1], t0=pp[4], truncdrifts=TRUE)
dt3 <- rlba(1000, b=pp[3]+pp[2], A=pp[2], vs=c(pp[7],1-pp[7]), s=pp[1], t0=pp[4], truncdrifts=TRUE)
dt <- rbind(dt1,dt2,dt3)
dt$diff <- gl(3,1000)
dt$resp <- ifelse(dt$resp==2,0,1)

# remove extreme RTs
dt <- dt[dt$rt>.180&dt$rt<2,]

# refit the model with the simualted data
m1rf <- lba(rt=rt,response=resp,drift=~diff-1,data=dt,startpars=pst)
}

}

\keyword{statmodel}
