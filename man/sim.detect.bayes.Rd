% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayes.R
\name{sim.detect.bayes}
\alias{sim.detect.bayes}
\title{Simulation using Bayesian Stopping Rule}
\usage{
sim.detect.bayes(GEN0, GEN1, th0, detect.bayes, nulower = NULL,
  nuupper = NULL, ...)
}
\arguments{
\item{GEN0}{A function that take no argument and return an observation from
the pre-change model.}

\item{GEN1}{A function that takes a value of the unknown parameter in the
post-change model and return an observation.}

\item{th0}{Numeric: True value of the unknown parameter in the post-change
model.}

\item{detect.bayes}{A function that implements the Bayesian stopping rule.
Currently only \code{\link{detect.bayes}} for continuous data is supported.}

\item{nulower, nuupper}{Optional nonnegative numerics: The earliest and latest
time of changepoint based on prior belief. The default is \code{nulower=0}
and \code{nuupper=18} which corresponds to the geometric prior distribution
with \code{p=0.1}.}

\item{...}{Additional arguments to be passed to \code{detect.stat}.}
}
\value{
A named numeric vector with components
\enumerate{
  \item{\code{is.FA}} {A numeric of 1 or 0 indicating whether a false alarm
  has been raised.}
  \item{\code{DD}} {A positive numeric: The delay to
  detection.}
  \item{\code{CT}} {A positive numeric: The computation time.}
  \item{\code{LER}} {A numeric in \code{(0,1)}: The low ESS rate, i.e., the proportion of iterations that ESS drops below \code{c*n}.}
  \item{\code{AAR}} {A numeric in \code{(0,1)}: The average acceptance rate of the Metropolis-Hastings sampling in the move step. \code{NaN} if ESS never drops below \code{c*n}.}
  }
}
\description{
Simulation experiment of changepoint detection with unknown post-change
distributions using the Bayesian stopping rule.
}
\examples{
##Change from N(0,1) to 2*N(0,1)+1 occurs between 10 and 25.
##The mean and standard deviation of the post-change normal distribution are unknown.

GEN0=function() rnorm(1)
GEN1=function(th) th[2]*rnorm(1)+th[1]
ULP1=function(x,th) -(x-th[1])^2/2/th[2]^2
GLP1=function(x,th) -(x-th[1])/th[2]^2
LLP1=function(x,th) -1/th[2]^2
ULP0=function(x) ULP1(x,c(0,1))
GLP0=function(x) GLP1(x,c(0,1))
LLP0=function(x) LLP1(x,c(0,1))
par0=log(2*pi)/2
par1=function(th) log(2*pi)/2+log(th[2])

#using hyvarinen score
sim.detect.bayes(GEN0=GEN0,GEN1=GEN1,th0=c(1,2),detect.bayes=detect.bayes,nulower=10,nuupper=25,alpha=0.1,thlower=c(-Inf,0),GLP0=GLP0,LLP0=LLP0,GLP1=GLP1,LLP1=LLP1)

#using log score, normalizing constant unknown
sim.detect.bayes(GEN0=GEN0,GEN1=GEN1,th0=c(1,2),detect.bayes=detect.bayes,nulower=10,nuupper=25,alpha=0.1,thlower=c(-Inf,0),score="log",ULP0=ULP0,ULP1=ULP1,lenx=1)

#using log score, normalizing constant known
sim.detect.bayes(GEN0=GEN0,GEN1=GEN1,th0=c(1,2),detect.bayes=detect.bayes,nulower=10,nuupper=25,alpha=0.1,thlower=c(-Inf,0),score="log",ULP0=ULP0,ULP1=ULP1,par0=par0,par1=par1)
}
