#' Statistic-based Stopping Rule for Continuous Data
#'
#' Changepoint detection for continuous data with known post-change
#' distributions using the statistic-based stopping rule.
#'
#' @param GEN A function of time that returns an observation.
#' @param alpha A numeric parameter in \code{(0,1)} that controls the probability of
#'   false alarm.
#' @param nulower,nuupper Optional nonnegative numerics: The earliest and latest
#'   time of changepoint based on prior belief. The default is \code{nulower=0}
#'   and \code{nuupper=18} which corresponds to the geometric prior distribution
#'   with \code{p=0.1}.
#' @param score An optional character specifying the type of score to be used:
#'   The default \code{"hyvarinen"} or the conventional \code{"logarithmic"}.
#'   Can be abbreviated. Case insensitive.
#' @param ULP0,GLP0,LLP0,ULP1,GLP1,LLP1 Functions of an observation: The log
#'   unnormalized probability function, its gradient and its laplacian for the
#'   pre-change (\code{ULP0,GLP0,LLP0}) and post-change (\code{ULP1,GLP1,LLP1}) models. If
#'   \code{score="hyvarinen"}, either \{\code{GLP0,LLP0,GLP1,LLP1}\} or
#'   \{\code{ULP0,ULP1}\} is required. The former is recommended. In the latter case,
#'   \{\code{GLP0,LLP0,GLP1,LLP1}\} will be computed via \code{\link[pracma]{grad}} and
#'   \code{\link[pracma]{laplacian}}. If \code{score="logarithmic"}, only
#'   \{\code{ULP0,ULP1}\} is required.
#' @param par0,par1 Optional numeric parameters for the pre-change (\code{par0}) and
#'   post-change (\code{par1}) models. If \code{score="hyvarinen"}, the positive
#'   tuning parameter with a default of 1. If \code{score="logarithmic"}, the
#'   negative log normalizing constant. If omitted, will be computed via
#'   \code{\link[stats]{integrate}} (if \code{lenx=1}) or \code{\link[cubature]{hcubature}} (if
#'   \code{lenx>1}).
#' @param lenx A positive numeric: The length of the variable of an
#'   obervation. Optional if \code{score="hyvarinen"} or if
#'   \code{score="logarithmic"} and \code{par0,par1} are specified.
#' @param lower0,upper0,lower1,upper1 Optional numeric vectors of length \code{lenx}:
#'   The lower and upper limits of an observation from the pre-change (\code{lower0,upper0})
#'   and post-change (\code{lower1,upper1}) models. The defaults are infinite.
#' @return A positive numeric: The stopping time.
#' @examples
#' ##Change from N(0,1) to N(1,1) at t=15.
#' ##Prior knowledge suggests change occurs between 10 and 25.
#'
#' GEN=function(t) { if(15>=t) rnorm(1) else rnorm(1)+1 }
#' ULP0=function(x) -x^2/2
#' ULP1=function(x) -(x-1)^2/2
#' GLP0=function(x) -x
#' GLP1=function(x) -(x-1)
#' LLP0=function(x) -1
#' LLP1=function(x) -1
#' par0=log(2*pi)/2;par1=par0
#'
#' #using hyvarinen score
#' detect.stat(GEN=GEN,alpha=0.1,nulower=10,nuupper=25,GLP0=GLP0,LLP0=LLP0,GLP1=GLP1,LLP1=LLP1)
#'
#' #using log score. normalizing constant is unknown
#' detect.stat(GEN=GEN,alpha=0.1,nulower=10,nuupper=25,score="log",ULP0=ULP0,ULP1=ULP1,lenx=1)
#'
#' #using log score. normalizing constant is known
#' detect.stat(GEN=GEN,alpha=0.1,nulower=10,nuupper=25,score="log",ULP0=ULP0,ULP1=ULP1,par0=par0,par1=par1)
#' @export
detect.stat=function(GEN,alpha,nulower=NULL,nuupper=NULL,score="hyvarinen",
                     ULP0=NULL,GLP0=NULL,LLP0=NULL,ULP1=NULL,GLP1=NULL,LLP1=NULL,
                     par0=NULL,par1=NULL,
                     lenx=NULL,lower0=NULL,upper0=NULL,lower1=NULL,upper1=NULL)
{
  GEN=match.fun(GEN)
  # check the false alarm rate
  if(alpha <=0 | alpha>=1) stop("'alpha' should be in (0,1)")

  # suppose the time of changepoint nu follows a geometric prior distribution, compute the parameter p
  # the default is nulower=0, nuupper=18, which corresponds to p=0.1
  if(is.null(nulower)) nulower=0
  if(is.null(nuupper)) p=0.1 else{
    if(nulower >= nuupper) stop("need 'nulower' < 'nuupper'")
    p=1/(1+(nulower+nuupper)/2)
  }

  # chose score and compute related functions
  score=match.arg(tolower(score),c("hyvarinen","logarithmic"))
  if(score=="hyvarinen")
  {
    # some checks, preparations, and transformations of the log probablity functions
    if(is.null(GLP0) | is.null(LLP0)) ULP0=match.fun(ULP0)
    if(is.null(GLP0)) GLP0=function(x) pracma::grad(ULP0,x) else GLP0=match.fun(GLP0)
    if(is.null(LLP0)) LLP0=function(x) pracma::laplacian(ULP0,x) else LLP0=match.fun(LLP0)
    if(is.null(GLP1) | is.null(LLP1)) ULP1=match.fun(ULP1)
    if(is.null(GLP1)) GLP1=function(x) pracma::grad(ULP1,x) else GLP1=match.fun(GLP1)
    if(is.null(LLP1)) LLP1=function(x) pracma::laplacian(ULP1,x) else LLP1=match.fun(LLP1)

    # check the optional parameters for the pre and post change models
    if(is.null(par0)) par0=1 else if(par0 <=0)
    {
      warning("'par0' should be positive for hyvarinen score. Use par0=1")
      par0=1
    }
    if(is.null(par1)) par1=1 else if(par1 <=0)
    {
      warning("'par1' should be positive for hyvarinen score. Use par1=1")
      par1=1
    }

    # compute the scores for the pre and post change models
    SC0=function(x) par0*(sum(GLP0(x)^2)/2+LLP0(x))
    SC1=function(x) par1*(sum(GLP1(x)^2)/2+LLP1(x))
  }else
  {
    # when choose logarithmic score
    ULP0=match.fun(ULP0);ULP1=match.fun(ULP1)

    if(is.null(par0) | is.null(par1))
    {
      if(is.null(lenx)) stop("'lenx' is missing")
      if(lenx <=0) stop("'lenx' should be a positive integer")
      if(is.null(lower0)) lower0=rep(-Inf,lenx)
      if(is.null(lower1)) lower1=rep(-Inf,lenx)
      if(is.null(upper0)) upper0=rep(Inf,lenx)
      if(is.null(upper1)) upper1=rep(Inf,lenx)
    }
    if(is.null(par0))
    {
      fn=function(x) exp(ULP0(x))
      if(lenx==1) {par0=log(integrate(fn,lower0,upper0)$value)
      }else par0=log(cubature::hcubature(fn,lower0,upper0)$integral)
    }
    if(is.null(par1))
    {
      fn=function(x) exp(ULP1(x))
      if(lenx==1) {par1=log(integrate(fn,lower1,upper1)$value)
      }else par1=log(cubature::hcubature(fn,lower1,upper1)$integral)
    }

    SC0=function(x) -ULP0(x)+par0
    SC1=function(x) -ULP1(x)+par1
  }

  # define a threshold A
  logA=log(1-alpha)-log(alpha)-log(p)
  # define a statistic R, R=0 when t=0
  logR=-Inf
  t=1

  # statistic-based stopping rule
  repeat
  {
    x=GEN(t)
    if (sum(is.na(x))>0) {
      warning("did not detect any change")
      break
    }
    z=SC0(x)-SC1(x)
    # compute R recursively
    logR=log(1+exp(logR))+z-log(1-p)
    if(t>nulower & logR>=logA) break
    t=t+1
  }
  return(t)
}

flip=function(x,ULP,tbin){
  # Function to do binary flip
  #
  # Args:
  #   x: a vector
  #   ULP: a log unnormalized probability function
  #   tbin: a scalar
  #
  # Returns:
  #   Returns a vector, the ith component of which is ULP(tbin-x[i]).
  ULP=match.fun(ULP)
  lenx=length(x)
  vec=numeric(lenx)
  for(i in 1:lenx){
    y=x;y[i]=tbin-x[i]
    vec[i]=ULP(y)
  }
  return(vec)
}

#' Statistic-based Stopping Rule for Binary Data
#'
#' Changepoint detection for binary data with known post-change distributions
#' using the statistic-based stopping rule.
#'
#' @param GEN A function of time that returns an observation.
#' @param alpha A numeric parameter in \code{(0,1)} that controls the probability
#'   of false alarm.
#' @param nulower,nuupper Optional nonnegative numerics: The earliest and latest
#'   time of changepoint based on prior belief. The default is \code{nulower=0}
#'   and \code{nuupper=18} which corresponds to the geometric prior distribution
#'   with \code{p=0.1}.
#' @param score An optional character specifying the type of score to be used:
#'   The default \code{"hyvarinen"} or the conventional \code{"logarithmic"}.
#'   Can be abbreviated. Case insensitive.
#' @param ULP0,ULP1 Functions of an observation: The log unnormalized
#'   probability function for the pre-change (\code{ULP0}) and post-change
#'   (\code{ULP1}) models.
#' @param par0,par1 Optional numeric parameters for the pre-change (\code{par0})
#'   and post-change (\code{par1}) models. If \code{score="hyvarinen"}, the
#'   positive tuning parameter with a default of 1. If
#'   \code{score="logarithmic"}, the negative log normalizing constant. If
#'   omitted, will be computed by summing over the sample space.
#' @param lenx A positive numeric: The length of the variable of an
#'   obervation. Optional if \code{score="hyvarinen"} or if
#'   \code{score="logarithmic"} and \code{par0,par1} are specified.
#' @param tbin Optional numeric specifying the binary type: The default
#'   \code{tbin=1} representing \{\code{1,0}\} or the alternative \code{tbin=2}
#'   representing \{\code{1,-1}\}.
#' @return A positive numeric: The stopping time.
#' @examples
#' ##Change from 3 iid Bernoulli(0.2) to 3 iid Bernoulli(0.8) at t=10.
#' ##Prior knowledge suggests change occurs before 20.
#'
#' GEN=function(t) { if(10>=t) rbinom(3,1,0.2) else rbinom(3,1,0.8)}
#' ULP0=function(x) sum(x)*(log(0.2)-log(1-0.2))
#' ULP1=function(x) sum(x)*(log(0.8)-log(1-0.8))
#' par0=-3*log(1-0.2)
#' par1=-3*log(1-0.8)
#'
#' #using hyvarinen score
#' detect.bin.stat(GEN=GEN,alpha=0.1,nuupper=20,ULP0=ULP0,ULP1=ULP1)
#'
#' #using log score. normalizing constant is unknown
#' detect.bin.stat(GEN=GEN,alpha=0.1,nuupper=20,score="log",ULP0=ULP0,ULP1=ULP1,lenx=3)
#'
#' #using log score. normalizing constant is known
#' detect.bin.stat(GEN=GEN,alpha=0.1,nuupper=20,score="log",ULP0=ULP0,ULP1=ULP1,par0=par0,par1=par1)
#' @export
detect.bin.stat=function(GEN,alpha,nulower=NULL,nuupper=NULL,score="hyvarinen",
                         ULP0,ULP1,
                         par0=NULL,par1=NULL,
                         lenx=NULL,tbin=1)
{
  GEN=match.fun(GEN)
  if(alpha <=0 | alpha>=1) stop("'alpha' should be in (0,1)")

  if(is.null(nulower)) nulower=0
  if(is.null(nuupper)) p=0.1 else{
    if(nulower >= nuupper) stop("need 'nulower' < 'nuupper'")
    p=1/(1+(nulower+nuupper)/2)
  }

  score=match.arg(tolower(score),c("hyvarinen","logarithmic"))
  ULP0=match.fun(ULP0);ULP1=match.fun(ULP1)
  if(score=="hyvarinen")
  {

    if(is.null(par0)) par0=1 else if(par0 <=0)
    {
      warning("'par0' should be positive for hyvarinen score. Use par0=1")
      par0=1
    }
    if(is.null(par1)) par1=1 else if(par1 <=0)
    {
      warning("'par1' should be positive for hyvarinen score. Use par1=1")
      par1=1
    }

    SC0=function(x) par0*sum((1+exp(ULP0(x)-flip(x,ULP0,tbin)))^(-2))
    SC1=function(x) par1*sum((1+exp(ULP1(x)-flip(x,ULP1,tbin)))^(-2))
  }else
  {

    if(is.null(par0) | is.null(par1))
    {
      dom=as.matrix(expand.grid(rep(list(0:1),lenx)))
      if(tbin==2) dom[dom==0]=-1
    }
    if(is.null(par0)) par0=log(sum(exp(apply(dom,1,ULP0))) )
    if(is.null(par1)) par1=log(sum(exp(apply(dom,1,ULP1))) )

    SC0=function(x) -ULP0(x)+par0
    SC1=function(x) -ULP1(x)+par1
  }
  logA=log(1-alpha)-log(alpha)-log(p)
  logR=-Inf
  t=1

  repeat
  {
    x=GEN(t)
    if (is.na(sum(x))) {
      warning("did not detect any change")
      break
      }
    z=SC0(x)-SC1(x)
    logR=log(1+exp(logR))+z-log(1-p)
    if(t>nulower & logR>=logA) break
    t=t+1
  }
  return(t)
}

#' Simulation using Statistic-based Stopping Rule
#'
#' Simulation experiment of changepoint detection with known post-change
#' distributions using the statistic-based stopping rule.
#'
#' @param GEN0,GEN1 Functions that take no argument and return an observation from the
#'   pre-change (\code{GEN0}) and post-change (\code{GEN1}) models.
#' @param detect.stat A function that implements the statistic-based stopping rule: \code{\link{detect.stat}} for continuous data or \code{\link{detect.bin.stat}} for binary data.
#' @param nulower,nuupper Optional nonnegative numerics: The earliest and latest
#'   time of changepoint based on prior belief. The default is \code{nulower=0}
#'   and \code{nuupper=18} which corresponds to the geometric prior distribution
#'   with \code{p=0.1}.
#' @param ... Additional arguments to be passed to \code{detect.stat}.
#' @return A named numeric vector with components
#' \enumerate{
#'   \item{\code{is.FA}} {A numeric of 1 or 0 indicating whether a false alarm has been raised.}
#'   \item{\code{DD}}  {A positive numeric: The delay to detection.}
#'   \item{\code{CT}}  {A positive numeric: The computation time.}
#' }
#' @examples
#' ##Change from N(0,1) to N(1,1) occurs between 10 and 25.
#'
#' GEN0=function() rnorm(1)
#' GEN1=function() rnorm(1)+1
#' ULP0=function(x) -x^2/2
#' ULP1=function(x) -(x-1)^2/2
#' GLP0=function(x) -x
#' GLP1=function(x) -(x-1)
#' LLP0=function(x) -1
#' LLP1=function(x) -1
#' par0=log(2*pi)/2;par1=par0
#'
#' #using hyvarinen score
#' sim.detect.stat(GEN0=GEN0,GEN1=GEN1,detect.stat=detect.stat,nulower=10,nuupper=25,alpha=0.1,GLP0=GLP0,LLP0=LLP0,GLP1=GLP1,LLP1=LLP1)
#'
#' #using log score. normalizing constant is unknown
#' sim.detect.stat(GEN0=GEN0,GEN1=GEN1,detect.stat=detect.stat,nulower=10,nuupper=25,alpha=0.1,score="log",ULP0=ULP0,ULP1=ULP1,lenx=1)
#'
#' #using log score. normalizing constant is known
#' sim.detect.stat(GEN0=GEN0,GEN1=GEN1,detect.stat=detect.stat,nulower=10,nuupper=25,alpha=0.1,score="log",ULP0=ULP0,ULP1=ULP1,par0=par0,par1=par1)
#'
#'
#' ###################################
#'
#' ##Change from 3 iid Bernoulli(0.2) to 3 iid Bernoulli(0.8) occurs before 20.
#'
#' GEN0=function() rbinom(3,1,0.2)
#' GEN1=function() rbinom(3,1,0.8)
#' ULP0=function(x) sum(x)*(log(0.2)-log(1-0.2))
#' ULP1=function(x) sum(x)*(log(0.8)-log(1-0.8))
#' par0=-3*log(1-0.2)
#' par1=-3*log(1-0.8)
#'
#' #using hyvarinen score
#' sim.detect.stat(GEN0=GEN0,GEN1=GEN1,detect.stat=detect.bin.stat,nuupper=20,alpha=0.1,ULP0=ULP0,ULP1=ULP1)
#'
#' #using log score. normalizing constant is unknown
#' sim.detect.stat(GEN0=GEN0,GEN1=GEN1,detect.stat=detect.bin.stat,nuupper=20,alpha=0.1,score="log",ULP0=ULP0,ULP1=ULP1,lenx=3)
#'
#' #using log score. normalizing constant is known
#' sim.detect.stat(GEN0=GEN0,GEN1=GEN1,detect.stat=detect.bin.stat,nuupper=20,alpha=0.1,score="log",ULP0=ULP0,ULP1=ULP1,par0=par0,par1=par1)
#' @export
sim.detect.stat=function(GEN0,GEN1,detect.stat,nulower=NULL,nuupper=NULL,...)
{
  GEN0=match.fun(GEN0);GEN1=match.fun(GEN1)
  detect.stat=match.fun(detect.stat)

  if(is.null(nulower)) nulower=0
  if(is.null(nuupper)) p=0.1 else{
    if(nulower >= nuupper) stop("need 'nulower' < 'nuupper'")
    p=1/(1+(nulower+nuupper)/2)
  }
  nu=rgeom(1,p)+nulower
  GEN=function(t)
  {
    if(nu>=t) GEN0() else GEN1()
  }
  CT0=proc.time()
  t=detect.stat(GEN=GEN,nulower=nulower,nuupper=nuupper,...)

  # compute the running time of function detect.stat
  CT=proc.time()-CT0
  out=c((nu>=t),max(t-nu,0),CT[1])
  names(out)=c("is.FA","DD","CT")
  return(out)
}
