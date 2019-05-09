#' Simulation using Bayesian Stopping Rule with unknown post-change
#'
#' Simulation experiment of changepoint detection with unknown post-change
#' distributions using the Bayesian stopping rule.
#'
#' @param GEN0 A function that takes a value of the unknown parameter in the
#'   pre-change model and return an observation.
#' @param GEN1 A function that takes a value of the unknown parameter in the
#'   post-change model and return an observation.
#' @param th0 Numeric: True value of the unknown parameter in the pre-change
#'   model.
#' @param th1 Numeric: True value of the unknown parameter in the post-change
#'  model.
#' @param detect.bayes.unknown.pre A function that implements the Bayesian stopping rule.
#' @param nulower,nuupper Optional nonnegative numerics: The earliest and latest
#'   time of changepoint based on prior belief. The default is \code{nulower=0}
#'   and \code{nuupper=18} which corresponds to the geometric prior distribution
#'   with \code{p=0.1}.
#' @param ... Additional arguments to be passed to \code{detect.stat.unknown.pre}.
#' @return A named numeric vector with components
#' \enumerate{
#'   \item{\code{is.FA}} {A numeric of 1 or 0 indicating whether a false alarm
#'   has been raised.}
#'   \item{\code{DD}} {A positive numeric: The delay to
#'   detection.}
#'   \item{\code{CT}} {A positive numeric: The computation time.}
#'   \item{\code{LER}} {A numeric in \code{(0,1)}: The low ESS rate, i.e., the proportion of iterations that ESS drops below \code{c*n}.}
#'   \item{\code{AAR}} {A numeric in \code{(0,1)}: The average acceptance rate of the Metropolis-Hastings sampling in the move step. \code{NaN} if ESS never drops below \code{c*n}.}
#'   }
#' @examples
#' #Change from N(0,1) to 2*N(0,1)+1 occurs between 10 and 25.
#' #The mean and standard deviation of the post-change normal distribution are unknown.
#'
#' GEN1=function(th) th[2]*rnorm(1)+th[1]
#' GEN0=function(th) th[2]*rnorm(1)+th[1]
#' ULP1=function(x,th) -(x-th[1])^2/2/th[2]^2
#' GLP1=function(x,th) -(x-th[1])/th[2]^2
#' LLP1=function(x,th) -1/th[2]^2
#' ULP0=function(x,th) -(x-th[1])^2/2/th[2]^2
#' GLP0=function(x,th) -(x-th[1])/th[2]^2
#' LLP0=function(x,th) -1/th[2]^2
#' par0=function(th) log(2*pi)/2+log(th[2])
#' par1=function(th) log(2*pi)/2+log(th[2])
#'
#' #using hyvarinen score
#' sim.detect.bayes.unknown.pre(GEN0=GEN0,GEN1=GEN1,th0=c(0,1),th1=c(1,2),detect.bayes.unknown.pre=detect.bayes.unknown.pre,alpha=0.3,nulower=5,nuupper=30,lenth1=2, lenth0=2,thlower1=c(-Inf,0),thlower0=c(-Inf,0),score="hyvarinen",GLP0=GLP0,LLP0=LLP0,GLP1=GLP1,LLP1=LLP1,par1=0.1,par0=0.1)
#'
#' #using log score, normalizing constant unknown
#' sim.detect.bayes.unknown.pre(GEN0=GEN0,GEN1=GEN1,th0=c(0,1),th1=c(1,2),detect.bayes.unknown.pre=detect.bayes.unknown.pre,nulower=20,nuupper=50,alpha=0.1,lenth1=2,lenth0=2,thlower1=c(-Inf,0),thlower0=c(-Inf,0),score="log",lenx=1,ULP0=ULP0,ULP1=ULP1)
#'
#'
#' #using log score, normalizing constant known
#' sim.detect.bayes.unknown.pre(GEN0=GEN0,GEN1=GEN1,th0=c(0,1),th1=c(1,2),detect.bayes.unknown.pre=detect.bayes.unknown.pre,nulower=20,nuupper=50,alpha=0.1,lenth1=2,lenth0=2,thlower1=c(-Inf,0),thlower0=c(-Inf,0),score="log",lenx=1,ULP0=ULP0,ULP1=ULP1,par0=par0,par1=par1)
#' @export
sim.detect.bayes.unknown.pre=function(GEN0,GEN1,th0,th1,detect.bayes.unknown.pre, nulower=NULL,nuupper=NULL,
                             lenth0, lenth1,...)
{
  GEN0=match.fun(GEN0);GEN1=match.fun(GEN1)
  detect.bayes.unknown.pre=match.fun(detect.bayes.unknown.pre)
  if(is.null(nulower)) nulower=0
  if(is.null(nuupper)) p=0.1 else{
    if(nulower >= nuupper) stop("need 'nulower' < 'nuupper'")
    p=1/(1+(nulower+nuupper)/2)
  }
  nu=rgeom(1,p)+nulower
  GEN=function(t)
  {
    if(nu>=t) GEN0(th0) else GEN1(th1)
  }

  CT0=proc.time()
  db=detect.bayes.unknown.pre(GEN=GEN,nulower=nulower,nuupper=nuupper,lenth1=lenth1,lenth0 = lenth0,...)
  CT=proc.time()-CT0
  out=c((nu>=db[1]),max(db[1]-nu,0),CT[1],db[-1])
  names(out)=c("is.FA","DD","CT","LER","AAR")
  return(out)
}
