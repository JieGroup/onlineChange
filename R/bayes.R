# transformation from finite to infinite
f2i=function(x,a,b){
  if(is.infinite(a) & is.infinite(b)) return(x)
  if(is.infinite(b)) return(log(x-a))
  if(is.infinite(a)) return(log(b-x))
  return(log(x-a)-log(b-x))
}

# transformation from infinite to finite
i2f=function(y,a,b){
  if(is.infinite(a) & is.infinite(b)) return(y)
  if(is.infinite(b)) return(exp(y)+a)
  if(is.infinite(a)) return(b-exp(y))
  return( a+(b-a)/(1+exp(-y)) )
}

# look at https://mc-stan.org/docs/2_18/reference-manual/variable-transforms-chapter.html
# for more details
logJ=function(y,a,b)
{
  if(is.infinite(a) & is.infinite(b)) return(0)
  if(is.infinite(b) | is.infinite(a)) return(y)
  return(-y-2*log(1+exp(-y)))
}

v2=function(xoy,a,b,FUN){
  # Function to apply transformation to three vectors.
  #
  # Args:
  #   xoy: a vector
  #   a: a vector that no shorter than xoy
  #   b: a vector that no shorter than xoy
  #   FUN: an already speified function that have 3 inputs
  #
  # Returns:
  #   Returns a vector, the ith component of which is FUN(xoy[i],a[i],b[i]).
  # copy to the local scope
  FUN=match.fun(FUN)
  n=length(xoy)
  if(n==1) return(FUN(xoy,a,b))
  yox=numeric(n)
  for(i in 1:n) yox[i]=FUN(xoy[i],a[i],b[i])
  return(yox)
}

m2=function(mat,a,b,FUN){
  # Function to apply transformation to a matrix.
  #
  # Args:
  #   mat: a matrix
  #   a: a scalar
  #   b: a scalar
  #   FUN: an already speified function that have 3 inputs
  #
  # Returns:
  #   Returns a matrix, the ith row of which is v2(mat[i,],a,b,FUN).
  FUN=match.fun(FUN)
  if(length(a)==1) return(sapply(mat,FUN,a,b))
  out=mat
  for(i in 1:nrow(mat)) out[i,]=v2(mat[i,],a,b,FUN)
  return(out)
}


f2=function(y,a,b,ULP){
  # Function to apply transformation to variables and give dditional Jacobian to prob. functions
  #
  # Args:
  #   y: a vector
  #   a: a vector that no shorter than y
  #   b: a vector that no shorter than y
  #   ULP: an already speified function, the log unnormalized probability function of a certain distribution
  #
  # Returns:
  #   Returns a vector that equals to the prob. functions of the transformed variables plus an additional Jacobian.
  ULP=match.fun(ULP)
  # the ith component of vector x is i2f(y[i],a[i],b[i]), i.e. a[i]+(b[i]-a[i])/(1+exp(-y[i]))
  x=v2(y,a,b,i2f)
  # the ith component of vector v2(y,a,b,logJ)) is logJ(y[i],a[i],b[i])
  return(ULP(x)+sum(v2(y,a,b,logJ)))
}

mult=function(w,n){
  # Function for multinomial resampling.
  #
  # Args:
  #   w: numeric non-negative vector of length K, specifying the probability for the K classes,
  #      is internally normalized to sum 1
  #   n: integer equals to K, specifying the total number of objects that are put into K boxes
  #      in the typical  multinomial experiment
  #
  # Returns:
  #   Returns a vector, which is the result of multinomial resampling.
  tms=c(rmultinom(1,n,w))
  return(rep(1:n,tms))
}


#' Bayesian Stopping Rule for Continuous Data
#'
#' Changepoint detection for continuous data with unknown post-change
#' distributions using the Bayesian stopping rule.
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
#' @param c,n Optional parameters of the Sequentital Monte Carlo algorithm: ESS
#'   threshold \code{0<c<1} and sample size \code{n>0}. Default is \code{c=0.5}
#'   and \code{n=1000}. The resample-move step is triggered whenever the ESS
#'   goes below \code{c*n}.
#' @param lenth A positive numeric: The length of the variable of the unknown
#'   parameter in the post-change model.
#' @param thlower,thupper Optional numeric vectors of length \code{lenth}: The
#'   lower and upper limits of the unknown parameter in the post-change model.
#'   The defaults are infinite.
#' @param GENTH An optional function that takes a sample size and returns a
#'   random sample from the prior distribution of the unknown parameter in the
#'   post-change model. Default is standard normal on the unconstrained
#'   space. Required if \code{ULPTH} is specified.
#' @param ULPTH An optional function: The log unnormalized probability function
#'   of the prior distribution of the unknown parameter in the post-change
#'   model. Default is standard normal on the unconstrained space. Required
#'   if \code{GENTH} is specified.
#' @param ULP0,GLP0,LLP0 Functions of an observation: The log unnormalized
#'   probability function, its gradient and its laplacian for the pre-change
#'   model. If \code{score="hyvarinen"}, either \{\code{GLP0,LLP0}\} or \code{ULP0}
#'   is required. The former is recommended. In the latter case,
#'   \{\code{GLP0,LLP0}\} will be computed via \code{\link[pracma]{grad}} and
#'   \code{\link[pracma]{laplacian}}. If \code{score="logarithmic"}, only
#'   \code{ULP0} is required.
#' @param ULP1,GLP1,LLP1 Functions of an observation and a numeric parameter:
#'   The log unnormalized probability function, its gradient and its laplacian
#'   for the post-change model. If \code{score="hyvarinen"}, either
#'   \{\code{GLP1,LLP1}\} or \code{ULP1} is required. The former is recommended. In
#'   the latter case, \{\code{GLP1,LLP1}\} will be computed via
#'   \code{\link[pracma]{grad}} and \code{\link[pracma]{laplacian}}. If
#'   \code{score="logarithmic"}, only \code{ULP1} is required.
#' @param par0,par1 Optional numeric parameters for the pre-change
#'   (\code{par0}) and post-change (\code{par1}) models, except if
#'   \code{score="logarithmic"} that \code{par1} is a function of the unknown
#'   parameter in the post-change model. If \code{score="hyvarinen"}, the positive
#'   tuning parameter with a default of 1. If \code{score="logarithmic"}, the
#'   negative log normalizing constant. If omitted, will be computed via
#'   \code{\link[stats]{integrate}} (if \code{lenx=1}) or
#'   \code{\link[cubature]{hcubature}} (if \code{lenx>1}).
#' @param lenx A positive numeric: The length of the variable of an obervation.
#'   Optional if \code{score="hyvarinen"} or if \code{score="logarithmic"} and
#'   \code{par0,par1} are specified.
#' @param lower0,upper0,lower1,upper1 Optional numeric vectors of length
#'   \code{lenx}: The lower and upper limits of an observation from the
#'   pre-change (\code{lower0,upper0}) and post-change (\code{lower1,upper1})
#'   models. The defaults are infinite.
#' @return A named numeric vector with components
#' \enumerate{
#'   \item{\code{t}} {A positive numeric: The stopping time.}
#'   \item{\code{LER}} {A numeric in \code{(0,1)}: The low ESS rate, i.e., the proportion of iterations that ESS drops below \code{c*n}.}
#'   \item{\code{AAR}} {A numeric in \code{(0,1)}: The average acceptance rate of the Metropolis-Hastings sampling in the move step. \code{NaN} if ESS never drops below \code{c*n}.}
#' }
#' @examples
#' ##Change from N(0,1) to 2*N(0,1)+1 at t=15.
#' ##Prior knowledge suggests change occurs between 10 and 25.
#' ##The mean and standard deviation of the post-change normal distribution are unknown.
#'
#' GEN=function(t){ if(15>=t) rnorm(1) else 2*rnorm(1)+1 }
#' ULP1=function(x,th) -(x-th[1])^2/2/th[2]^2
#' GLP1=function(x,th) -(x-th[1])/th[2]^2
#' LLP1=function(x,th) -1/th[2]^2
#' ULP0=function(x) ULP1(x,c(0,1))
#' GLP0=function(x) GLP1(x,c(0,1))
#' LLP0=function(x) LLP1(x,c(0,1))
#' par0=log(2*pi)/2
#' par1=function(th) log(2*pi)/2+log(th[2])
#'
#' #using hyvarinen score
#' detect.bayes(GEN=GEN,alpha=0.1,nulower=10,nuupper=25,lenth=2,thlower=c(-Inf,0),GLP0=GLP0,LLP0=LLP0,GLP1=GLP1,LLP1=LLP1)
#'
#' #using log score, normalizing constant unknown
#' detect.bayes(GEN=GEN,alpha=0.1,nulower=10,nuupper=25,lenth=2,thlower=c(-Inf,0),score="log",ULP0=ULP0,ULP1=ULP1,lenx=1)
#'
#' #using log score, normalizing constant known
#' detect.bayes(GEN=GEN,alpha=0.1,nulower=10,nuupper=25,lenth=2,thlower=c(-Inf,0),score="log",ULP0=ULP0,ULP1=ULP1,par0=par0,par1=par1)
#' @export
detect.bayes=function(GEN,alpha,nulower=NULL,nuupper=NULL,score="hyvarinen",c=0.5,n=1000,
                      lenth,thlower=NULL,thupper=NULL,
                      GENTH=NULL,ULPTH=NULL,
                      ULP0=NULL,GLP0=NULL,LLP0=NULL,ULP1=NULL,GLP1=NULL,LLP1=NULL,
                      par0=NULL,par1=NULL,
                      lenx=NULL,lower0=NULL,upper0=NULL,lower1=NULL,upper1=NULL)
{
  GEN=match.fun(GEN)
  # check the false alarm rate
  if(alpha <=0 | alpha>=1) stop("'alpha' should be in (0,1)")
  # choose score
  score=match.arg(tolower(score),c("hyvarinen","logarithmic"))

  # suppose the time of changepoint nu follows a geometric prior distribution, compute the parameter p
  # the default is nulower=0, nuupper=18, which corresponds to p=0.1
  if(is.null(nulower)) nulower=0
  if(is.null(nuupper)) p=0.1 else{
    if(nulower >= nuupper) stop("need 'nulower' < 'nuupper'")
    p=1/(1+(nulower+nuupper)/2)
  }

  # suppose the unknown parameters in the post-change model th = c(th[1], th[2], ...) follows some known prior distribution
  # if not specified, the defaults of the lower and upper limits of th are -Inf and Inf
  if(is.null(thlower)) thlower=rep(-Inf,lenth)
  if(is.null(thupper)) thupper=rep(Inf,lenth)

  # if not specified, the default prior distribution of th is standard normal N(0,1) on the unconstrained space
  if(is.null(GENTH) & is.null(ULPTH)){
    # give the log unnormalized probability function of the default prior distribution of th
    ULPTH=function(th) -sum(th^2)/2
    # give a random sample of size n from the default prior distribution of th, a matrix with n rows and lenth columns
    GENTH=function(n) matrix(rnorm(n*lenth),nrow=n)

    # create a matrix nuth that combines 2 parts by column:
    # a random sample of the time of changepoint nu, a vector of length n
    # a random sample of the unknown parameters in the post-change model th, a matrix with n rows and lenth columns
    nuth=cbind(nulower+rgeom(n,p),GENTH(n))
  }
  else{
    # GENTH and ULPTH should always exist or be  missing at the same time
    if(is.null(GENTH)) stop("'GENTH' is missing")
    if(is.null(ULPTH)) stop("'ULPTH' is missing")

    # when both GENTH and ULPTH are specified
    fULPTH=match.fun(ULPTH)
    # apply transformation f2 to ULPTH, which adds an additional Jacobian to ULPTH
    ULPTH=function(th) f2(th,thlower,thupper,fULPTH)
    GENTH=match.fun(GENTH)

    # when building nuth, apply transformation m2 to GENTH(n), which makes finite to infinite transformations to GENTH(n)
    nuth=cbind(nulower+rgeom(n,p),m2(GENTH(n),thlower,thupper,f2i))
  }

  # compute the log of Metropols-Hastings kernal with independent proposal from Geom(1/(1+numean)) * N(thmean, thsd^2)
  # the input nuth here is a vector of length lenth+1
  LPQ=function(nuth,numean,thmean,thsd){
    # sample nu from its prior distribution
    nu=nuth[1]
    # sample th from its prior distribution
    th=nuth[-1]
    (nu-nulower)*(log(1-p)-log(1-1/(1+numean)))+ULPTH(th)+sum((th-thmean)^2/thsd^2)/2
  }

  # when choose hyvarinen score
  if(score=="hyvarinen")
  {
    # some checks, preparations, and transformations of the log probablity functions
    # functions for pre-change:
    if(is.null(GLP0) | is.null(LLP0)) ULP0=match.fun(ULP0)
    if(is.null(GLP0)) GLP0=function(x) pracma::grad(ULP0,x) else GLP0=match.fun(GLP0)
    if(is.null(LLP0)) LLP0=function(x) pracma::laplacian(ULP0,x) else LLP0=match.fun(LLP0)

    # functions for post-change:
    if(is.null(GLP1) | is.null(LLP1)){
      fULP1=match.fun(ULP1)
      ULP1=function(x,th) fULP1(x,v2(th,thlower,thupper,i2f))
    }
    if(is.null(GLP1)) GLP1=function(x,th) pracma::grad(ULP1,x,th=th) else{
      fGLP1=match.fun(GLP1)
      GLP1=function(x,th) fGLP1(x,v2(th,thlower,thupper,i2f))
    }
    if(is.null(LLP1)) LLP1=function(x,th) pracma::laplacian(ULP1,x,th=th) else{
      fLLP1=match.fun(LLP1)
      LLP1=function(x,th) fLLP1(x,v2(th,thlower,thupper,i2f))
    }

    # check the optional parameters for the pre and post change models
    if(is.null(par0)) par0=1 else if(par0 <=0){
      warning("'par0' should be positive for hyvarinen score. Use par0=1")
      par0=1
    }
    if(is.null(par1)) par1=1 else if(par1 <=0){
      warning("'par1' should be positive for hyvarinen score. Use par1=1")
      par1=1
    }

    # compute the scores for the pre and post change models
    SC0=function(x) par0*(sum(GLP0(x)^2)/2+LLP0(x))
    SC1=function(x,th) par1*(sum(GLP1(x,th)^2)/2+LLP1(x,th))
    rSC1=function(th,x) SC1(x,th)
    hSC1=function(x,th){
      if(lenth==1) sapply(th,rSC1,x) else apply(th,1,rSC1,x)
    }

    # compute the difference of the pre and post change scores
    hLZSUM=function(nuth,t,x){
      out=numeric(n)
      # the vector of random samples of nu
      nu=nuth[,1]
      # choose the random samples of nu that smaller than t, i.e. change occurs before t
      id=which(nu < t)
      for(i in id){
        xtail=x
        if(nu[i]>0) xtail=x[-(1:nu[i])]
        out[i]=sum(sapply(xtail,SC0))-sum(sapply(xtail,SC1,nuth[i,-1]))
      }
      return(out)
    }
  }else
  {
    # when choose logarithmic score
    ULP0=match.fun(ULP0)
    fULP1=match.fun(ULP1)
    ULP1=function(x,th) fULP1(x,v2(th,thlower,thupper,i2f))
    if(is.null(par0) | is.null(par1)){
      if(is.null(lenx)) stop("'lenx' is missing")
      if(lenx <=0) stop("'lenx' should be a positive integer")
      if(is.null(lower0)) lower0=rep(-Inf,lenx)
      if(is.null(lower1)) lower1=rep(-Inf,lenx)
      if(is.null(upper0)) upper0=rep(Inf,lenx)
      if(is.null(upper1)) upper1=rep(Inf,lenx)
    }
    if(is.null(par0)){
      fn=function(x) exp(ULP0(x))
      if(lenx==1) {par0=log(integrate(fn,lower0,upper0)$value)
      }else par0=log(cubature::hcubature(fn,lower0,upper0)$integral)
    }
    if(is.null(par1)){
      fn=function(x,th) exp(ULP1(x,th))
      if(lenx==1) {par1=function(th) log(integrate(fn,lower1,upper1,th)$value)
      }else par1=function(th) log(cubature::hcubature(fn,lower1,upper1,th)$integral)
    }else{
      fpar1=match.fun(par1)
      par1=function(th) fpar1(v2(th,thlower,thupper,i2f))
    }

    if(lenth==1) lpar1=sapply(nuth[,-1],par1) else lpar1=apply(nuth[,-1],1,par1)

    SC0=function(x) -ULP0(x)+par0
    SC1=function(x,th,par1th) -ULP1(x,th)+par1th

    #### ?? potential bug ?? ####
    lSC1=function(x,th,lpar1){
      out=numeric(length(lpar1))
      #### when length(lpar1)==1, th is a vector, th[i,] would return error.
      #### SOL: add if(length(lpar1)==1) part
      if(length(lpar1)==1){
        out=SC1(x,th,lpar1)
      }else{
        for(i in 1:length(lpar1)){
          if(lenth==1) out[i]=SC1(x,th[i],lpar1[i]) else out[i]=SC1(x,th[i,],lpar1[i])
        }
      }
      return(out)
    }

    lLZSUM=function(nuth,t,x,lpar1){
      out=numeric(n)
      nu=nuth[,1]
      id=which(nu < t)
      for(i in id)
      {
        xtail=x
        if(nu[i]>0) xtail=x[-(1:nu[i])]
        out[i]=sum(sapply(xtail,SC0))-sum(sapply(xtail,SC1,th=nuth[i,-1],par1th=lpar1[i]))
      }
      return(out)
    }
  }

  # At step t=0:
  t=0
  # initialize some useful variables
  # w: weights
  w=rep(1/n,n)
  # lw: the log of weights
  lw=rep(-log(n),n)
  # lrat: the log hastings ratio
  lrat=numeric(n)
  # ngsc: the negative scores
  ngsc=numeric(n)
  # lowESS: low Effective Sample Size
  lowESS=0
  # ar: acceptance rate
  ar=numeric()
  # x: an observation sampled from the pre-change model
  x=numeric()

  # At step t = 2, 3, ...:
  repeat
  {
    # if ESS = 1/sum(w^2) < c*n, do the resample and move step and draw (nu,th) for time = t
    if(1/sum(w^2) < c*n)
    {
      lowESS=lowESS+1

      # resample step:
      # use a multinomial distribution to obtain equally weighted samples
      id=mult(w,n)
      nuth=nuth[id,]
      if(score=="logarithmic") lpar1=lpar1[id]
      # set the weights to be equal
      w=rep(1/n,n);lw=rep(-log(n),n)

      # move step:
      # consider an MCMC kernel targeting the joint posterior of (nu,th) conditional on the t-1 observations
      # compute the mean of randomly sampled nu
      numean=mean(nuth[,1])
      # draw a new nu from the MCMC kernel
      nu.p=nulower+rgeom(n,1/(1+(numean-nulower)))
      # draw a new th from the MCMC kernel
      if(lenth==1){
        thmean=mean(nuth[,-1]);thsd=sd(nuth[,-1])
        th.p=rnorm(n)*thsd+thmean
      }else{
        thmean=apply(nuth[,-1],2,mean);thsd=apply(nuth[,-1],2,sd)
        th.p=matrix(rnorm(n*lenth)*rep(thsd,n)+rep(thmean,n),nrow=n,byrow=T)
      }
      # combine the new nu and th together
      nuth.p=cbind(nu.p,th.p)
      # compute the log Hastings ratio
      if(score=="logarithmic"){
        if(lenth==1) lpar1.p=sapply(th.p,par1) else lpar1.p=apply(th.p,1,par1)
      }
      lrat=apply(nuth.p,1,LPQ,numean,thmean,thsd)-apply(nuth,1,LPQ,numean,thmean,thsd)
      if(score=="hyvarinen"){ lrat=lrat+hLZSUM(nuth.p,t,x)-hLZSUM(nuth,t,x)
      } else lrat=lrat+lLZSUM(nuth.p,t,x,lpar1.p)-lLZSUM(nuth,t,x,lpar1)

      # accept the new (nu,th) if the Hastings ratio smaller than a random sample from Unif(0,1)
      acc=(lrat >= log(runif(n)))
      if(sum(acc)>0){
        nuth[acc,]=nuth.p[acc,]
        if(score=="logarithmic") lpar1[acc]=lpar1.p[acc]
      }
      ar=c(ar,mean(acc))
    }

    t=t+1
    xt=GEN(t)

    if(is.na(xt)) stop("Did not detect any change.")

    x=c(x,xt)
    xg=(nuth[,1] <t)

    if(sum(!xg)>0) ngsc[!xg]=-SC0(xt)
    if(sum(xg)>0){
      if(score=="hyvarinen") {ngsc[xg]=-hSC1(xt,nuth[xg,-1])
      }else ngsc[xg]=-lSC1(xt,nuth[xg,-1],lpar1[xg])   #### ?? lSC1 ??####
    }
    lw.ngsc=lw+ngsc
    diff=lw.ngsc-max(lw.ngsc)
    # calclute the new weights w
    lw=diff-log(sum(exp(diff)))
    w=exp(lw)
    # now the new (w, nu, th) approximate the joint posterior of (nu, th) conditional on the t observations

    #### ?? potential bug: sum(w[xg]) could be NA, i.e. all nu > t
    #### SOL: when sum(w[xg])=NA, set it to 0
    sum.w = sum(w[xg])
    if(is.na(sum.w)==TRUE){
      sum.w = 0
    }

    if(t>nulower & sum.w>=1-alpha) break
  }


  out=c(t,lowESS/t,mean(ar))
  names(out)=c("t","LER","AAR")
  return(out)
}

#' Simulation using Bayesian Stopping Rule
#'
#' Simulation experiment of changepoint detection with unknown post-change
#' distributions using the Bayesian stopping rule.
#'
#' @param GEN0 A function that take no argument and return an observation from
#'   the pre-change model.
#' @param GEN1 A function that takes a value of the unknown parameter in the
#'   post-change model and return an observation.
#' @param th0 Numeric: True value of the unknown parameter in the post-change
#'   model.
#' @param detect.bayes A function that implements the Bayesian stopping rule.
#'   Currently only \code{\link{detect.bayes}} for continuous data is supported.
#' @param nulower,nuupper Optional nonnegative numerics: The earliest and latest
#'   time of changepoint based on prior belief. The default is \code{nulower=0}
#'   and \code{nuupper=18} which corresponds to the geometric prior distribution
#'   with \code{p=0.1}.
#' @param ... Additional arguments to be passed to \code{detect.stat}.
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
#' ##Change from N(0,1) to 2*N(0,1)+1 occurs between 10 and 25.
#' ##The mean and standard deviation of the post-change normal distribution are unknown.
#'
#' GEN0=function() rnorm(1)
#' GEN1=function(th) th[2]*rnorm(1)+th[1]
#' ULP1=function(x,th) -(x-th[1])^2/2/th[2]^2
#' GLP1=function(x,th) -(x-th[1])/th[2]^2
#' LLP1=function(x,th) -1/th[2]^2
#' ULP0=function(x) ULP1(x,c(0,1))
#' GLP0=function(x) GLP1(x,c(0,1))
#' LLP0=function(x) LLP1(x,c(0,1))
#' par0=log(2*pi)/2
#' par1=function(th) log(2*pi)/2+log(th[2])
#'
#' #using hyvarinen score
#' sim.detect.bayes(GEN0=GEN0,GEN1=GEN1,th0=c(1,2),detect.bayes=detect.bayes,nulower=10,nuupper=25,alpha=0.1,thlower=c(-Inf,0),GLP0=GLP0,LLP0=LLP0,GLP1=GLP1,LLP1=LLP1)
#'
#' #using log score, normalizing constant unknown
#' sim.detect.bayes(GEN0=GEN0,GEN1=GEN1,th0=c(1,2),detect.bayes=detect.bayes,nulower=10,nuupper=25,alpha=0.1,thlower=c(-Inf,0),score="log",ULP0=ULP0,ULP1=ULP1,lenx=1)
#'
#' #using log score, normalizing constant known
#' sim.detect.bayes(GEN0=GEN0,GEN1=GEN1,th0=c(1,2),detect.bayes=detect.bayes,nulower=10,nuupper=25,alpha=0.1,thlower=c(-Inf,0),score="log",ULP0=ULP0,ULP1=ULP1,par0=par0,par1=par1)
#' @export
sim.detect.bayes=function(GEN0,GEN1,th0,detect.bayes,nulower=NULL,nuupper=NULL,...)
{
  GEN0=match.fun(GEN0);GEN1=match.fun(GEN1)
  detect.bayes=match.fun(detect.bayes)
  if(is.null(nulower)) nulower=0
  if(is.null(nuupper)) p=0.1 else{
    if(nulower >= nuupper) stop("need 'nulower' < 'nuupper'")
    p=1/(1+(nulower+nuupper)/2)
  }
  nu=rgeom(1,p)+nulower
  GEN=function(t)
  {
    if(nu>=t) GEN0() else GEN1(th0)
  }
  CT0=proc.time()
  db=detect.bayes(GEN=GEN,nulower=nulower,nuupper=nuupper,lenth=length(th0),...)
  CT=proc.time()-CT0
  out=c((nu>=db[1]),max(db[1]-nu,0),CT[1],db[-1])
  names(out)=c("is.FA","DD","CT","LER","AAR")
  return(out)
}

