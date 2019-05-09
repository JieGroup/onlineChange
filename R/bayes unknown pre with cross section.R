detect.bayes.unknown.pre.cross=function(GEN,alpha,nulower=NULL,nuupper=NULL,score="hyvarinen",c=0.5,n=1000,
                         lenth1, lenth0, thlower1=NULL,thupper1=NULL, thlower0=NULL,thupper0=NULL,
                         GENTH1=NULL,ULPTH1=NULL,GENTH0=NULL,ULPTH0=NULL,
                         ULP0=NULL,GLP0=NULL,LLP0=NULL,ULP1=NULL,GLP1=NULL,LLP1=NULL,
                         par0=NULL,par1=NULL,
                         lenx=NULL,lower0=NULL,upper0=NULL,lower1=NULL,upper1=NULL)
{
  GEN=match.fun(GEN)
  if(alpha <=0 | alpha>=1) stop("'alpha' should be in (0,1)")
  score=match.arg(tolower(score),c("hyvarinen","logarithmic"))

  # suppose the time of changepoint nu follows a geometric prior distribution, compute the parameter p
  # the default is nulower=0, nuupper=18, which corresponds to p=0.1
  if(is.null(nulower)) nulower=0
  if(is.null(nuupper)) p=0.1 else{
    if(nulower >= nuupper) stop("need 'nulower' < 'nuupper'")
    p=1/(1+(nulower+nuupper)/2)
  }

  if(is.null(thlower1)) thlower1=rep(-Inf,lenth1)
  if(is.null(thupper1)) thupper1=rep(Inf,lenth1)

  if(is.null(thlower0)) thlower0=rep(-Inf,lenth0)
  if(is.null(thupper0)) thupper0=rep(Inf,lenth0)

  # create a matrix nuth that combines 3 parts by column:
  # a random sample of the time of changepoint nu, a vector of length n
  # a random sample of the transformed unknown parameters in the post-change model th1, a matrix with n rows and lenth1 columns
  # a random sample of the transformed unknown parameters in the pre-change model th0, a matrix with n rows and lenth0 columns

  # because we assume normal kernels for the unknown pre-change and post-change parameters in the MCMC procedure,
  # we transformed the bounded  parameters into unbounded ones here
  if(is.null(GENTH1) & is.null(ULPTH1)){

    # if not specified, the default assumption is that th1 follows N(0,1)
    ULPTH1=function(th1) -sum((th1-0)^2)/(2*1)
    GENTH1=function(n) matrix(rnorm(n*lenth1,mean=0,sd=1),nrow=n)
    nuth1=cbind(nulower+rgeom(n,p),GENTH1(n))
  } else{
    if(is.null(GENTH1)) stop("'GENTH1' is missing")
    if(is.null(ULPTH1)) stop("'ULPTH1' is missing")
    fULPTH1=match.fun(ULPTH1)
    ULPTH1=function(th1) f2(th1,thlower1,thupper1,fULPTH1)
    GENTH1=match.fun(GENTH1)
    nuth1=cbind(nulower+rgeom(n,p),m2(GENTH1(n),thlower1,thupper1,f2i))
  }
  if(is.null(GENTH0) & is.null(ULPTH0)){
    ULPTH0=function(th0) -sum(th0^2)/(2*1)
    GENTH0=function(n) matrix(rnorm(n*lenth0,sd=1),nrow=n)
    nuth=cbind(nuth1,GENTH0(n))
  } else{
    if(is.null(GENTH0)) stop("'GENTH0' is missing")
    if(is.null(ULPTH0)) stop("'ULPTH0' is missing")
    fULPTH0=match.fun(ULPTH0)
    ULPTH0=function(th0) f2(th0,thlower0,thupper0,fULPTH0)
    GENTH0=match.fun(GENTH0)
    nuth=cbind(nuth1,m2(GENTH0(n),thlower0,thupper0,f2i))
  }

  LPQ=function(nuth,numean,thmean1,thsd1,thmean0,thsd0){
    # sample nu from its prior distribution
    nu=nuth[1]
    # sample th1 and th0 from their prior distributions
    th1=nuth[2:(1+lenth1)]
    th0=nuth[(length(nuth)+1-lenth0) : length(nuth)]
    (nu-nulower)*(log(1-p)-log(1-1/(1+numean)))+ ULPTH1(th1)+sum((th1-thmean1)^2/thsd1^2)/2 +
      ULPTH0(th0)+sum((th0-thmean0)^2/thsd0^2)/2
  }
  if(score=="hyvarinen")
  {
    # some checks, preparations of the log probablity functions
    # to get the log probablity functions, we have to revert the transformation of the unknown parameters
    if(is.null(GLP0) | is.null(LLP0)){
      fULP0=match.fun(ULP0)
      ULP0=function(x,th) fULP0(x,v2(th,thlower0,thupper0,i2f))
    }
    if(is.null(GLP0)) GLP0=function(x,th) pracma::grad(ULP0,x,th=th) else{
      fGLP0=match.fun(GLP0)
      GLP0=function(x,th) fGLP0(x,v2(th,thlower0,thupper0,i2f))
    }
    if(is.null(LLP0)) LLP0=function(x,th) pracma::laplacian(ULP0,x,th=th) else{
      fLLP0=match.fun(LLP0)
      LLP0=function(x,th) fLLP0(x,v2(th,thlower0,thupper0,i2f))
    }

    if(is.null(GLP1) | is.null(LLP1)){
      fULP1=match.fun(ULP1)
      ULP1=function(x,th) fULP1(x,v2(th,thlower1,thupper1,i2f))
    }
    if(is.null(GLP1)) GLP1=function(x,th) pracma::grad(ULP1,x,th=th) else{
      fGLP1=match.fun(GLP1)
      GLP1=function(x,th) fGLP1(x,v2(th,thlower1,thupper1,i2f))
    }
    if(is.null(LLP1)) LLP1=function(x,th) pracma::laplacian(ULP1,x,th=th) else{
      fLLP1=match.fun(LLP1)
      LLP1=function(x,th) fLLP1(x,v2(th,thlower1,thupper1,i2f))
    }

    if(is.null(par0)) par0=1 else if(par0 <=0){
      warning("'par0' should be positive for hyvarinen score. Use par0=1")
      par0=1
    }
    if(is.null(par1)) par1=1 else if(par1 <=0){
      warning("'par1' should be positive for hyvarinen score. Use par1=1")
      par1=1
    }

    # compute the scores for the pre-change and post-change models
    SC0=function(x,th) par0*(sum(GLP0(x,th)^2)/2+LLP0(x,th))
    rSC0=function(th,x) SC0(x,th)
    hSC0=function(x,th){
      if(lenth0==1) sapply(th,rSC0,x) else apply(th,1,rSC0,x)
    }

    SC1=function(x,th) par1*(sum(GLP1(x,th)^2)/2+LLP1(x,th))
    rSC1=function(th,x) SC1(x,th)
    hSC1=function(x,th){
      if(lenth1==1) sapply(th,rSC1,x) else apply(th,1,rSC1,x)
    }

    # compute the log of the difference of the pre-change and post-change scores
    hLZSUM=function(nuth,t,x){
      out=numeric(n)
      nu=nuth[,1]
      id=which(nu < t)
      for(i in id){
        xtail=x
        if(nu[i]>0) xtail=x[-(1:nu[i])]
        out[i]=sum(sapply(xtail,SC0,nuth[i,(ncol(nuth)+1-lenth0):ncol(nuth)]))-
          sum(sapply(xtail,SC1,nuth[i,2:(1+lenth1)]))
      }
      return(out)
    }
  }else {
    # when choose logarithmic score
    fULP0=match.fun(ULP0)
    ULP0=function(x,th) fULP0(x,v2(th,thlower0,thupper0,i2f))

    fULP1=match.fun(ULP1)
    ULP1=function(x,th) fULP1(x,v2(th,thlower1,thupper1,i2f))

    if(is.null(par0) | is.null(par1)){
      if(is.null(lenx)) stop("'lenx' is missing")
      if(lenx <=0) stop("'lenx' should be a positive integer")
      if(is.null(lower0)) lower0=rep(-Inf,lenx)
      if(is.null(lower1)) lower1=rep(-Inf,lenx)
      if(is.null(upper0)) upper0=rep(Inf,lenx)
      if(is.null(upper1)) upper1=rep(Inf,lenx)
    }

    if(is.null(par1)){
      fn1=function(x,th) exp(ULP1(x,th))
      if(lenx==1) {par1=function(th) log(integrate(fn1,lower1,upper1,th)$value)
      }else par1=function(th) log(cubature::hcubature(fn1,lower1,upper1,th)$integral)
    } else{
      fpar1=match.fun(par1)
      par1=function(th) fpar1(v2(th,thlower1,thupper1,i2f))
    }

    if(is.null(par0)){
      fn0=function(x,th) exp(ULP0(x,th))
      if(lenx==1) {par0=function(th) log(integrate(fn0,lower0,upper0,th)$value)
      }else par0=function(th) log(cubature::hcubature(fn0,lower0,upper0,th)$integral)
    } else{
      fpar0=match.fun(par0)
      par0=function(th) fpar0(v2(th,thlower0,thupper0,i2f))
    }


    if(lenth1==1) lpar1=sapply(nuth[,2:(1+lenth1)],par1) else lpar1=apply(nuth[,2:(1+lenth1)],1,par1)
    if(lenth0==1) lpar0=sapply(nuth[,(ncol(nuth)+1-lenth0):ncol(nuth)],par0) else lpar0=apply(nuth[,(ncol(nuth)+1-lenth0):ncol(nuth)],1,par0)

    #if((sum(lpar1==(-Inf))+sum(lpar1==Inf))>0) stop("bad prior for post-change model,please check GENTH1 and ULPTH1")
    #if((sum(lpar0==(-Inf))+sum(lpar0==Inf))>0) stop("bad prior for pre-change model,please check GENTH0 and ULPTH0")

    SC1=function(x,th1,par1th) -ULP1(x,th1)+par1th
    SC0=function(x,th0,par0th) -ULP0(x,th0)+par0th

    lSC1=function(x,th,lpar1){
      out=numeric(length(lpar1))
      if(length(lpar1)==1){
        out=SC1(x,th,lpar1)
      }else{
        for(i in 1:length(lpar1)){
          if(lenth1==1) out[i]=SC1(x,th[i],lpar1[i]) else out[i]=SC1(x,th[i,],lpar1[i])
        }
      }
      return(out)
    }


    lSC0=function(x,th,lpar0){
      out=numeric(length(lpar0))
      if(length(lpar0)==1){
        out=SC0(x,th,lpar0)
      }else{
        for(i in 1:length(lpar0)){
          if(lenth0==1) out[i]=SC0(x,th[i],lpar0[i]) else out[i]=SC0(x,th[i,],lpar0[i])
        }
      }
      return(out)
    }

    lLZSUM=function(nuth,t,x,lpar1,lpar0){
      out=numeric(n)
      nu=nuth[,1]
      id=which(nu < t)
      for(i in id)
      {
        xtail=x
        if(nu[i]>0) xtail=x[-(1:nu[i])]
        out[i]=sum(sum(sapply(xtail,SC0,th0=nuth[i,(ncol(nuth)+1-lenth0):ncol(nuth)],par0th=lpar0[i]))-
                     sum(sapply(xtail,SC1,th1=nuth[i,2:(1+lenth1)],par1th=lpar1[i])))
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

  length.unique.id = 1000

  # At step t = 2, 3, ...:
  repeat
  {
    if(1/sum(w^2) < c*n)
    {

      # resample step:
      # use a multinomial distribution to obtain equally weighted samples
      id=mult(w,n)
      length.unique.id=length(unique(id))
      if ( length.unique.id==1) {
        warning("The multiple resampling is too concentrated. Please check the prior belief of both pre-change model and post-change model.")
        break
      }

      lowESS=lowESS+1
      # set the weights to be equal
      w=rep(1/n,n);lw=rep(-log(n),n)

      nuth=nuth[id,]
      if(score=="logarithmic") {
        lpar1=lpar1[id]
        lpar0=lpar0[id]
      }


      # move step:
      # consider an MCMC kernel targeting the joint posterior of (nu,th1,th0) conditional on the t-1 observations
      # compute the mean of randomly sampled nu
      numean=mean(nuth[,1])
      # draw a new nu from the MCMC kernel
      nu.p=nulower+rgeom(n,1/(1+(numean-nulower)))
      # draw new th1 and th0 from the MCMC kernel
      if(lenth1==1){
        thmean1=mean(nuth[,2:(1+lenth1)]);thsd1=sd(nuth[,2:(1+lenth1)])
        th.p1=rnorm(n)*thsd1+thmean1
      }else{
        thmean1=apply(nuth[,2:(1+lenth1)],2,mean);thsd1=apply(nuth[,2:(1+lenth1)],2,sd)
        th.p1=matrix(rnorm(n*lenth1)*rep(thsd1,n)+rep(thmean1,n),nrow=n,byrow=T)
      }


      if(lenth0==1){
        thmean0=mean(nuth[,(ncol(nuth)+1-lenth0):ncol(nuth)]);thsd0=sd(nuth[,(ncol(nuth)+1-lenth0):ncol(nuth)])
        th.p0=rnorm(n)*thsd0+thmean0
      }else{
        thmean0=apply(nuth[,(ncol(nuth)+1-lenth0):ncol(nuth)],2,mean);thsd0=apply(nuth[,(ncol(nuth)+1-lenth0):ncol(nuth)],2,sd)
        th.p0=matrix(rnorm(n*lenth0)*rep(thsd0,n)+rep(thmean0,n),nrow=n,byrow=T)
      }

      # combine the new nu and th1, th0 together
      nuth.p=cbind(nu.p,th.p1,th.p0)

      # compute the new normalizing constants for the logarithmic score
      if(score=="logarithmic"){
        if(lenth1==1) lpar1.p=sapply(th.p1,par1) else lpar1.p=apply(th.p1,1,par1)
        if(lenth0==1) lpar0.p=sapply(th.p0,par0) else lpar0.p=apply(th.p0,1,par0)
      }

      lrat=apply(nuth.p,1,LPQ,numean,thmean1,thsd1,thmean0,thsd0)-
        apply(nuth,1,LPQ,numean,thmean1,thsd1,thmean0,thsd0)

      if(score=="hyvarinen"){ lrat=lrat+hLZSUM(nuth.p,t,x)-hLZSUM(nuth,t,x)
      } else lrat=lrat+lLZSUM(nuth.p,t,x,lpar1.p,lpar0.p)-lLZSUM(nuth,t,x,lpar1,lpar0)

      # accept the new (nu,th1,th0) if the accept ratio smaller than a random sample from Unif(0,1)
      acc=(lrat >= log(runif(n)))
      if(sum(acc)>0){
        nuth[acc,]=nuth.p[acc,]
        if(score=="logarithmic") {
          lpar1[acc]=lpar1.p[acc]
          lpar0[acc]=lpar0.p[acc]
        }
      }
      ar=c(ar,mean(acc))
    }

    if (length.unique.id==1) t=unique(nuth[,1])
    if (length.unique.id==1) break

    t=t+1
    xt=GEN(t)
    if(is.na(xt)) {
      warning("Did not detect any change.")
      break
    }
    x=c(x,xt)
    xg=(nuth[,1] <t)

    if(sum(!xg)>0){
      if(score=="hyvarinen") {ngsc[!xg]=-hSC0(xt,nuth[!xg,(ncol(nuth)+1-lenth0):ncol(nuth)])
      }else ngsc[!xg]=-lSC0(xt,nuth[!xg,(ncol(nuth)+1-lenth0):ncol(nuth)],lpar0[!xg])
    }

    if(sum(xg)>0){
      if(score=="hyvarinen") {ngsc[xg]=-hSC1(xt,nuth[xg,2:(1+lenth1)])
      }else ngsc[xg]=-lSC1(xt,nuth[xg,2:(1+lenth1)],lpar1[xg])
    }

    lw.ngsc=lw+ngsc
    diff=lw.ngsc-max(lw.ngsc)
    # calclute the new weights w
    lw=diff-log(sum(exp(diff)))
    w=exp(lw)

    sum.w = sum(w[xg])
    if(is.na(sum.w)==TRUE) sum.w = 0

    if(t>nulower & sum.w >= 1-alpha) break
  }

  out=c(t,lowESS/t,mean(ar))
  names(out)=c("t","LER","AAR")

  out2 = nuth[,-1]

  mylist = list(out,out2)

  names(mylist) = c("detected change time","transformed cross section")
  return(mylist)
}
