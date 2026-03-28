#' @title LikMargIG
#'
#' @description This function calculates neg-loglikelihood for a proportional hazards model with inverse gaussian frailty.
#'
#' @param y Initial vector of parameters in the form:
#' \deqn{y = (ln(a), \beta , ln(\sigma ^2))}
#' for the exponential baseline hazard function,
#' \deqn{y = (ln(a), ln(b), \beta , ln(\sigma ^2))}
#' for the Weibull baseline hazard function and
#' \deqn{y = (ln(10^3a), ln(10b), \beta , ln(\sigma ^2))}
#' for the Gompertz baseline hazard function, where a and b are the slope and the shape parameters,
#' \eqn{\beta } are the Cox-regression parameters, and \eqn{\sigma ^2} is the variance of frailty.
#' @param D  A data.frame in which to interpret the variables named in the formula.
#' The data set includes the following fields:
#' \enumerate{
#' \item time-to-failure 'time';
#' \item censoring 'event' (censoring must be either 0 for no event or 1 for event);
#' \item  covariates (continuous or categorical, first 'nf' columns) used in a study (can be empty set);
#' \item cluster - the name of a cluster variable in data (is equal to NULL for the fixed-effect model).}
#' @param nf The number of covariates under study.
#' @param ncl The number of clusters.
#' @param dist Baseline hazard function ('Exponential', 'Weibull' or 'Gompertz').
#'
#' @details Three kinds of the baseline hazards are used in this function:
#' \enumerate{
#' \item Exponential with baseline cumulative hazard function \deqn{H_{0}(t;a)=at;}
#' \item Weibull with baseline cumulative hazard function \deqn{H_{0}(t;a,b)=at^b;}
#' \item Gompertz with baseline cumulative hazard function \deqn{H_{0}(t;a,b)=a(exp(bt)-1)).}
#'
#' @return Neg-loglikelihood
#'
#' @examples
#' \dontrun{
#' LikMargIG(y,D,nf,ncl,dist)
#' }
#' @export
#'
LikMargIG=function(y,D,nf,ncl,dist){
  #######This function calculates negloglikelihood for Weibull or Gompertz survival model
   k=1
   if (dist=='Weibull'){
    lambda0=exp(y[1])
    k0=exp(y[2])
    } else if (dist=='Gompertz') {
      lambda0=1e-4*exp(y[1])
      k0=1e-1*exp(y[2])
    } else if (dist=='Exponential'){
      lambda0=exp(y[1])
    }
  Cox=c(rep(0,nrow(D)))
  if (nf>0){
    for (i in 1:nf){
      if (dist=='Exponential'){
      Cox=Cox+D[,i]*y[1+i]
      } else {
      Cox=Cox+D[,i]*y[2+i]
      }
    }
  }
  LCox=Cox
  Cox=exp(Cox)
  if (ncl>0) {
    if (dist=='Exponential'){
      g=y[2+nf]
      G2=exp(y[2+nf])
    } else {
      g=y[3+nf]
    G2=exp(y[3+nf])
    }
    list=unique(D$cluster)
    nl=length(list)
  }
  Cens=D$event
  x1=D$time
  if (dist=='Exponential'){
    if (G2>1e-8 & ncl>0){
      Lik=0


      for (i in 1:nl){
        ID=list[i]
        ind=which(D$cluster==ID)
        nn=length(ind)
        ic1=1*(Cens[ind]==1)
        d=sum(ic1)
        Hfull1= Cox[ind]*lambda0*x1[ind]
        sH=sum(Hfull1)
        z=sqrt(1+2*G2*sH)/G2
        J=1+2*G2*sH
        mufull1=Cox[ind]*lambda0
        Lmufull1=log(lambda0)+LCox[ind]
        if (d>0) {
        as=log(besselK(z,d-0.5))
        Lik=Lik-(d-0.5)*log(z)+as-d*g+1/G2
        } else {
        Lik=Lik+(1-sqrt(1+2*G2*sH))/G2
       }
        Lik=Lik+sum(Lmufull1*ic1)
     }
    } else {
      ic=1*(Cens==1)
      Hfull1= Cox*lambda0*x1
      Lmufull1=LCox+log(lambda0)
      Lik=sum(Lmufull1*ic)-sum(Hfull1)
    }
  }

  if (dist=='Gompertz'){
    if (G2>1e-8 & ncl>0){
      Lik=0
      for (i in 1:nl){
        ID=list[i]
        ind=which(D$cluster==ID)
        nn=length(ind)
        ic1=1*(Cens[ind]==1)
        d=sum(ic1)
        Hfull1= Cox[ind]*(lambda0/k0)*(exp(k0*x1[ind])-1)
        sH=sum(Hfull1)
        z=sqrt(1+2*G2*sH)/G2
        J=1+2*G2*sH
        mufull1=Cox[ind]*lambda0*exp(k0*x1[ind])
        Lmufull1=LCox[ind]+log(lambda0)+k0*x1[ind]
        if (d>0) {
        as=log(besselK(z,d-0.5))
        Lik=Lik-(d-0.5)*log(z)+as-d*g+1/G2
        } else {
        Lik=Lik+(1-sqrt(1+2*G2*sH))/G2
      }
        Lik=Lik+sum(Lmufull1*ic1)
     }
    } else {
      ic=1*(Cens==1)
      Hfull1= Cox*(lambda0/k0)*(exp(k0*x1)-1)
      Lmufull1=LCox+log(lambda0)+k0*x1
      Lik=sum(Lmufull1*ic)-sum(Hfull1)

    }
  }
  if (dist=='Weibull'){
    if (G2>1e-8 & ncl>0){
      Lik=0
      for (i in 1:nl){
        ID=list[i]
        ind=which(D$cluster==ID)
        nn=length(ind)
        ic1=1*(Cens[ind]==1)
        d=sum(ic1)
        Hfull1= Cox[ind]*lambda0*x1[ind]^k0
        sH=sum(Hfull1)
        z=sqrt(1+2*G2*sH)/G2
        J=1+2*G2*sH
        mufull1=Cox[ind]*lambda0*k0*x1[ind]^(k0-1)
        Lmufull1=log(lambda0)+LCox[ind]+log(k0)+(k0-1)*log(x1[ind])
        if (d>0) {
        as=log(besselK(z,d-0.5))
        Lik=Lik-(d-0.5)*log(z)+as-d*g+1/G2
        } else {
        Lik=Lik+(1-sqrt(1+2*G2*sH))/G2
       }
        Lik=Lik+sum(Lmufull1*ic1)
     }
    } else {
      ic=1*(Cens==1)
      Hfull1= Cox*lambda0*x1^k0
      Lmufull1=LCox+log(lambda0)+log(k0)+(k0-1)*log(x1)
      Lik=sum(Lmufull1*ic)-sum(Hfull1)
    }
  }
  Lik=-Lik
  if (is.na(Lik) | is.infinite(Lik)) Lik=1e+50
  return(Lik)
}

