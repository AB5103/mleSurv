#' @title BiasCorrectionGamma
#'
#' @description   Returns the bias correction for ML estimates in the model with gamma frailty.
#'
#' @param D The data frame. The first 'nf' columns are the covariates under study. Other columns must include fields
#' 'time' for time-to-failure, 'Cens' for censoring and 'cluster' for cluster identification. The number of rows is equal to the number of subjects.
#' @param nf The number of covariates under study.
#' @param ncl The number of clusters.
#' @param para  The vector of the ML estimates of unknown parameters. This vector has a form
#' \deqn{y = (ln(a), \beta , ln(\sigma ^2))}
#' for the exponential baseline hazard function,
#' \deqn{y = (ln(a), ln(b), \beta , ln(\sigma ^2))}
#' for the Weibull baseline hazard function and
#' \deqn{y = (ln(10^3a), ln(10b), \beta , ln(\sigma ^2))}
#' for the Gompertz baseline hazard function, where a and b are the slope and the shape parameters,
#' \eqn{\beta } are the Cox-regression parameters, and \eqn{\sigma ^2} is the variance of frailty.
#' @param nMonte  The number of the resamplings for the Monte Carlo integrations
#' 
#' @details This function calculates the bias correction for the proportional hazards model with gamma frailty
#' in accordance with the Cox-Snell method.  The baseline cumulative hazard functions are:
#' \enumerate{
#' \item Exponential with baseline cumulative hazard function \deqn{H_{0}(t;a)=at;}
#' \item Weibull with baseline cumulative hazard function \deqn{H_{0}(t;a,b)=at^b;}
#' \item Gompertz with baseline cumulative hazard function \deqn{H_{0}(t;a,b)=a(exp(bt)-1)).}
#' }
#'
#' @return Two fields - the bias correction \eqn{Bias} and its history \eqn{Biasdiag} for ML estimates
#' (for \eqn{nMonte} replacements) in the form:
#' \deqn{y1 = (ln(a), \beta , ln(\sigma ^2))}
#' for the exponential baseline hazard function,
#' \deqn{y1 = (ln(a), \beta , ln(\sigma ^2), b)}
#' for the Weibull baseline hazard function and
#' \deqn{y1 = (ln(10^3a), \beta , ln(\sigma ^2), 10b)}

#' If the calculation of the bias correction is not possible returns NAs.
#'
#' @import MASS
#'
#' @examples
#' \dontrun{
#' BiasCorrectionGamma(D,nf,ncl,para,nMonte)
#'
#'}
#' @import MASS
#' @export
#'

BiasCorrectionGamma=function(D0,nf,ncl,dist,para,nMonte){
 if (nf>0 & dist=="Exponential")  parnam=c("a",colnames(D0)[1:nf],"g")
 if (nf>0 & dist!="Exponential")  parnam=c("a",colnames(D0)[1:nf],"g","b")
 if (nf==0 & dist=="Exponential")  parnam=c("a","g")
 if (nf==0 & dist!="Exponential")  parnam=c("a","g","b")
 nr=nrow(D0)
 ord=c(1:nr)
 nfg=nf+2

 if (dist=='Exponential'){
 np=nf+2
 np1=np 
 a=para[1]
 g=para[nf+2]
 if (nf!=0) {par=c(a,para[2:(1+nf)],g)} else {par=c(a,g)}
 } else if (dist=='Weibull') {
 np=nf+3
 np1=nf+2 
 a=para[1]
 b=exp(para[2])
 g=para[nf+3]
 if (nf!=0) {par=c(a,para[3:(2+nf)],g,b)} else {par=c(a,g,b)}
 } else if (dist=='Gompertz') {
 np=nf+3
 np1=nf+2 
 a=para[1]-para[2]+log(1e-3)
 b=0.1*exp(para[2])
 g=para[nf+3]
 if (nf!=0) {par=c(a,para[3:(2+nf)],g,b)} else {par=c(a,g,b)}
 } 
 G2=exp(g)
 RC=  array(0,c(np))
 RCC=  array(0,c(np,np))
 RCCC= array(0,c(np,np,np))
 RCC_C=array(0,c(np,np,np))
 Biasdiag={}
 iMonte=0 
  while (iMonte<=nMonte) {
  iMonte=iMonte+1
  iord=sort(sample(ord,nr,replace=TRUE))
  D=D0[iord,]
  a.scale=rep(1,nr)
  g.scale=rep(1,nr)
  if (nf!=0) {uR_beta=D[,1:nf]} else {uR_beta={}}
  UR=cbind(a.scale,uR_beta,g.scale)

  Cens=D$event
  Tcens=c(D$time)
  pow=a*c(UR[,1]) 
  for (i in 2:np1){
  pow=  pow+  par[i]*UR[,i]
   }
  if (dist=='Exponential') {
  HR=exp(pow)*Tcens}                        else if (dist=='Gompertz') {
  HR=exp(pow)*(exp(b*Tcens)-1)
  URb=Tcens*exp(b*Tcens)/(exp(b*Tcens)-1)
  URb_b=Tcens*URb-URb^2
  URb_bb=Tcens^2*URb-3*Tcens*URb^2+2*URb^3
  URbb=Tcens^2*exp(b*Tcens)/(exp(b*Tcens)-1) 
  URbbb=Tcens^3*exp(b*Tcens)/(exp(b*Tcens)-1)} else if (dist=='Weibull') {

  HR=exp(pow)*Tcens^b
  }
  if (dist=='Gompertz')    UR=cbind(UR,URb) 
  if (dist=='Weibull')     UR=cbind(UR,log(Tcens))
  list=sort(unique(D$cluster))
################################################# begin cycle by clusters
  for (lj in 1:ncl){
    ind=which(D$cluster==list[lj])
    n=length(ind)
####################################################
    rC=array(0,c(np))
    rCC=array(0,c(np,np))
#################################################### begin check for zero sample in a cluster
    if (n!=0) {
#################################################### begin preliminaries
    Censd=Cens[ind]
    d=sum(Censd)
    if (dist=='Gompertz') {
    Ub=matrix(as.numeric(URb[ind]),length(ind),1)
    Ub_b=matrix(as.numeric(URb_b[ind]),length(ind),1)
    Ub_bb=matrix(as.numeric(URb_bb[ind]),length(ind),1)

    Ubb=matrix(as.numeric(URbb[ind]),length(ind),1)
    Ubbb=matrix(as.numeric(URbbb[ind]),length(ind),1)
    }
    if (n==1)  {
    U=matrix(as.numeric(UR[ind,]),length(ind),ncol(UR))
    H=matrix(as.numeric(HR[ind]),length(ind),1)
    Tc=matrix(as.numeric(Tcens[ind]),length(ind),1)} else {
    U=UR[ind,]
    H=HR[ind]
    Tc=Tcens[ind]
    }
    sH=sum(H)
    nul=0
    uH=c(rep(NA,np))
    uuH=matrix(NA,np,np)
    uuH=matrix(NA,np,np)
    uuuH=array(NA,c(np,np,np))
    for (j1 in 1:np){
      uH[j1]            =sum(c(U[,j1])*H)
      for (j2 in 1:np){
        if (dist!='Gompertz' || (dist=='Gompertz' && (j1!=np || j2!=np))) uuH[j1,j2]=sum(c(U[,j1]*U[,j2])*H)
        if (dist=='Gompertz' && j1==np && j2==np)   uuH[np,np]=sum(c(Ubb)*H)
        for (j3 in 1:np){
        if (dist!='Gompertz' || (dist=='Gompertz' && ((j1!=np && j2!=np)|| (j1!=np && j3!=np) || (j2!=np && j3!=np)))) uuuH[j1,j2,j3]=sum(c(U[,j1]*U[,j2]*U[,j3])*H)
        if (dist=='Gompertz' && ((j1==np && j2==np && j3!=np)))    uuuH[j1,j2,j3]=sum(c(Ubb)*c(U[,j3])*H)
        if (dist=='Gompertz' && ((j1==np && j2!=np && j3==np)))    uuuH[j1,j2,j3]=sum(c(Ubb)*c(U[,j2])*H)
        if (dist=='Gompertz' && ((j1!=np && j2==np && j3==np)))    uuuH[j1,j2,j3]=sum(c(Ubb)*c(U[,j1])*H)
        if (dist=='Gompertz' && ((j1==np && j2==np && j3==np)))    uuuH[j1,j2,j3]=sum(c(Ubbb)*H)
        }
      }
    }
    Part3=expression(lgamma(exp(-g)+d)-lgamma(exp(-g)))
    Part3_1=as.expression(DD(Part3,"g",1))
    Part3_2=as.expression(DD(Part3,"g",2))
    Part3_3=as.expression(DD(Part3,"g",3))
##################################################    end preliminaries
    #First, second and third order derivatives F1+F2+F3
################################################## begin F1+F2+F3
    #First, second and third order derivatives
################################################## begin F1

    if (dist=='Weibull') {
    rC[np]=rC[np]+d/b
    }
    if (dist=='Gompertz') {
    rC[np]=rC[np]+d/b+sum(Tc*Censd)
    }
    for (j1 in 1:np){
      if (dist=='Exponential' || dist=='Weibull') {
      rC[j1] =  rC[j1]   +sum(U[,j1]*Censd)
      }
      if (dist=='Gompertz' && j1!=np) {
      rC[j1] = rC[j1]    +sum(U[,j1]*Censd)
      }
    }

    rC[np1]=        rC[np1]    +eval(Part3_1)+exp(-g)*log(1+sH)
    for (i in 1:np){
    rC[i]=          rC[i]   -(exp(-g)+d)*uH[i]/(1+sH)
    }
    ########################################## end F1
    #Second order derivative
    ########################################## begin F2
    if (dist=='Weibull' || dist=='Gompertz') {
    rCC[np,np] =   rCC[np,np]    -d/b^2
    }
    rCC[np1,np1]=    rCC[np1,np1]    +eval(Part3_2)-exp(-g)*log(1+sH)
    ###########################F2.1
    for (i in 1:np){
    rCC[i,np1]=    rCC[i,np1]   +exp(-g)*uH[i]/(1+sH)
    rCC[np1,i]=    rCC[np1,i]   +exp(-g)*uH[i]/(1+sH)
    }
    ###########################F2.2
    for (j1 in 1:np){
      for (j2 in 1:np){
      rCC[j1,j2]=    rCC[j1,j2]    -(exp(-g)+d)*uuH[j1,j2]/(1+sH)
      } 
    }
    ###########################F2.3
    for (j1 in 1:np){
      for (j2 in 1:np){
      rCC[j1,j2] =   rCC[j1,j2]    +(exp(-g)+d)*uH[j1]*uH[j2]/(1+sH)^2
      } 
    }

######################################
    for (j1 in 1:np){
      RC[j1]    =  RC[j1]+rC[j1]
      for (j2 in 1:np){
      RCC[j1,j2] =  RCC[j1,j2]+rCC[j1,j2]
        for (j3 in 1:np){
        RCC_C[j1,j2,j3]=RCC_C[j1,j2,j3]+rCC[j1,j2]*rC[j3]
        }
      } 
    }
###################################### end F2
    #Third order derivative
    ########################################## begin F3
    ###########################F3.1
    if (dist=='Weibull' || dist=='Gompertz') {
    RCCC[np,np,np]= RCCC[np,np,np]+2*d/b^3
    }
    RCCC[np1,np1,np1]=RCCC[np1,np1,np1]+eval(Part3_3)+exp(-g)*log(1+sum(H))
    for (i in 1:np){
    RCCC[i,np1,np1]=RCCC[i,np1,np1]-exp(-g)*uH[i]/(1+sH)
    RCCC[np1,i,np1]=RCCC[np1,i,np1]-exp(-g)*uH[i]/(1+sH)
    RCCC[np1,np1,i]=RCCC[np1,np1,i]-exp(-g)*uH[i]/(1+sH)
    }

    ###########################F3.2
    for (j1 in 1:np){
      for (j2 in 1:np){
      RCCC[np1,j1,j2]=RCCC[np1,j1,j2]+exp(-g)*uuH[j1,j2]/(1+sH)
      RCCC[j1,np1,j2]=RCCC[j1,np1,j2]+exp(-g)*uuH[j1,j2]/(1+sH)
      RCCC[j1,j2,np1]=RCCC[j1,j2,np1]+exp(-g)*uuH[j1,j2]/(1+sH)
      } 
    }
    ###########################F3.3
    for (j1 in 1:np){
      for (j2 in 1:np){
      RCCC[np1,j1,j2]=RCCC[np1,j1,j2]-exp(-g)*uH[j1]*uH[j2]/(1+sH)^2
      RCCC[j1,np1,j2]=RCCC[j1,np1,j2]-exp(-g)*uH[j1]*uH[j2]/(1+sH)^2
      RCCC[j1,j2,np1]=RCCC[j1,j2,np1]-exp(-g)*uH[j1]*uH[j2]/(1+sH)^2
      } 
    }
    ###########################F3.4
    for (j1 in 1:np){
      for (j2 in 1:np){
         for (j3 in 1:np){
         RCCC[j1,j2,j3]=RCCC[j1,j2,j3]-(exp(-g)+d)*uuuH[j1,j2,j3]/(1+sH)
         }
      } 
    }
    ####################################F3.5
    for (j1 in 1:np){
      for (j2 in 1:np){
         for (j3 in 1:np){
         RCCC[j1,j2,j3]=RCCC[j1,j2,j3]+(exp(-g)+d)*uuH[j1,j2]*uH[j3]/(1+sH)^2
         RCCC[j1,j2,j3]=RCCC[j1,j2,j3]+(exp(-g)+d)*uuH[j1,j3]*uH[j2]/(1+sH)^2
         RCCC[j1,j2,j3]=RCCC[j1,j2,j3]+(exp(-g)+d)*uuH[j2,j3]*uH[j1]/(1+sH)^2
         }
      } 
    }
    ####################################F3.6
    for (j1 in 1:np){
      for (j2 in 1:np){
         for (j3 in 1:np){
         RCCC[j1,j2,j3]=RCCC[j1,j2,j3]-2*(exp(-g)+d)*uH[j1]*uH[j2]*uH[j3]/(1+sH)^3
         }
      } 
    }
########################################################### end F3
  } #n!=0
#################################################### end check for zero sample in a cluster
  } #lj
################################################# end cycle by clusters
  RCi=      RC/iMonte
  RCCi=    RCC/iMonte
  RCCCi=  RCCC/iMonte
  RCC_Ci=RCC_C/iMonte

  if (!(any(is.na(RCCi)) | any(is.na(RCCCi)))){
  RCCi_=ginv(RCCi)
  krai=0.5*RCCCi+RCC_Ci
  Biasi={}
  for (ir in 1:np){
    worki=0
    for (ir1 in 1:np){
      for (ir2 in 1:np){
        for (ir3 in 1:np){
         worki=worki+krai[ir1,ir2,ir3]*RCCi_[ir,ir1]*RCCi_[ir2,ir3]
        }
      }
    }
    Biasi=c(Biasi,worki)
 }
  } else {
    Biasi=rep(NA,np)
    stop("Bias cannot be calculated")
  }
  Biasdiag=rbind(Biasdiag,Biasi)
}#iMonte
BiasR=NULL
BiasR$Bias=Biasi
BiasR$diag=Biasdiag
return(BiasR)
}
