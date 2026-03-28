#' @title BiasCorrectionIG
#'
#' @description   Returns the bias correction for ML estimates in the model with inverse gaussian frailty.
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
#' BiasCorrectionIG(D,nf,ncl,para,nMonte)
#'
#'}
#' @import MASS
#' @export
#'
BiasCorrectionIG=function(D0,nf,ncl,dist,para,nMonte){ 
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
 RC=array(0,c(np))
 RCC=array(0,c(np,np))
 RCCC=array(0,c(np,np,np))
 RCC_C=array(0,c(np,np,np))

 RC1=array(0,c(np))
 RCC1=array(0,c(np,np))
 RC2=array(0,c(np))
 RCC2=array(0,c(np,np))
 RC3=array(0,c(np))
 RCC3=array(0,c(np,np))

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
  if (dist=='Exponential') {
    if (nf!=0) {
    parnam=c("a",paste0("beta_",as.character(1:nf)),"g")} else {
    parnam=c("a","g")}
  } else {
    if (nf!=0) {
    parnam=c("a",paste0("beta_",as.character(1:nf)),"g","b")} else {
    parnam=c("a","g","b")}
    }
  list=sort(unique(D$cluster))
  for (lj in 1:ncl){
    ind=which(D$cluster==list[lj])
    n=length(ind)
##################################################
    if (n!=0) {
################################################## begin preliminaries
    rC=array(0,c(np))
    rCC=array(0,c(np,np))
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
        if (dist!='Gompertz' || (dist=='Gompertz' && (j1!=np || j2!=np))) uuH[j1,j2]      =sum(c(U[,j1]*U[,j2])*H)
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
##################################################    end preliminaries
    #First, second and third order derivatives
################################################## begin F2
    if (dist=='Weibull')  {
    rC[np]=rC[np]+d/b
    }
    if (dist=='Gompertz') {
    rC[np]=rC[np]+d/b+sum(Tc*Censd)
    }
    for (j1 in 1:np){
      if ((dist=='Exponential' || dist=='Weibull') && j1!=np1) {
      rC[j1] =  rC[j1]   +sum(U[,j1]*Censd)
      }
      if (dist=='Gompertz' && j1!=np && j1!=np1) {
      rC[j1] = rC[j1]    +sum(U[,j1]*Censd)
      }
    }
    if (dist=='Weibull' || dist=='Gompertz') {
    rCC[np,np]    = -d/b^2
    }
################################################## end F2
    J=2*sH+1
    J05=J^0.5
    J_05=J^(-0.5)
    J_15=J^(-1.5)
    J_25=J^(-2.5)

    z=J^0.5/G2
    z_g=-z
    z_gg=z
    z_ggg=-z
    Kq_05z=besselK(z,d-0.5) #integrate(dbesl,0.001,100,nu=d-0.5,w=z,rel.tol = .Machine$double.eps^0.88)$value
    Kq_15z=besselK(z,d-1.5) #integrate(dbesl,0.001,100,nu=d-1.5,w=z,rel.tol = .Machine$double.eps^0.88)$value
    Kq_25z=besselK(z,d-2.5) #integrate(dbesl,0.001,100,nu=d-2.5,w=z,rel.tol = .Machine$double.eps^0.88)$value
    Kq_35z=besselK(z,d-3.5) #integrate(dbesl,0.001,100,nu=d-3.5,w=z,rel.tol = .Machine$double.eps^0.88)$value
    Kq05z =besselK(z,d+0.5) #integrate(dbesl,0.001,100,nu=d+0.5,w=z,rel.tol = .Machine$double.eps^0.88)$value
    Kq15z =besselK(z,d+1.5) #integrate(dbesl,0.001,100,nu=d+1.5,w=z,rel.tol = .Machine$double.eps^0.88)$value
    Kq25z =besselK(z,d+2.5) #integrate(dbesl,0.001,100,nu=d+2.5,w=z,rel.tol = .Machine$double.eps^0.88)$value
    K0=Kq_05z
    K1=-0.5*(Kq_15z+Kq05z)
    K2=0.25*(Kq_25z+2*Kq_05z+Kq15z)
    K3=-0.125*(Kq_35z+3*Kq_15z+3*Kq05z+Kq25z)
################################################## begin derivative of J
    rC[np1]=          rC[np1] -0.5     -1/G2
    rCC[np1,np1]=     rCC[np1,np1]     +1/G2

    RCCC[np1,np1,np1]=RCCC[np1,np1,np1]-1/G2
    for (j1 in 1:np) {
    J_p1=2*uH[j1]
    lnJ_p1=J_p1/J
    rC[j1]=rC[j1]-0.5*(d-0.5)*lnJ_p1
       for (j2 in 1:np) {
       J_p2=  2*uH[j2]
       J_p1p2=2*uuH[j1,j2]
       lnJ_p2=J_p2/J
       lnJ_p1p2=J_p1p2/J-(J_p1/J)*(J_p2/J)
       rCC[j1,j2]=rCC[j1,j2]-0.5*(d-0.5)*lnJ_p1p2
          for (j3 in 1:np) {
          J_p3=    2*uH[j3]
          J_p1p3=  2*uuH[j1,j3]
          J_p2p3=  2*uuH[j2,j3]
          J_p1p2p3=2*uuuH[j1,j2,j3]
          lnJ_p1p2p3=J_p1p2p3/J-(J_p1p2/J)*(J_p3/J)-(J_p1p3/J)*(J_p2/J)-(J_p2p3/J)*(J_p1/J)+2*(J_p1/J)*(J_p2/J)*(J_p3/J)  
          RCCC[j1,j2,j3]=RCCC[j1,j2,j3]-0.5*(d-0.5)*lnJ_p1p2p3
          }
       }
    }
##################################################   end derivative of J
##################################################   begin derivative of K
    
    K0_g=K1*z_g                                  
    K0_gg=K2*z_g^2+K1*z_gg                       
    K0_ggg=K3*z_g^3+3*K2*z_g*z_gg+K1*z_ggg       
    lnK_g=K0_g/K0
    lnK_gg=K0_gg/K0-(K0_g/K0)^2
    lnK_ggg=K0_ggg/K0-3*(K0_g/K0)*(K0_gg/K0)+2*(K0_g/K0)^3
    rC[np1]=rC[np1]+lnK_g
    rCC[np1,np1]=rCC[np1,np1]+lnK_gg
    RCCC[np1,np1,np1]=RCCC[np1,np1,np1]+lnK_ggg
    for (j1 in 1:np) {
    z_p1=uH[j1]/(z*G2^2)   
    z_gp1=-z_p1
    z_ggp1=z_p1
    K0_p1=K1*z_p1
    K0_gp1=K2*z_g*z_p1+K1*z_gp1                                                                       
    K0_ggp1=K3*z_g^2*z_p1+K2*z_gg*z_p1+2*K2*z_g*z_gp1+K1*z_ggp1                                       
    lnK_p1=K0_p1/K0
    lnK_gp1=K0_gp1/K0-(K0_p1/K0)*(K0_g/K0)
    lnK_ggp1=K0_ggp1/K0-(K0_gg/K0)*(K0_p1/K0)-2*(K0_g/K0)*(K0_gp1/K0)+2*(K0_g/K0)^2*(K0_p1/K0)#K0_ggp1/K0-(K0_gp1/K0)*(K0_g/K0)-(K0_gp1/K0-(K0_p1/K0)*(K0_g/K0))*(K0_g/K0)-(K0_p1/K0)*(K0_gg/K0-(K0_g/K0)^2)
    rC[j1]=rC[j1]+lnK_p1
    rCC[j1,np1]=rCC[j1,np1]+lnK_gp1
    rCC[np1,j1]=rCC[np1,j1]+lnK_gp1
    RCCC[j1,np1,np1]=RCCC[j1,np1,np1]+lnK_ggp1
    RCCC[np1,j1,np1]=RCCC[np1,j1,np1]+lnK_ggp1
    RCCC[np1,np1,j1]=RCCC[np1,np1,j1]+lnK_ggp1
        for (j2 in 1:np) {
        z_p2=uH[j2]/(z*G2^2)   
        z_gp2=-z_p2
        z_ggp2=z_p2
        z_p1p2=uuH[j1,j2]/(z*G2^2)-uH[j1]*uH[j2]/(z^3*G2^4)
        z_gp1p2=-z_p1p2
        K0_p2=K1*z_p2
        K0_p1p2=K2*z_p1*z_p2+K1*z_p1p2
        K0_gp2=K2*z_g*z_p2+K1*z_gp2
        K0_ggp2=K3*z_g^2*z_p2+K2*z_gg*z_p2+2*K2*z_g*z_gp2+K1*z_ggp2
        K0_gp1p2=K3*z_g*z_p1*z_p2+K2*(z_gp1*z_p2+z_gp2*z_p1+z_g*z_p1p2)+K1*z_gp1p2
        lnK_p2=K0_p2/K0
        lnK_gp2=K0_gp2/K0-(K0_p2/K0)*(K0_g/K0)
        lnK_p1p2=K0_p1p2/K0-(K0_p1/K0)*(K0_p2/K0)       
        lnK_gp1p2=K0_gp1p2/K0-(K0_p1p2/K0)*(K0_g/K0)-(K0_gp1/K0)*(K0_p2/K0)-(K0_gp2/K0)*(K0_p1/K0)+2*(K0_g/K0)*(K0_p1/K0)*(K0_p2/K0)
        rCC[j1,j2]=rCC[j1,j2]+lnK_p1p2
        RCCC[np1,j1,j2]=RCCC[np1,j1,j2]+lnK_gp1p2
        RCCC[j1,np1,j2]=RCCC[j1,np1,j2]+lnK_gp1p2
        RCCC[j1,j2,np1]=RCCC[j1,j2,np1]+lnK_gp1p2

           for (j3 in 1:np) {
           z_p3=uH[j3]/(z*G2^2)   
           z_p1p3=uuH[j1,j3]/(z*G2^2)-uH[j1]*uH[j3]/(z^3*G2^4)
           z_p2p3=uuH[j2,j3]/(z*G2^2)-uH[j2]*uH[j3]/(z^3*G2^4)
           z_p1p2p3=uuuH[j1,j2,j3]/(z*G2^2)-(uuH[j1,j2]*uH[j3]+uuH[j1,j3]*uH[j2]+uuH[j2,j3]*uH[j1])/(z^3*G2^4)+3*uH[j1]*uH[j2]*uH[j3]/(z^5*G2^6)
           K0_p3=K1*z_p3               
           K0_p1p3=K2*z_p1*z_p3+K1*z_p1p3                    
           K0_p2p3=K2*z_p2*z_p3+K1*z_p2p3                    
           K0_p1p2p3=K3*z_p1*z_p2*z_p3+K2*(z_p1p2*z_p3+z_p1p3*z_p2+z_p2p3*z_p1)+K1*z_p1p2p3
           lnK_p1p2p3=K0_p1p2p3/K0-(K0_p1p2/K0)*(K0_p3/K0)-(K0_p1p3/K0)*(K0_p2/K0)-(K0_p2p3/K0)*(K0_p1/K0)+2*(K0_p1/K0)*(K0_p2/K0)*(K0_p3/K0)
           RCCC[j1,j2,j3]=RCCC[j1,j2,j3]+lnK_p1p2p3
           }
        }
    }
#################################################### end derivative of K
    for (j1 in 1:np){
      RC[j1]    =  RC[j1]+rC[j1]
      for (j2 in 1:np){
      RCC[j1,j2] =  RCC[j1,j2]+rCC[j1,j2]
        for (j3 in 1:np){
        RCC_C[j1,j2,j3]=RCC_C[j1,j2,j3]+rCC[j1,j2]*rC[j3]
        }
      } 
    }
##################################################   end end
    if (dist=='Weibull' || dist=='Gompertz') {
    RCCC[np,np,np]= RCCC[np,np,np]+2*d/b^3
    }
##################################################
  }#n!=0
#################################################
  }#lj
#################################################
################################## begin diagnostic
  RCi=      RC/iMonte
  RCCi=RCC/iMonte
  RCC_Ci=RCC_C/iMonte
  RCCCi=RCCC/iMonte
  if (!(any(is.na(RCCi)) || any(is.na(RCCCi)))){
  RCCi_=ginv(RCCi)
  krai=RCC_Ci+0.5*RCCCi
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
    print("Bias cannot be calculated")
  }
  Biasi
  Biasdiag=rbind(Biasdiag,Biasi)
 }#iMonte

BiasR=NULL
BiasR$Bias=Biasi
BiasR$diag=Biasdiag
return(BiasR)
}

