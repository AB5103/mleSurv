#' @title opt.model
#'
#' @description   Calculates AIC and BIC for ML estimator applied to the proportional hazards model with gamma frailty
#' and different hazard functions.
#'
#' @param formula A formula object, with the response on the left of a ~ operator, and the terms on
#' the right. The response must be a survival object as returned by the Surv() function.
#' The status indicator corresponding to censoring in the Surv object must be 0 for censored objects and 1 for non-censored ones.
#' This object describes the effect of several factors in the Cox proportional hazards model.
#' @param data The data frame. Includes fields for:
#' \enumerate{
#' \item Time-to-failure;
#' \item Censoring;
#' \item Possible covariates (optional);
#' \item Cluster (if any).
#' }
#' The names of these fields must be included in 'formula'.
#' @param cluster The name of a cluster variable in data (is equal to NULL for the fixed-effect model).
#'
#' @return Retutns the latex table with AIC and BIC values for three hazard functions: exponential, Gompertz and Weibull.
#'
#' @import survival
#' @import MASS
#' @import ucminf
#' @import stats
#' @import xtable
#' @export
#'
#' @examples
#' \dontrun{
#' D=veteran
#' D$prior <- factor(as.character(D$prior), labels = c(0, 1))
#' D$trt <- factor(as.character(D$trt), labels = c(0, 1))
#' cluster="cluster"
#' D=data.frame(D$karno,D$trt,D$diagtime,D$prior,D$age,D$celltype,D$time,D$status)
#' colnames(D)=c("karno","trt","dtime","prior","age","cluster","time","status")
#' formula=as.formula("Surv(time, status) ~ karno")
#' RRR=opt.model(formula,D,cluster)
#' }
#'
#' @export
#'
opt.model=function(formula,data,cluster)
{
  nr=nrow(data)
  obsdata <- NULL
  obsdata$time <- eval(formula[[2]][[2]], envir = data)
  obsdata$event <- eval(formula[[2]][[3]], envir = data)
  obsdata$x <- as.data.frame(model.matrix.lm(formula, data = data,na.action='na.pass'))
  names(obsdata$x)=names(obsdata$x)

  ind.x=which(is.na(c(apply(obsdata$x,1,sum))))

  if (is.null(cluster)) {
    obsdata$ncl <- 0
    obsdata$di <- sum(obsdata$event)
    obsdata$cluster <- c(rep(0,nrow(data)))
    ind.cl <- as.numeric({})
  }   else {
    if (!cluster %in% names(data)) {
      stop(paste0("object '", cluster, "' not found"))
    }
    obsdata$cluster <- eval(as.name(cluster), envir = data)
    obsdata$ncl <- length(levels(as.factor(obsdata$cluster)))
    obsdata$di <- aggregate(obsdata$event, by = list(obsdata$cluster),FUN = sum)[, , drop = FALSE]
    cnames <- obsdata$di[, 1]
    obsdata$di <- as.vector(obsdata$di[, 2])
    names(obsdata$di) <- cnames
    ind.cl=which(is.na(obsdata$cl))
  }
  ind=sort(unique(c(ind.x,ind.cl)))
  nx=rep(0,nr)
  nx[ind]=1
  if (is.factor(obsdata$cluster)) obsdata$cluster=as.character(obsdata$cluster)
  obs=data.frame(obsdata$x[-1],obsdata$event,obsdata$time,obsdata$cluster)[nx!=1,]
  nr=nrow(obs)
  namesf=names(obsdata$x)[-1]
  names(obs)=c(names(obsdata$x)[-1],names(obsdata$event),"event","time","cluster")
  ncl=obsdata$ncl
  nf=length(namesf)
  Dist=c('Exponential','Gompertz','Weibull')
  Frailty=c('gamma','inverse gaussian')
  BBIC={}
  for (j in 1:2) {
  ABIC={}
  frailty=Frailty[j]
  for (i in 1:3){
    dist=Dist[i]
  par0=c(0,0)
  if (nf>0)  par0=c(0,c(rep(0,nf)),0)
  if (!(dist=='Exponential')){
    if (dist=='Weibull') {
    par0=c(0,0,1)
    if (nf>0)  par0=c(0,0,c(rep(0,nf)),0)} else {
    par0=c(0,-10,1)
    if (nf>0)  par0=c(0,-10,c(rep(0,nf)),0)
    }
  }
  if (frailty=='gamma')            ResultMarg=ucminf(par=par0,fn=LikMargGamma,gr=NULL,D=obs,nf=nf,ncl=ncl,dist=dist,hessian=3)
  if (frailty=='inverse gaussian') ResultMarg=ucminf(par=par0,fn=LikMargIG,gr=NULL,D=obs,nf=nf,ncl=ncl,dist=dist,hessian=3)
  para=ResultMarg$par
  if (any(!is.finite(as.matrix(ResultMarg$hessian))) | any(suppressWarnings(diag(ginv(ResultMarg$hessian)))<0)) {
    AIC=NA
    BIC=NA
  } else {
    Lik=-ResultMarg$value
    AIC=2*(length(para)-Lik)
    BIC=length(para)*log(nrow(obs))-2*Lik
  }
  ABIC=rbind(ABIC,c(AIC,BIC))
  }
  BBIC=cbind(BBIC,ABIC)
  }

  
  BBIC=data.frame(Dist,BBIC)
  rownames(BBIC)=c("Exponential","Gompertz","Weibull")
  colnames(BBIC)=c("Distance","AIC gamma","BIC gamma","AIC inverse gaussian","BIC inverse gaussian")
  xtab=xtable(BBIC)
  print(xtable(xtab, caption = "AIC and BIC values"),include.rownames=FALSE)
    return(xtab)
}
