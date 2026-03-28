library(mlesurv)
library(ucminf)
library(MASS)
library(xtable)
library(stats)
library(graphics)
library(survival)

 D=veteran
 D$prior <- factor(as.character(D$prior), labels = c(0, 1))
 D$trt <- factor(as.character(D$trt), labels = c(0, 1))
 cluster="cluster"
 D=data.frame(D$karno,D$trt,D$diagtime,D$prior,D$age,D$celltype,D$time,D$status)
 colnames(D)=c("karno","trt","dtime","prior","age","celltype","time","status")
 dist="Exponential"
 frailty="inverse gaussian" 
 cluster="celltype"
 nMonte=1000
 formula=as.formula("Surv(time, status) ~ karno")
 R.model=bcmle(formula,D,cluster,dist,frailty,nMonte)
 plotSurv(R.model)
 RRR=opt.model(formula,D,cluster)

