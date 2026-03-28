#' @title plotSurv
#'
#' @description   Returns the plot for empirical, estimated and corrected survivals.
#'
#' @param R.model Output of function 'bcmle'.
#'
#' @return A plot.
#'
#' @examples
#' \dontrun{
#' library(survival)
#' D=veteran
#' D$prior <- factor(as.character(D$prior), labels = c(0, 1))
#' D$trt <- factor(as.character(D$trt), labels = c(0, 1))
#' cluster="cluster"
#' D=data.frame(D$karno,D$trt,D$diagtime,D$prior,D$age,D$celltype,D$time,D$status)
#' colnames(D)=c("karno","trt","dtime","prior","age","celltype","time","status")
#' dist="Exponential"
#' frailty="gamma"
#' cluster="celltype"
#' nMonte=1000
#' formula=as.formula("Surv(time, status) ~ karno")
#' R.model=bcmle(formula,D,cluster,dist,frailty,nMonte)
#' plotSurv(R.model)
#'
#' }
#'
#' @import graphics
#' @export
#'
plotSurv=function(R.model) {
 ing=R.model$sEmp[(R.model$sEmp)>0]
 lg=0.9*min(c(ing,R.model$sEst,R.model$sCor))
 pl=0.6*max(R.model$time)
 suppressWarnings(plot(R.model$time,R.model$sEmp,type="l",xlab="time (days)",log="y",ylab="survival",
 lty=1,ylim=c(lg,1.),lwd=2))
 lines(R.model$time,R.model$sEst,col="black",lty=2,lwd=2)
 lines(R.model$time,R.model$sCor,col="black",lty=3,lwd=2)
 legend(pl,1,c("empirical","estimated","corrected"),lty=1:3,bty="n",lwd=2)
 title('Empirical, estimated and corrected survivals.')
 mtext(side = 3, line=0.5,paste("Model with",R.model$dist,"baseline hazard function and",R.model$frailty,"frailty",sep=" "))
}
