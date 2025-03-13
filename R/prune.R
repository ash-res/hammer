#' Pruning function
#'
#' @description Selects relevant variants for all exposures in a way that maximises collective instrument strength for all exposures in a multivariable model
#'
#' @param by beta coefficients for the outcome (a J-vector for J genetic variants)
#' @param bx beta coefficients for the exposures (a J x K matrix, where K is the number of exposures, and the k-th column are the J-vector of beta coefficients for the k-th exposure)
#' @param sy standard errors for the outcome
#' @param sx standard errors for the exposures (a J x K matrix, where the k-th column are the J-vector of standard errors for the k-th exposure)
#' @param ld genetic variant (LD) correlation matrix
#' @param prune.thres (optional) R^2 pruning threshold for variant correlations; if unspecified the default value is R^2 = 0.95
#'
#' @details The method leaves a set of genetic variants that are collectively relevant for all exposures in the model. \cr
#'
#' @return Output is a list containing:
#' \itemize{
#'  \item \code{by} \cr
#'  Vector of outcome beta-coefficients for the genetic variants retained after pruning
#'  \item \code{sy} \cr
#'  Vector of outcome standard errors for the genetic variants retained after pruning
#'  \item \code{bx} \cr
#'  Matrix of exposure beta-coefficients for the genetic variants retained after pruning (the k-th column corresponds to the k-th exposure)
#'  \item \code{sx} \cr
#'  Matrix of exposure standard errors for the genetic variants retained after pruning (the k-th column corresponds to the k-th exposure)
#'  \item \code{ld} \cr
#'  matrix of genetic variant (LD) correlations for the genetic variants retained after pruning
#'  \item \code{prune.thres} \cr
#'  the R^2 pruning threshold used
#' }
#'
#'
#' @author Jack Bowden, Stephen Burgess, Ashish Patel, Dmitry Shungin
#'
#' @export

prune <- function(by,sy,bx,sx,ld,prune.thres){
  bx <- as.matrix(bx); sx <- as.matrix(sx); K=ncol(bx)
  # calculate t-statistics and re-order traits
  ts.y <- abs(by)/sy
  ts.x <- as.matrix(abs(bx)/sx)

  # prune variants
  sel.y <- order(ts.y,decreasing=TRUE)
  sel.x <- function(k){order(ts.x[,k],decreasing=TRUE)}; sel.x <- sapply(1:K,sel.x)
  sel <- as.vector((t(cbind(sel.x,sel.y))))
  rm(sel.y,sel.x)
  colnames(ld) <- as.character(1:nrow(ld)); rownames(ld) <- as.character(1:nrow(ld)); snps <- as.character(1:nrow(ld))[sel]
  i=1
  while(i<=length(snps)){
    if(snps[i] %in% colnames(ld)){
      del <- which(as.vector((ld[which(colnames(ld)==snps[i]),-which(colnames(ld)==snps[i])]^2))>prune.thres)
      del <- colnames(ld[-which(colnames(ld)==snps[i]),-which(colnames(ld)==snps[i])])[del] # snps to delete
      if(length(del)>0){
        snps <- snps[-which(snps %in% del)]
        bx <- as.matrix(bx[-which(colnames(ld) %in% del),]); by <- by[-which(colnames(ld) %in% del)]
        sx <- as.matrix(sx[-which(colnames(ld) %in% del),]); sy <- sy[-which(colnames(ld) %in% del)]
        ts.x <- as.matrix(ts.x[-which(colnames(ld) %in% del),]); ts.y <- ts.y[-which(colnames(ld) %in% del)]
        ld <- ld[-which(colnames(ld) %in% del),-which(colnames(ld) %in% del)]
      }
    }
    i=i+1
  }
  return(list("by"=by,"sy"=sy,"bx"=bx,"sx"=sx,"ld"=ld,"prune.thres"=prune.thres))
}
