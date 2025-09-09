#' MVsumstats function
#'
#' @description A method for transforming one-at-a-time uncorrelated univariable regression summary statistics into multivariate regression summary statistics
#'
#' @param by beta coefficients for the outcome (a J-vector for J genetic variants)
#' @param bx beta coefficients for the exposures (a J x K matrix, where K is the number of exposures, and the k-th column are the J-vector of beta coefficients for the k-th exposure)
#' @param sy standard errors for the outcome
#' @param sx standard errors for the exposures (a J x K matrix, where the k-th column are the J-vector of standard errors for the k-th exposure)
#' @param N sample size; a one-sample analysis is performed if only one sample size is entered. If two sample sizes N=c(Ny,Nx) are entered, then a two-sample analysis is performed under the assumption that genetic associations with traits with the outcome (Ny-sample) and exposures (Nx-sample) are measured from non-overlapping samples
#' @param ld genetic variant (LD) correlation matrix
#' @param cor.x (optional) K x K exposure correlation matrix (where K is the number of exposures); default is a diagonal matrix so the exposures are assumed to be unconditionally uncorrelated (results are not too sensitive to this choice)
#' @param cor.xy (optional) K x 1 exposure-outcome correlation vector (where K is the number of exposures); default is a vector of zeros, so that the outcome and exposures are unconditionally uncorrelated (again, results are not too sensitive to this choice)
#'
#' @details This method transforms one-at-a-time uncorrelated univariable regression summary statistics into multivariate regression summary statistics for the joint vector of instrument associations with the exposure and outcome.
#'
#' @return Output is a list containing:
#' \itemize{
#'  \item \code{GamY} \cr
#'  Vector of regression coefficients from regressing the outcome on the joint J-vector of instruments
#'  \item \code{GamX} \cr
#'  Matrix of regression coefficients where the k-th column corresponds to regressing the k-th exposure on the joint J-vector of instruments
#'  \item \code{SigY} \cr
#'  a J-vector of variances, where the j-th element corresponds to the variance of the j-th element of GamY
#'  \item \code{SigX} \cr
#'  a J-list of K x K matrices covariance matrices where the j-th list corresponds to the variance-covariance matrix for the j-th row of GamX
#'  \item \code{SigXY} \cr
#'  a J-list of K-vectors, where the j-th element corresponds to the covariance of the j-th row of GamX with j-th element of GamY
#' }
#'
#'
#' @author Jack Bowden, Stephen Burgess, Ashish Patel, Dmitry Shungin
#'
#' @export

MVsumstats <- function(by,bx,sy,sx,N,cor.x=NULL,cor.xy=NULL){
if(missing(cor.x)){cor.x <- diag(ncol(bx))}
if(missing(cor.xy)){cor.xy <- rep(0,ncol(bx))}
if(length(N)==1){
# variance-covariance matrix of exposure associations
GamX <- bx; GamY <- by; J <- nrow(GamX); K <- ncol(GamX)
ax <- matrix(NA,nrow=J,ncol=K)
for (k in 1:K){ax[,k] <- 1/((N*sx[,k]^2)+bx[,k]^2)}
Bx <- function(k){ax[,k]*bx[,k]}; Bx <- sapply(1:K,Bx)
SigX <- function(j){
  SigX.mat <- matrix(NA,nrow=K,ncol=K)
  for (k1 in 1:K){
    for (k2 in 1:K){
      SigX.mat[k1,k2] <- ((1/sqrt(ax[j,k1]*ax[j,k2])))*(cor.x[k1,k2]-(sum((Bx[,k1]*Bx[,k2])*(1/sqrt(ax[,k1]*ax[,k2])))))
    }
  }
  return(SigX.mat/N)
}
SigX <- lapply(1:J,SigX)

# variance of genetic-outcome associations
ay <- 1/((N*sy^2)+by^2)
By <- ay*by
SigY <- function(j){(1-sum((By^2)/ay))/(N*ay[j])}
SigY <- sapply(1:J,SigY)

# covariance vector of genetic-exposure and genetic-outcome associations
SigXY <- function(j){
  SigXY.vector <- vector(,length=K)
  for (k in 1:K){
    SigXY.vector[k] <- ((1/sqrt(ax[j,k]*ay[j])))*(cor.xy[k]-(sum((Bx[,k]*By)*(1/sqrt(ax[,k]*ay)))))
  }
  return(SigXY.vector/N)
}
SigXY <- lapply(1:J,SigXY)
return(list("GamY"=GamY,"GamX"=GamX,"SigY"=SigY,"SigX"=SigX,"SigXY"=SigXY))
}
if(length(N)==2){
  # variance-covariance matrix of exposure associations
  GamX <- bx; GamY <- by; J <- nrow(GamX); K <- ncol(GamX)
  ax <- matrix(NA,nrow=J,ncol=K)
  for (k in 1:K){ax[,k] <- 1/((N[2]*sx[,k]^2)+bx[,k]^2)}
  Bx <- function(k){ax[,k]*bx[,k]}; Bx <- sapply(1:K,Bx)
  SigX <- function(j){
    SigX.mat <- matrix(NA,nrow=K,ncol=K)
    for (k1 in 1:K){
      for (k2 in 1:K){
        SigX.mat[k1,k2] <- ((1/sqrt(ax[j,k1]*ax[j,k2])))*(cor.x[k1,k2]-(sum((Bx[,k1]*Bx[,k2])*(1/sqrt(ax[,k1]*ax[,k2])))))
      }
    }
    return(SigX.mat/N[2])
  }
  SigX <- lapply(1:J,SigX)

  # variance of genetic-outcome associations
  ay <- 1/((N[1]*sy^2)+by^2)
  By <- ay*by
  SigY <- function(j){(1-sum((By^2)/ay))/(N[1]*ay[j])}
  SigY <- sapply(1:J,SigY)

  # covariance vector of genetic-exposure and genetic-outcome associations
  SigXY <- function(j){
    SigXY.vector <- rep(0,K)
    return(SigXY.vector)
  }
  SigXY <- lapply(1:J,SigXY)
  return(list("GamY"=GamY,"GamX"=GamX,"SigY"=SigY,"SigX"=SigX,"SigXY"=SigXY))
}
}
