#' PCsumstats function
#'
#' @description A method for calculating summary association data using principal components of many highly correlated variants
#'
#' @param by beta coefficients for the outcome (a J-vector for J genetic variants)
#' @param bx beta coefficients for the exposures (a J x K matrix, where K is the number of exposures, and the k-th column are the J-vector of beta coefficients for the k-th exposure)
#' @param sy standard errors for the outcome
#' @param sx standard errors for the exposures (a J x K matrix, where the k-th column are the J-vector of standard errors for the k-th exposure)
#' @param N sample size; a one-sample analysis is performed if only one sample size is entered. If two sample sizes N=c(Ny,Nx) are entered, then a two-sample analysis is performed under the assumption that genetic associations with traits with the outcome (Ny-sample) and exposures (Nx-sample) are measured from non-overlapping samples
#' @param ld genetic variant (LD) correlation matrix
#' @param pc.thres (optional) the proportion of variance in genetic associations that the principal components of genetic variants used as instruments are required to explain; default is 0.99 so that at least 99% of variation in the genetic variant LD matrix will be explained by the instruments
#'
#' @details This method transforms univariable genetic associations based a set of many highly correlated genetic variants into univariable genetic associations based on a smaller set of uncorrelated principal components of the many correlated variants. \cr
#'
#' @return Output is a list containing:
#' \itemize{
#'  \item \code{by} \cr
#'  Vector of outcome beta-coefficients for the the principal components of genetic variants
#'  \item \code{sy} \cr
#'  Vector of outcome standard errors for the the principal components of genetic variants
#'  \item \code{bx} \cr
#'  Matrix of exposure beta-coefficients for the the principal components of genetic variants (the k-th column corresponds to the k-th exposure)
#'  \item \code{sx} \cr
#'  Matrix of exposure standard errors for the the principal components of genetic variants (the k-th column corresponds to the k-th exposure)
#'  \item \code{pc.thres} \cr
#'  the proportion of variance in genetic associations that the principal components of genetic variants used as instruments explain
#' }
#'
#'
#' @author Jack Bowden, Stephen Burgess, Ashish Patel, Dmitry Shungin
#'
#' @export

PCsumstats <- function(by,sy,bx,sx,ld,N,pc.thres){
  bx <- as.matrix(bx); sx <- as.matrix(sx); J = nrow(bx); K = ncol(bx); N <- as.vector(N)
  if(length(N)==1){
    # calculate principal components
    ay <- 1/((N*sy^2)+by^2)
    Ay <- (sqrt(ay)%*%t(sqrt(ay)))*ld
    Phi <- Ay
    r <- which(cumsum(prcomp(Phi,scale=FALSE)$sdev^2/sum((prcomp(Phi,scale=FALSE)$sdev^2)))>pc.thres)[1]
    lambda <- sqrt(J)*prcomp(Phi,scale=FALSE)$rotation[,1:r]
    evec <- eigen((t(lambda)%*%lambda))$vectors
    eval <- eigen((t(lambda)%*%lambda))$values
    lambda <- lambda%*%(solve(evec%*%diag(sqrt(eval))%*%t(evec)))
    dim(lambda) <- c(J,r)

    # transform univariable summary data
    transform.sumstats <- function(beta,se,ld,N,lambda){
      a <- 1/(((N-1)*se^2)+beta^2)
      A <- (sqrt(a)%*%t(sqrt(a)))*ld
      vA <- diag(t(lambda)%*%A%*%lambda)
      B <- a*beta; vB <- as.vector(t(lambda)%*%B)
      beta.til <- vB/vA
      se.til <- sqrt(((1/vA)-(beta.til^2))/(N-1))
      return(list("beta"=beta.til, "se"=se.til))
    }
    by.t <- transform.sumstats(beta=by,se=sy,ld=ld,N=N,lambda=lambda)
    se.betaY <- by.t$se; betaY <- by.t$beta; rm(by.t)
    betaX <- matrix(NA,nrow=r,ncol=K); se.betaX <- matrix(NA,nrow=r,ncol=K)
    for (k in 1:K){
      bx.t <- transform.sumstats(beta=bx[,k],se=sx[,k],ld=ld,N=N,lambda=lambda)
      se.betaX[,k] <- bx.t$se; betaX[,k] <- bx.t$beta; rm(bx.t)
    }

    # return transformed summary data
    return(list("by"=betaY,"sy"=se.betaY,"bx"=betaX,"sx"=se.betaX,"pc.thres"=pc.thres))
  }
  if(length(N)==2){

    # calculate principal components
    ay <- 1/((N[1]*sy^2)+by^2)
    Ay <- (sqrt(ay)%*%t(sqrt(ay)))*ld
    Phi <- Ay
    r <- which(cumsum(prcomp(Phi,scale=FALSE)$sdev^2/sum((prcomp(Phi,scale=FALSE)$sdev^2)))>pc.thres)[1]
    lambda <- sqrt(J)*prcomp(Phi,scale=FALSE)$rotation[,1:r]
    evec <- eigen((t(lambda)%*%lambda))$vectors
    eval <- eigen((t(lambda)%*%lambda))$values
    lambda <- lambda%*%(solve(evec%*%diag(sqrt(eval))%*%t(evec)))
    dim(lambda) <- c(J,r)

    # transform univariable summary data
    transform.sumstats <- function(beta,se,ld,N,lambda){
      a <- 1/(((N-1)*se^2)+beta^2)
      A <- (sqrt(a)%*%t(sqrt(a)))*ld
      vA <- diag(t(lambda)%*%A%*%lambda)
      B <- a*beta; vB <- as.vector(t(lambda)%*%B)
      beta.til <- vB/vA
      se.til <- sqrt(((1/vA)-(beta.til^2))/(N-1))
      return(list("beta"=beta.til, "se"=se.til))
    }
    by.t <- transform.sumstats(beta=by,se=sy,ld=ld,N=N[1],lambda=lambda)
    se.betaY <- by.t$se; betaY <- by.t$beta; rm(by.t)
    betaX <- matrix(NA,nrow=r,ncol=K); se.betaX <- matrix(NA,nrow=r,ncol=K)
    for (k in 1:K){
      bx.t <- transform.sumstats(beta=bx[,k],se=sx[,k],ld=ld,N=N[2],lambda=lambda)
      se.betaX[,k] <- bx.t$se; betaX[,k] <- bx.t$beta; rm(bx.t)
    }

    # return transformed summary data
    return(list("by"=betaY,"sy"=se.betaY,"bx"=betaX,"sx"=se.betaX,"pc.thres"=pc.thres))
  }
}
