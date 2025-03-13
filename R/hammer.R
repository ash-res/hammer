#' Heterogeneity-Adjusted Multivariable MEndelian Randomization
#'
#' @description A method for multivariable cis-Mendelian randomization analysis with summary data.
#'
#' @param by beta coefficients for the outcome (a J-vector for J genetic variants)
#' @param bx beta coefficients for the exposures (a J x K matrix, where K is the number of exposures, and the k-th column are the J-vector of beta coefficients for the k-th exposure)
#' @param sy standard errors for the outcome
#' @param sx standard errors for the exposures (a J x K matrix, where the k-th column are the J-vector of standard errors for the k-th exposure)
#' @param N sample size; a one-sample analysis is performed if only one sample size is entered. If two sample sizes N=c(Ny,Nx) are entered, then a two-sample analysis is performed under the assumption that genetic associations with traits with the outcome (Ny-sample) and exposures (Nx-sample) are measured from non-overlapping samples
#' @param ld genetic variant (LD) correlation matrix
#' @param cor.x (optional) K x K exposure correlation matrix (where K is the number of exposures); default is a diagonal matrix so the exposures are assumed to be unconditionally uncorrelated (results are not too sensitive to this choice)
#' @param cor.xy (optional) K x 1 exposure-outcome correlation vector (where K is the number of exposures); default is a vector of zeros, so that the outcome and exposures are unconditionally uncorrelated (again, results are not too sensitive to this choice)
#' @param alpha (optional) (1-alpha)*100% is the nominal coverage of confidence intervals reported in the table of resuls; default is 0.05
#' @param robust (optional) logical: if \code{TRUE} then the results are adjusted for overdispersion heterogeneity in genetic variant associations with the outcome if such heterogeneity is detected; default is \code{TRUE}
#' @param pc.thres (optional) the proportion of variance in genetic associations that the principal components of genetic variants used as instruments are required to explain; default is 0.99 so that at least 99% of variation in genetic variants will be explained by the instruments
#' @param prune.thres (optional) R^2 pruning threshold for variant correlations; if unspecified the default value is R^2 = 0.95
#' @param exposure.names (optional) a J x 1 character vector of the names of exposures included in the model
#'
#' @details This method uses the principal components of genetic variants in a single gene region to instrument multiple exposures for Mendelian randomization randomization analyses.
#'
#' @return Output is a list containing:
#' \itemize{
#'  \item \code{est} \cr
#'  K x 1 vector of exposure causal effect estimates on the outcome
#'  \item \code{se} \cr
#'  K x 1 vector of standard errors relating to the exposure causal effect estimates
#'  \item \code{condF} \cr
#'  K x 1 vector of conditional F-statistics (only reported if number of exposures is K > 1), otherwise a scalar F-statistic in the univariable exposure-outcome model is reported
#'  \item \code{kappa} \cr
#'  the estimated overdispersion heterogeneity parameter (a parameter estimate > 0 means the model was adjusted for overdispersion heterogeneity); only reported if the input \code{robust=TRUE}
#' }
#'
#'
#' @author Jack Bowden, Stephen Burgess, Ashish Patel, Dmitry Shungin
#'
#' @export

hammer <- function(by,bx,sy,sx,N,alpha=0.05,ld=NULL,cor.x=NULL,cor.xy=NULL,robust=TRUE,pc.thres=NULL,prune.thres=NULL,exposure.names=NULL){
bx <- as.matrix(bx); sx <- as.matrix(sx); J = nrow(bx); K = ncol(bx); N <- as.vector(N)
if(J<=K){stop("We require more instruments than the number of exposures.")}
cor.x0 <- 1; cor.xy0 <- 1; ld0 <- 1
if(missing(cor.x)){cor.x0 <- 0; cor.x <- diag(K)}
if(missing(cor.xy)){cor.xy0 <- 0; cor.xy <- rep(0,K)}
if(missing(ld)){ld0 <- 0}
if(missing(pc.thres)){pc.thres=0.99}
if(missing(prune.thres)){prune.thres=0.95}
if(missing(exposure.names)){exposure.names <- as.character(1:K)}
if(length(N)>2){stop("The sample size N must be a single (scalar) sample size for a one-sample analysis, or a vector of size 2 for a two-sample analysis, where the first element N[1] is the sample size for the outcome, and N[2] is the sample size for exposures.")}

if(ld0==1){
  prune <- function(by,sy,bx,sx,ld){
    # calculate t-statistics and re-order traits
    ts.y <- abs(by)/sy
    ts.x <- abs(bx)/sx

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
  if(length(N)==1){pca.sumstats <- function(by,sy,bx,sx,ld,N,pc.thres){
    J = nrow(bx); K = ncol(bx); bx <- as.matrix(bx); sx <- as.matrix(sx)

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
  }}
  if(length(N)==2){pca.sumstats <- function(by,sy,bx,sx,ld,N,pc.thres){
    J = nrow(bx); K = ncol(bx); bx <- as.matrix(bx); sx <- as.matrix(sx)

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

  ### prune variants
  pruned <- prune(by=by,sy=sy,bx=bx,sx=sx,ld=ld)
  by=pruned$by; sy=pruned$sy; bx=pruned$bx; sx=pruned$sx; ld=pruned$ld; prune.thres=pruned$prune.thres; rm(pruned)

  ### transform univariable summary statistics into PCA-transformed univariable statistics
  pca <- pca.sumstats(by=by,sy=sy,bx=bx,sx=sx,ld=ld,N,pc.thres=pc.thres)
  by=pca$by; sy=pca$sy; bx=pca$bx; sx=pca$sx; pc.thres=pca$pc.thres; rm(pca)
  bx <- as.matrix(bx); sx <- as.matrix(sx); J = nrow(bx); K = ncol(bx)
  if(J<=K){stop("We require more instruments than the number of exposures. Consider using a higher threshold for pc.thres in order to use more instruments.")}
}

### transform univariable summary data into multivariable summary data
if(length(N)==1){multi.dat <- function(by,bx,sy,sx,N){
# variance-covariance matrix of exposure associations
GamX <- bx; GamY <- by
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
}}
if(length(N)==2){multi.dat <- function(by,bx,sy,sx,N){
  # variance-covariance matrix of exposure associations
  GamX <- bx; GamY <- by
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
}}
multi.dat <- multi.dat(by=by,bx=bx,sy=sy,sx=sx,N=N)
GamY=multi.dat$GamY; GamX=multi.dat$GamX; SigY=multi.dat$SigY; SigX=multi.dat$SigX; SigXY=multi.dat$SigXY; rm(multi.dat)

### compute conditional F-statistics
condF <- function(k){
g.f <- function(delta){GamX[,k]-as.vector(as.matrix(GamX[,-k])%*%delta)}
Q.f <- function(delta){as.numeric(t(g.f(delta))%*%g.f(delta))}
delta.f <- nlminb(rep(0,(K-1)),objective=Q.f)$par
Om.f <- function(delta){
  Om0.f <- vector(,length=J)
  for (j in 1:J){Om0.f[j] <- as.numeric(SigX[[j]][k,k]-(2*(t(delta)%*%as.matrix(SigX[[j]][k,-k])))+(t(delta)%*%as.matrix(SigX[[j]][-k,-k])%*%delta))}
  return(Om0.f)
}
Q.f <- function(delta){as.numeric(t(g.f(delta)/Om.f(delta))%*%g.f(delta))}
DQ.f <- function(delta){-2*as.matrix(t(GamX[,-k])%*%(g.f(delta)/Om.f(delta)))}
condf <- nlminb(delta.f,objective=Q.f,gradient=DQ.f)$objective/(J-K+1)
return(condf)
}
if(K>1){condF <- sapply(1:K,condF)}
if(K==1){condF <- sum(as.vector(GamX^2)/unlist(SigX))/J}

### Step 1: preliminary theta estimator assuming no overdispersion heterogeneity
Om <- function(tet){
  Om0 <- vector(,length=J)
  for (j in 1:J){Om0[j] <- (SigY[j]-2*as.numeric(t(tet)%*%SigXY[[j]])+as.numeric(t(tet)%*%SigX[[j]]%*%tet))}
  return(Om0)
}
sig.epgam <- function(tet){
  sig.epgam0 <- matrix(,nrow=K,ncol=J)
  for (j in 1:J){sig.epgam0[,j] <- (SigXY[[j]]-(SigX[[j]]%*%tet))}
  return(sig.epgam0)
}
g <- function(tet){GamY-as.vector(GamX%*%tet)}
Q.gg <- function(tet){as.numeric(t(g(tet))%*%g(tet))}
est.init <- nlminb(rep(0,K),Q.gg,lower=rep(-5,K),upper=rep(5,K))$par
Q <- function(tet){as.numeric(t(g(tet)/Om(tet))%*%g(tet))}
est.init <- nlminb(est.init,Q,lower=rep(-5,K),upper=rep(5,K))$par

### if robust == FALSE
if(robust==FALSE){
est <- est.init
Om <- Om(est)
sig.epgam <- matrix(,nrow=K,ncol=J)
for (j in 1:J){sig.epgam[,j] <- (SigXY[[j]]-(SigX[[j]]%*%est))}
V1 <- (t(GamX/Om)%*%GamX)
V2A <- list()
for (j in 1:J){V2A[[j]] <- SigX[[j]]/(Om[j])}
V2A <- Reduce('+',V2A)
V2B <- sig.epgam%*%((t(sig.epgam)/(Om^2)))
V2 <- V2A+V2B
V0 <- V1+V2
D0 <- -V1
var.est <- solve(D0)%*%V0%*%t(solve(D0))
se.mwi <- sqrt(diag(var.est))

decimals <- function(number, places) format(round(number, places), nsmall = places)
Interval_type <- paste(100*(1-alpha), "% CI", sep = "")
if(K>1){
  Statistic <- c("Exposure", "Estimate", "Std Error", Interval_type, "", "P-value", "Cond F-stat")
}
if(K==1){
  Statistic <- c("Exposure", "Estimate", "Std Error", Interval_type, "", "P-value", "F-stat")
}
  CILower <- est-(qnorm(1-alpha/2)*se.mwi); CIUpper <- est+(qnorm(1-alpha/2)*se.mwi)
  Pvalue <- 2*(1-pnorm(abs(est/se.mwi)))
  Value <- cbind(exposure.names, decimals(est, 3), decimals(se.mwi,3),
               paste(decimals(CILower, 3), ",", sep = ""), decimals(CIUpper,3),
               decimals(Pvalue, 3), decimals(condF, 3))
output.table <- data.frame(matrix(Value, nrow = K))
colnames(output.table) <- Statistic

cat("\nHeterogeneity-Adjusted Multivariable MEndelian Randomization (HAMMER) method\n")
cat(c("\nOne-sample results","\nTwo-sample results")[length(N)])
cat("\n------------------------------------------------------------------\n")
print(output.table, quote = F, row.names = FALSE, justify = "left")
cat("------------------------------------------------------------------\n")

cat("\nNon-robust model with no overdispersion heterogeneity was selected by the user.")
cat("\nNumber of instruments used :", nrow(bx))
if(ld0==1){cat("\nGenetic variants were correlated. The principal components of genetic variants were used as instruments. \nPrincipal components that explained",pc.thres*100,"% of genetic variation were used as instruments.")}
if(cor.x0==0){cat("\nExposure correlation matrix was not supplied. Exposures were assumed to be mutually uncorrelated.")}
if(cor.xy0==0){cat("\nExposure-outcome correlations were not supplied. The exposures were assumed to be uncorrelated with the outcome.")}

return(list("est"=est,"se"=se.mwi,"condF"=condF))
}

### if robust == TRUE
if(robust==TRUE){
### Step 2: overdispersion parameter estimator
g.kap <- GamY-as.vector(GamX%*%est.init)
if(length(N)==1){Om.kap <- function(kap){
  Om0.kap <- vector(,length=J)
  for (j in 1:J){Om0.kap[j] <- (SigY[j]+(kap/N)-2*as.numeric(t(est.init)%*%SigXY[[j]])+as.numeric(t(est.init)%*%SigX[[j]]%*%est.init))}
  return(Om0.kap)
}}
if(length(N)==2){Om.kap <- function(kap){
  Om0.kap <- vector(,length=J)
  for (j in 1:J){Om0.kap[j] <- (SigY[j]+(kap/N[1])-2*as.numeric(t(est.init)%*%SigXY[[j]])+as.numeric(t(est.init)%*%SigX[[j]]%*%est.init))}
  return(Om0.kap)
}}

V1.kap <- function(kap){(t(GamX/Om.kap(kap))%*%GamX)}
adj1 <- function(kap){
adj10 <- vector(,length=J)
for (j in 1:J){adj10[j] <- -2*as.numeric(t(GamX[j,])%*%solve(V1.kap(kap))%*%GamX[j,])}
return(adj10)
}
adj2 <- function(kap){
adj20 <- vector(,length=J)
for (j in 1:J){adj20[j] <- (as.numeric(t(GamX[j,])%*%solve(V1.kap(kap))%*%GamX[j,])^2)*(1/Om.kap(kap)[j])}
return(adj20)
}
adj1 <- adj1(0); adj2 <- adj2(0)
Q.kap <- function(kap){as.numeric(t(((g.kap^2)-Om.kap(kap)-adj1-adj2))%*%((g.kap^2)-Om.kap(kap)-adj1-adj2))}
init.val <- seq(-0.5,2,0.5)
Q.kap.init <- function(l){optim(init.val[l], Q.kap, method="Brent",lower=-5,upper=50)$value}
Q.kap.init <- sapply(1:length(init.val),Q.kap.init)
kap.est <- optim(init.val[which.min(Q.kap.init)[[1]]], Q.kap, method="Brent",lower=-5,upper=50)$par

### Step 3: Robust theta estimator
if(length(N)==1){Om <- function(tet){
  Om0 <- vector(,length=J)
  for (j in 1:J){Om0[j] <- (SigY[j]+max((kap.est/N),0)-2*as.numeric(t(tet)%*%SigXY[[j]])+as.numeric(t(tet)%*%SigX[[j]]%*%tet))}
  return(Om0)
}}
if(length(N)==2){Om <- function(tet){
  Om0 <- vector(,length=J)
  for (j in 1:J){Om0[j] <- (SigY[j]+max((kap.est/N[1]),0)-2*as.numeric(t(tet)%*%SigXY[[j]])+as.numeric(t(tet)%*%SigX[[j]]%*%tet))}
  return(Om0)
}}
sig.epgam <- function(tet){
  sig.epgam0 <- matrix(,nrow=K,ncol=J)
  for (j in 1:J){sig.epgam0[,j] <- (SigXY[[j]]-(SigX[[j]]%*%tet))}
  return(sig.epgam0)
}
g <- function(tet){GamY-as.vector(GamX%*%tet)}
Q <- function(tet){as.numeric(t(g(tet)/Om(tet))%*%g(tet))}
est <- nlminb(est.init,Q,lower=rep(-5,K),upper=rep(5,K))$par
# inference
Om <- Om(est)
sig.epgam <- sig.epgam(est)
V1 <- (t(GamX/Om)%*%GamX)
V2A <- list()
for (j in 1:J){V2A[[j]] <- SigX[[j]]/(Om[j])}
V2A <- Reduce('+',V2A)
V2B <- sig.epgam%*%((t(sig.epgam)/(Om^2)))
V2 <- V2A+V2B
V0 <- V1+V2
D0 <- -V1
var.est <- solve(D0)%*%V0%*%t(solve(D0))
se.mwi <- sqrt(diag(var.est))

decimals <- function(number, places) format(round(number, places), nsmall = places)
Interval_type <- paste(100*(1-alpha), "% CI", sep = "")
if(K>1){
  Statistic <- c("Exposure", "Estimate", "Std Error", Interval_type, "", "P-value", "Cond F-stat")
}
if(K==1){
  Statistic <- c("Exposure", "Estimate", "Std Error", Interval_type, "", "P-value", "F-stat")
}
CILower <- est-(qnorm(1-alpha/2)*se.mwi); CIUpper <- est+(qnorm(1-alpha/2)*se.mwi)
Pvalue <- 2*(1-pnorm(abs(est/se.mwi)))
Value <- cbind(exposure.names, decimals(est, 3), decimals(se.mwi,3),
               paste(decimals(CILower, 3), ",", sep = ""), decimals(CIUpper,3),
               decimals(Pvalue, 3), decimals(condF, 3))

output.table <- data.frame(matrix(Value, nrow = K))
colnames(output.table) <- Statistic

cat("\nHeterogeneity-Adjusted Multivariable MEndelian Randomization (HAMMER) method\n")
cat(c("\nOne-sample results","\nTwo-sample results")[length(N)])
cat("\n------------------------------------------------------------------\n")
print(output.table, quote = F, row.names = FALSE, justify = "left")
cat("------------------------------------------------------------------\n")

cat("\nOverdispersion heterogeneity parameter estimate =", decimals(max(kap.est,0),3))
cat("\nNumber of instruments used :", nrow(bx))
if(ld0==1){cat("\nGenetic variants were correlated. The principal components of genetic variants were used as instruments. \nPrincipal components that explained",pc.thres*100,"% of genetic variation were used as instruments.")}
if(cor.x0==0){cat("\nExposure correlation matrix was not supplied. Exposures were assumed to be mutually uncorrelated.")}
if(cor.xy0==0){cat("\nExposure-outcome correlations were not supplied. The exposures were assumed to be uncorrelated with the outcome.")}

return(list("est"=est,"se"=se.mwi,"kappa"=kap.est,"condF"=condF))
}
}
