#' Heterogeneity-Adjusted Multivariable MEndelian Randomization
#'
#' @description A method for multivariable cis-Mendelian randomization analysis with summary data.
#'
#' @param by beta coefficients for the outcome (a J-vector for J genetic variants)
#' @param bx beta coefficients for the exposures (a J x K matrix, where K is the number of exposures, and the k-th column are the J-vector of beta coefficients for the k-th exposure)
#' @param sy standard errors for the outcome
#' @param sx standard errors for the exposures (a J x K matrix, where the k-th column are the J-vector of standard errors for the k-th exposure)
#' @param N K+1 vector of sample sizes; the first element is the sample size for outcome associations. The (k+1)-th element is the sample size for the k-th exposure associations
#' @param ld genetic variant (LD) correlation matrix
#' @param cor.x (optional) K x K exposure correlation matrix (where K is the number of exposures); default is a diagonal matrix so the exposures are assumed to be unconditionally uncorrelated (results are not too sensitive to this choice)
#' @param cor.xy (optional) K x 1 exposure-outcome correlation vector (where K is the number of exposures); default is a vector of zeros, so that the outcome and exposures are unconditionally uncorrelated (again, results are not too sensitive to this choice)
#' @param alpha (optional) (1-alpha)*100% is the nominal coverage of confidence intervals reported in the table of resuls; default is 0.05
#' @param robust (optional) logical: if \code{TRUE} then the results are adjusted for overdispersion heterogeneity in genetic variant associations with the outcome if such heterogeneity is detected; default is \code{TRUE}
#' @param pc.thres (optional) the proportion of variance in genetic associations that the principal components of genetic variants used as instruments are required to explain; default is 0.99 so that at least 99% of variation in genetic variants will be explained by the instruments
#' @param clump.thres (optional) R^2 clumping threshold for variant correlations; if unspecified the default value is R^2 = 0.95
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
#'  \item \code{numberIVs} \cr
#'  the number of instruments used by the method
#' }
#'
#'
#' @author Jack Bowden, Stephen Burgess, Ashish Patel, Dmitry Shungin
#'
#' @export

hammer <-  function(by,bx,sy,sx,N,alpha=0.05,ld=NULL,cor.x=NULL,cor.xy=NULL,robust=TRUE,pc.thres=NULL,clump.thres=NULL,exposure.names=NULL){
bx <- as.matrix(bx); sx <- as.matrix(sx); J = nrow(bx); K = ncol(bx); N <- as.vector(N)
if(J<=K){stop("We require more instruments than the number of exposures.")}
cor.x0 <- 1; cor.xy0 <- 1; ld0 <- 1
if(missing(cor.x)){cor.x0 <- 0; cor.x <- diag(K)}
if(missing(cor.xy)){cor.xy0 <- 0; cor.xy <- rep(0,K)}
if(missing(ld)){ld0 <- 0}
if(missing(pc.thres)){pc.thres=0.99}
if(missing(clump.thres)){clump.thres=0.95}
if(missing(exposure.names)){exposure.names <- as.character(1:K)}
if(length(N)!=(K+1)){stop("N is required to be a K+1 vector of sample sizes; the first element is the sample size for outcome associations. The (k+1)-th element is the sample size for the k-th exposure associations.")}

if(ld0==1){

  ### clump variants
  clumped <- clump(by=by,sy=sy,bx=bx,sx=sx,ld=ld,clump.thres=clump.thres)
  by=clumped$by; sy=clumped$sy; bx=clumped$bx; sx=clumped$sx; ld=clumped$ld; clump.thres=clumped$clump.thres; rm(clumped)

  ### transform univariable summary statistics into PCA-transformed univariable statistics
  pca <- PCsumstats(by=by,sy=sy,bx=bx,sx=sx,ld=ld,N,pc.thres=pc.thres)
  by=pca$by; sy=pca$sy; bx=pca$bx; sx=pca$sx; pc.thres=pca$pc.thres; rm(pca)
  bx <- as.matrix(bx); sx <- as.matrix(sx); J = nrow(bx); K = ncol(bx)
  if(J<=K){stop("We require more instruments than the number of exposures. Consider using a higher threshold for pc.thres in order to use more instruments.")}
}

### transform univariable summary data into multivariable summary data
multi.dat <- MVsumstats(by=by,bx=bx,sy=sy,sx=sx,N=N,cor.x=cor.x,cor.xy=cor.xy)
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

### Step 1: preliminary GMM theta estimator assuming no overdispersion heterogeneity
Om <- function(tet){
  Om0 <- vector(,length=J)
  for (j in 1:J){Om0[j] <- (SigY[j]-2*as.numeric(t(tet)%*%SigXY[[j]])+as.numeric(t(tet)%*%SigX[[j]]%*%tet))}
  return(Om0)
}
g <- function(tet){GamY-as.vector(GamX%*%tet)}
Q.gg <- function(tet){as.numeric(t(g(tet))%*%g(tet))}
est.init <- nlminb(rep(0,K),Q.gg,lower=rep(-5,K),upper=rep(5,K))$par
Q <- function(tet){as.numeric(t(g(tet)/Om(tet))%*%g(tet))}
gmm <- nlminb(est.init,Q,lower=rep(-5,K),upper=rep(5,K))$par; rm(Q.gg,est.init)
V1 <- t(GamX/Om(gmm))%*%GamX

### if robust == FALSE
if(robust==FALSE){
est <- gmm
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
cat("\n------------------------------------------------------------------\n")
print(output.table, quote = F, row.names = FALSE, justify = "left")
cat("------------------------------------------------------------------\n")

cat("\nNon-robust model with no overdispersion heterogeneity was selected by the user.")
cat("\nNumber of instruments used :", nrow(bx))
if(ld0==1){cat("\nGenetic variants were correlated. The principal components of genetic variants were used as instruments. \nPrincipal components that explained",pc.thres*100,"% of genetic variation were used as instruments.")}
if(cor.x0==0){cat("\nExposure correlation matrix was not supplied. Exposures were assumed to be either mutually uncorrelated, or genetic associations with each exposure were assumed to be measured from non-overlapping samples.")}
if(cor.xy0==0){cat("\nExposure-outcome correlations were not supplied. The exposures were assumed to be either uncorrelated with the outcome, or genetic associations with each exposure were assumed to be measured from an non-overlapping sample with genetic associations with the outcome.")}

return(list("est"=est,"se"=se.mwi,"condF"=condF,"numberIVs"=J))
}

### if robust == TRUE
if(robust==TRUE){
### Step 2: overdispersion heterogeneity parameter inference
M <- as.vector((GamX/Om(gmm))%*%t(t(colMeans(GamX))%*%solve(V1)))
M0 <- (1-M)^2
kappa.est <- N[1]*(sum(g(gmm)^2)-sum(M0*Om(gmm)))/sum(M0)
Om.kap <- Om(gmm)+(kappa.est/N[1])
kappa.var <- (2*sum((M0^2)*(Om.kap^2)))/(sum(M0)^2)
kappa.se <- sqrt(kappa.var)

### Step 3: Heterogeneity-Adjusted Multivariable MEndelian Randomization
sig.epgam <- function(tet){
  sig.epgam0 <- matrix(,nrow=K,ncol=J)
  for (j in 1:J){sig.epgam0[,j] <- (SigXY[[j]]-(SigX[[j]]%*%tet))}
  return(sig.epgam0)
}
Om.star <- function(tet){Om(tet)+(max(kappa.est,0)/N[1])}
Q.star <- function(tet){as.numeric(t(g(tet)/Om.star(tet))%*%g(tet))}
gmm.star <- nlminb(gmm,Q.star,lower=rep(-5,K),upper=rep(5,K))$par
Om.star <- function(kappa){Om(gmm.star)+(max(kappa,0)/N[1])}

sig.epgam <- sig.epgam(gmm.star)
V1.star <- function(kappa){(t(GamX/Om.star(kappa))%*%GamX)}
V2.star <- function(kappa){
  V2A <- list()
  for (j in 1:J){V2A[[j]] <- SigX[[j]]/(Om.star(kappa)[j])}
  V2A <- Reduce('+',V2A)
  V2B <- sig.epgam%*%((t(sig.epgam)/(Om.star(kappa)^2)))
  return(V2A-V2B)
}
var.star <- function(kappa){solve(V1.star(kappa))%*%(V1.star(kappa)+V2.star(kappa))%*%t(solve(V1.star(kappa)))}
se.star <- sqrt(diag(var.star(kappa.est)))

decimals <- function(number, places) format(round(number, places), nsmall = places)
Interval_type <- paste(100*(1-alpha), "% CI", sep = "")
if(K>1){
  Statistic <- c("Exposure", "Estimate", "Std Error", Interval_type, "", "P-value", "Cond F-stat")
}
if(K==1){
  Statistic <- c("Exposure", "Estimate", "Std Error", Interval_type, "", "P-value", "F-stat")
}
CILower <- gmm.star-(qnorm(1-alpha/2)*se.star); CIUpper <- gmm.star+(qnorm(1-alpha/2)*se.star)
Pvalue <- 2*(1-pnorm(abs(gmm.star/se.star)))
Value <- cbind(exposure.names, decimals(gmm.star, 3), decimals(se.star,3),
               paste(decimals(CILower, 3), ",", sep = ""), decimals(CIUpper,3),
               decimals(Pvalue, 3), decimals(condF, 3))

output.table <- data.frame(matrix(Value, nrow = K))
colnames(output.table) <- Statistic

cat("\nHeterogeneity-Adjusted Multivariable MEndelian Randomization (HAMMER) method\n")
cat("\n------------------------------------------------------------------\n")
print(output.table, quote = F, row.names = FALSE, justify = "left")
cat("------------------------------------------------------------------\n")

cat("\nOverdispersion heterogeneity parameter estimate =", decimals(max(kappa.est,0),3))
cat("\nNumber of instruments used :", nrow(bx))
if(ld0==1){cat("\nGenetic variants were correlated. The principal components of genetic variants were used as instruments. \nPrincipal components that explained",pc.thres*100,"% of genetic variation were used as instruments.")}
if(cor.x0==0){cat("\nExposure correlation matrix was not supplied. Exposures were assumed to be either mutually uncorrelated, or genetic associations with each exposure were assumed to be measured from non-overlapping samples.")}
if(cor.xy0==0){cat("\nExposure-outcome correlations were not supplied. The exposures were assumed to be either uncorrelated with the outcome, or genetic associations with each exposure were assumed to be measured from an non-overlapping sample with genetic associations with the outcome.")}

return(list("est"=gmm.star,"se"=se.star,"kappa"=kappa.est,"kappa.se"=kappa.se,"condF"=condF,"numberIVs"=J))
}
}
