remove(list = ls())

## data details:
## Y: the dyspnoea statuses in 1992
##      Y = 1: a worker had dyspnoea
##covariates: fexpd(X1), the gender(sex,X2), the age(X3), the dyspnoea statuses in 1981(X4) 
##      X1:the exposure to cotton dust in the factory,1 being exposure and 0 otherwise
##      X2:sex,1 for male and 0 for female
##      X3:age
##      X4:dyspnoea statuses in 1981(IV)

data_cotton = read.csv(file = "P:/MLC/GOF_0131/real_data/new_cotton/cotton.csv",
                       header = TRUE,sep = ",")
# data_cotton = read.csv(file = "D:/20241031GOF-testing/boot-test-v2/2026年重新做的realdata/cotton/cotton.csv",
#                        header = TRUE,sep = ",")


yy0 = data_cotton$y1 ##the dyspnoea statuses in 1981
yy1 = data_cotton$y2 ##the dyspnoea statuses in 1986
yy2 = data_cotton$y3 ##the dyspnoea statuses in 1981
yy = yy2
yy[is.na(yy2)] = -1
xx1= as.vector(scale(data_cotton$age))
xx2= data_cotton$expd
xx3= data_cotton$sex
xx4= data_cotton$y1 ##instrument varoable 
rr = (yy != -1)+0
cotton = cbind(yy,xx1,xx2,xx3,xx4,rr)
###1. goodness of fit to find the suitable regression model
library(rms)
ryy = yy[rr==1]
rxx1= xx1[rr==1]
rxx2= xx2[rr==1]
rxx3= xx3[rr==1]
rxx4= xx4[rr==1]
o1 = lrm(ryy ~rxx1+rxx2+rxx3+rxx4,x = T, y = T )
out= resid(o1,"gof")
out
# o1
# Logistic Regression Model
# 
# lrm(formula = ryy ~ rxx1 + rxx2 + rxx3 + rxx4, x = T, y = T)
# 
# Model Likelihood      Discrimination    Rank Discrim.    
# Ratio Test             Indexes          Indexes    
# Obs           781    LR chi2      33.85      R2       0.087    C       0.673    
# 0            700    d.f.             4      R2(4,781)0.038    Dxy     0.333    
# 1             81    Pr(> chi2) <0.0001    R2(4,217.8)0.128    gamma   0.348    
# max |deriv| 1e-06                            Brier    0.088    tau-a   0.064    
# 
# Coef    S.E.   Wald Z Pr(>|Z|)
# Intercept -2.5222 0.2212 -11.40 <0.0001
# rxx1       0.3610 0.1305   2.77 0.0057
# rxx2       0.2281 0.2533   0.90 0.3677
# rxx3      -0.0399 0.2531  -0.16 0.8747
# rxx4       1.3062 0.3080   4.24 <0.0001
# out
# Sum of squared errors     Expected value|H0 
# 68.5946073            68.5480746 
# SD                     Z 
# 0.2764991             0.1682922 
# P 
# 0.8663534 


startpar <- function(yy, xx1, xx2, xx3, xx4, rr, inip = rep(1, 5)){
  # ----------------------------------------------------------
  # Step 1: estimate xi from the regression model for Y | X, R=1
  # This part is unchanged
  # ----------------------------------------------------------
  glom <- glm(
    yy[rr == 1] ~ xx1[rr == 1] + xx2[rr == 1] + xx3[rr == 1] + xx4[rr == 1],
    family = binomial(link = "logit")
  )
  
  glom.coef <- coef(glom)
  hatxi <- as.numeric(glom.coef)
  
  # linear predictor based on hatxi
  xixx <- hatxi[1] + hatxi[2] * xx1 + hatxi[3] * xx2 + hatxi[4] * xx3 + hatxi[5] * xx4
  px <- 1 / (1 + exp(-xixx))
  qx <- 1 / (1 + exp(xixx))
  
  # ----------------------------------------------------------
  # Step 2: estimate theta
  # Modified part:
  # remove xx1, and include xx4
  # tt = theta0 + theta1*xx2 + theta2*xx3 + theta3*xx4
  #      + log(qx + exp(theta4)*px)
  # ----------------------------------------------------------
  ell <- function(par){
    tt <- par[1] + par[2] * xx2 + par[3] * xx3 + par[4] * xx4 +
      log(qx + exp(par[5]) * px)
    
    res <- sum(rr * tt + log(1 + exp(-tt)))
    return(res)
  }
  
  Grad <- function(par){
    tt <- par[1] + par[2] * xx2 + par[3] * xx3 + par[4] * xx4 +
      log(qx + exp(par[5]) * px)
    
    tem0 <- 1 - rr - 1 / (1 + exp(-tt))
    
    # derivative w.r.t. intercept
    tem1 <- sum(tem0)
    
    # derivative w.r.t. coefficients of xx2, xx3, xx4
    tem2 <- colSums((matrix(tem0) %*% matrix(1, 1, 3)) * cbind(xx2, xx3, xx4))
    tem2 <- as.numeric(tem2)
    
    # derivative w.r.t. par[5]
    tem3 <- sum(tem0 * px * exp(par[5]) / (qx + exp(par[5]) * px))
    
    c(-tem1, -tem2, -tem3)
  }
  
  hattheta <- optim(
    par = inip,
    fn = ell,
    gr = Grad,
    method = "L-BFGS-B",
    lower = rep(-5, 5),
    upper = rep(5, 5)
  )$par
  
  return(c(as.numeric(hattheta), as.numeric(hatxi)))
}


estpar <- function(yy, xx1, xx2, xx3, xx4, rr, inip = rep(1, 5))
{
  start <- startpar(yy, xx1, xx2, xx3, xx4, rr, inip)
  
  ell <- function(par){
    alp  <- par[1]
    beta <- par[2:4]
    gam  <- par[5]
    
    # xi part: unchanged
    ss <- par[6] + par[7] * xx1 + par[8] * xx2 + par[9] * xx3 + par[10] * xx4
    ppy <- 1 / (1 + exp(-ss))
    
    # theta part: remove xx1, include xx4
    tt <- alp + beta[1] * xx2 + beta[2] * xx3 + beta[3] * xx4 +
      log(1 + (exp(gam) - 1) * ppy)
    
    tem1 <- rr * ((1 - yy) * ss + log(1 + exp(-ss)))
    tem2 <- rr * tt + log(1 + exp(-tt))
    res <- sum(tem1 + tem2)
    return(res)
  }
  
  Grad <- function(par){
    alp  <- par[1]
    beta <- par[2:4]
    gam  <- par[5]
    
    # xi part: unchanged
    ss <- par[6] + par[7] * xx1 + par[8] * xx2 + par[9] * xx3 + par[10] * xx4
    px <- 1 / (1 + exp(-ss))
    qx <- 1 / (1 + exp(ss))
    
    # theta part: remove xx1, include xx4
    tt <- par[1] + par[2] * xx2 + par[3] * xx3 + par[4] * xx4 +
      log(qx + exp(par[5]) * px)
    
    # derivatives for theta part
    tem0 <- 1 - rr - 1 / (1 + exp(-tt))
    tem1 <- sum(tem0)
    
    tem2 <- colSums((matrix(tem0) %*% matrix(1, 1, 3)) * cbind(xx2, xx3, xx4))
    tem2 <- as.numeric(tem2)
    
    tem3 <- sum(tem0 * px * exp(par[5]) / (1 + (exp(par[5]) - 1) * px))
    
    # derivatives for xi part
    tt4 <- rr * (yy - px) +
      tem0 * (exp(gam) - 1) * px * qx / (1 + (exp(gam) - 1) * px)
    
    tem4 <- colSums((matrix(tt4) %*% matrix(1, 1, 5)) * cbind(1, xx1, xx2, xx3, xx4))
    tem4 <- as.numeric(tem4)
    
    c(-tem1, -tem2, -tem3, -tem4)
  }
  
  eres <- optim(
    par = start,
    fn = ell,
    gr = Grad,
    method = "L-BFGS-B",
    lower = rep(-5, 10),
    upper = rep(5, 10),
    hessian = TRUE
  )
  
  epar <- eres$par
  eval <- eres$value
  
  # xi part: unchanged
  ess <- epar[6] + epar[7] * xx1 + epar[8] * xx2 + epar[9] * xx3 + epar[10] * xx4
  ppy <- 1 / (1 + exp(-ess))
  
  # theta part: remove xx1, include xx4
  nlp <- epar[1] + epar[2] * xx2 + epar[3] * xx3 + epar[4] * xx4 +
    log(1 - ppy + exp(epar[5]) * ppy)
  
  pix <- 1 / (1 + exp(nlp))
  
  HH <- (rr - pix)^2 - pix * (1 - pix)
  nn <- length(rr)
  tstat <- mean(HH) * sqrt(nn)
  
  return(list(
    hatxi = as.numeric(epar[6:10]),
    hattheta = as.numeric(epar[1:5]),
    tstat = tstat,
    likelihood_value = eval,
    hes = eres$hessian
  ))
}



boot <- function(data, nboot = 100, inip = rep(1, 5)){
  yy  <- data[, 1]
  xx1 <- data[, 2]
  xx2 <- data[, 3]
  xx3 <- data[, 4]
  xx4 <- data[, 5]
  rr  <- data[, 6]
  
  nn <- length(rr)
  
  out <- estpar(yy, xx1, xx2, xx3, xx4, rr, inip)
  hatxi <- out$hatxi
  hattheta <- out$hattheta
  tstat <- out$tstat
  likval <- out$likelihood_value
  
  bootstat <- c()
  
  for(i in 1:nboot){ 
    set.seed(137 * i + 537)
    
    indx <- sample(1:nn, nn, replace = TRUE)
    bxx1 <- xx1[indx]
    bxx2 <- xx2[indx]
    bxx3 <- xx3[indx]
    bxx4 <- xx4[indx]
    
    # xi part: unchanged
    bhatolp <- hatxi[1] + hatxi[2] * bxx1 + hatxi[3] * bxx2 +
      hatxi[4] * bxx3 + hatxi[5] * bxx4
    bhatppy <- 1 / (1 + exp(-bhatolp))
    
    # theta part: remove bxx1, include bxx4
    bnlp <- hattheta[1] + hattheta[2] * bxx2 + hattheta[3] * bxx3 +
      hattheta[4] * bxx4 + log(1 + (exp(hattheta[5]) - 1) * bhatppy)
    
    bp <- 1 / (1 + exp(bnlp))
    brr <- rbinom(nn, 1, bp)
    byy <- rbinom(nn, 1, bhatppy)
    
    bout <- estpar(byy, bxx1, bxx2, bxx3, bxx4, brr, inip)
    bootstat <- c(bootstat, bout$tstat)
  }
  
  list(
    tstat = tstat,
    bootstat = bootstat,
    likelihood_value = likval
  )
}



### directly estimate asymptotic variance
asyVar_fun <- function(data, inip = rep(1, 5)){
  ## data = cotton
  yy  <- data[, 1]
  xx1 <- data[, 2]
  xx2 <- data[, 3]
  xx3 <- data[, 4]
  xx4 <- data[, 5]
  rr  <- data[, 6]
  nn  <- length(rr)
  
  out <- estpar(yy, xx1, xx2, xx3, xx4, rr, inip)
  hatxi <- out$hatxi
  hattheta <- out$hattheta
  tstat <- out$tstat
  likval <- out$likelihood_value
  
  ## xi part: unchanged
  hatolp <- hatxi[1] + hatxi[2] * xx1 + hatxi[3] * xx2 + hatxi[4] * xx3 + hatxi[5] * xx4
  hatppy <- 1 / (1 + exp(-hatolp))
  
  ## theta part: remove xx1, include xx4
  tt <- hattheta[1] + hattheta[2] * xx2 + hattheta[3] * xx3 +
    hattheta[4] * xx4 + log(1 + (exp(hattheta[5]) - 1) * hatppy)
  pix <- 1 / (1 + exp(tt))
  
  ## theta covariates: now xx2, xx3, xx4
  rx <- cbind(xx2, xx3, xx4)
  dim_rx <- ncol(rx)
  
  ## xi covariates remain unchanged
  cX <- cbind(1, xx1, xx2, xx3, xx4)
  colnames(cX) <- NULL
  
  alp  <- hattheta[1]
  beta <- hattheta[2:(1 + dim_rx)]   # now beta corresponds to xx2, xx3, xx4
  gam  <- hattheta[dim_rx + 2]
  
  Rat10 <- exp(gam + hatolp) / (1 + exp(hatolp))  ## A11/A01
  
  Delta <- lapply(1:nn, function(i){
    term1 <- cX[i, ] * (Rat10[i] - hatppy[i])
    term  <- c(1, rx[i, ], Rat10[i], term1)
    return(matrix(term))
  })
  Delta <- do.call("cbind", Delta)   ## dim_theta * nn
  
  Lambda <- lapply(1:nn, function(i){
    term1 <- cX[i, ] * (yy[i] - hatppy[i])
    term  <- c(rep(0, (dim_rx + 2)), term1)
    return(matrix(term))
  })
  Lambda <- do.call("cbind", Lambda) ## dim_theta * nn
  
  dim_theta <- 5 + 5
  Rmat  <- matrix(rep(1, dim_theta)) %*% matrix(rr,  nrow = 1)
  Pimat <- matrix(rep(1, dim_theta)) %*% matrix(pix, nrow = 1)
  coEmat <- matrix(rep(1, dim_theta)) %*% matrix((1 + 2 * rr - 4 * pix) * pix * (1 - pix), nrow = 1)
  
  Smat <- Rmat * Lambda + (Pimat - Rmat) * Delta
  HH <- (rr - pix)^2 - pix * (1 - pix)
  EH <- matrix(rowMeans(coEmat * Delta))   ## dim_theta * 1
  JJ <- Smat %*% t(Smat) / nn
  inv_JJ <- solve(JJ)
  KK <- HH + t(t(EH) %*% inv_JJ %*% Smat)
  asyV <- var(KK)
  
  ## asymptotic variance for theta-hat
  SigmaS <- cov(t(Smat))   ## dim_theta x dim_theta
  asyVar_theta <- inv_JJ %*% SigmaS %*% inv_JJ
  
  ## CI for theta only: theta = (alp, beta_xx2, beta_xx3, beta_xx4, gam)
  ci_theta <- matrix(0, dim_rx + 2, 6)
  hes_inv_diag <- diag(solve(out$hes))
  
  for(i in 1:(dim_rx + 2)){
    SE <- sqrt(hes_inv_diag[i])
    waldZ <- hattheta[i] / SE
    pvalue <- 2 * (1 - pnorm(abs(waldZ)))
    ci_theta[i, ] <- c(
      hattheta[i], SE, waldZ, pvalue,
      hattheta[i] - 1.96 * SE,
      hattheta[i] + 1.96 * SE
    )
  }
  
  rownames(ci_theta) <- c("alp", "beta_xx2", "beta_xx3", "beta_xx4", "gam")
  colnames(ci_theta) <- c("est", "SE", "Wald Z", "p-value", "low", "upp")
  
  return(list(
    out = c(tstat, asyV, likval),
    ci_theta = ci_theta,
    asyVar_theta = asyVar_theta
  ))
}


### to calculate the P-value
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(doSNOW))
suppressPackageStartupMessages(library(pracma))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
numCores = 50
nboot = 500


set.seed(010101)
lenini = 50
inipmat = matrix(runif(lenini*5,min = -2,max = 2), ncol = 5)
#######################################################
####### calculate P-value for different initial point
cll <- makeSOCKcluster(numCores)
registerDoSNOW(cll)
parallel::clusterExport(cll, varlist = c("startpar","estpar","boot","asyVar_fun"))
parallel::clusterExport(cl=cll,varlist = c('cotton','nboot','inipmat'))

pb <- txtProgressBar(min=1, max=lenini, style=3) ### set progress bar
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

pval_out <- foreach(i=1:lenini,.options.snow=opts,.combine='rbind') %dopar% {
  bout = boot(cotton, nboot, inipmat[i,])
  pval_boot=mean(abs(bout$bootstat) > abs(bout$tstat) )
  likval_boot =  bout$likelihood_value
  
  asyres=asyVar_fun(cotton,inipmat[i,])
  asyout=asyres$out
  
  
  pval_asy = 1 - pnorm(asyout[1],mean = 0, sd =sqrt(asyout[2]))
  likval_asy = asyout[3]
  
  return(c(pval_boot,pval_asy,likval_boot,likval_asy))
}
close(pb)
stopCluster(cll)

colnames(pval_out) = c("pval_boot","pval_asym", "likeval_boot","likeval_asym")
rownames(pval_out) = NULL


minidx_boot = which.min(pval_out[,3])
minidx_asym = which.min(pval_out[,4])
print(minidx_boot)
print(minidx_asym)

optim_pval <- rbind(pval_out[minidx_boot,],
                    pval_out[minidx_asym,])
rownames(optim_pval) = c("minidx_boot","minidx_asym")
asyres=asyVar_fun(cotton,inipmat[minidx_asym,])
tab_ci_theta = asyres$ci_theta



#### save resulta 
DGP_type_model = "IV_age/results_pvalue"
###
#path <- "D:/20241031GOF-testing/boot-test/"
path <- "P:/MLC/GOF_0131/real_data/new_cotton/"
setwd(path)

sub_path = paste(path,DGP_type_model,sep = "")
if (!dir.exists(sub_path)){
  dir.create(sub_path, recursive = TRUE)  #recursive = TRUE 确保所有父目录也被创建
}

filename_pvals <- paste("res_pvalues.csv", sep="")
full_res_pvals_path <- file.path(sub_path,  filename_pvals)
write.csv(pval_out,file = full_res_pvals_path, row.names = FALSE)
write.csv(optim_pval,file = file.path(sub_path, "optimal_pvalues.csv"), row.names = TRUE)
write.csv(tab_ci_theta,file = file.path(sub_path, "optimal_tab_ci_theta.csv"), row.names = TRUE)

###
> optim_pval
pval_boot pval_asym likeval_boot likeval_asym
minidx_boot     0.776 0.4594713     616.0254     616.0254
minidx_asym     0.776 0.4594713     616.0254     616.0254
> tab_ci_theta
              est        SE     Wald Z      p-value        low         upp
alp      -1.80693984 0.3131874 -5.7695164 7.949932e-09 -2.4207871 -1.19309254
beta_xx2 -0.45107992 0.2102453 -2.1454937 3.191340e-02 -0.8631607 -0.03899917
beta_xx3  0.07207887 0.1999054  0.3605649 7.184247e-01 -0.3197357  0.46389340
beta_xx4 -0.23222870 0.5835981 -0.3979257 6.906850e-01 -1.3760810  0.91162364
gam       1.22543039 1.3795198  0.8883022 3.743782e-01 -1.4784284  3.92928916



# ##--------------------------------
# ##only for plug-in test 
# asy_pval_likvals <- matrix(0, lenini,2)
# for (i in 1:lenini) {
#   asyres=asyVar_fun(cotton,inipmat[i,])
#   asyout=asyres$out
#   pval_asy = 1 - pnorm(asyout[1],mean = 0, sd =sqrt(asyout[2]))
#   asy_likvals = asyout[3]
#   asy_pval_likvals[i,] = c(pval_asy, asy_likvals)
# }
# 
# minidx_asym = which.min(asy_pval_likvals[,2])
# print(minidx_asym) #19
# 
# asyoptim =asyVar_fun(cotton,inipmat[19,])
# tab_ci_theta = asyoptim$ci_theta
# 
# 
# write.csv(tab_ci_theta,file = "D:/20241031GOF-testing/boot-test-v2/2026年重新做的realdata/cotton/results_pvalue/new_optimal_tab_ci_theta.csv", row.names = TRUE)