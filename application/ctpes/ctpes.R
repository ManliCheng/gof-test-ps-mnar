remove(list = ls())
## dataset details:
## Y: teacher’s report of the psychopathology status of the student
##      Y = 0: normal
##      Y = 1: borderline or clinical psychopathology; Y = -1: missing
##covariates: father figure present(X1), physical health(X2), Parents’ report(X3)
##      X1 = 0:father figure present
##      X2 = 0:no health problems
##      X3 = 0:normal(IV)
##Model:
## logit(p(Y=1|X1,X2))
## logit(P(R=1|Y,X2))


data_ctpes = matrix(c(rep(c(1,0,-1,1),27),
                      rep(c(1,0,-1,0),100),
                      rep(c(1,0,1,1),13),
                      rep(c(1,0,1,0),16),
                      rep(c(1,0,0,1),13),
                      rep(c(1,0,0,0),88),
                      rep(c(1,1,-1,1),33),
                      rep(c(1,1,-1,0),85),
                      rep(c(1,1,1,1),25),
                      rep(c(1,1,1,0),15),
                      rep(c(1,1,0,1),28),
                      rep(c(1,1,0,0),71),
                      rep(c(0,0,-1,1),57),
                      rep(c(0,0,-1,0),407),
                      rep(c(0,0,1,1),15),
                      rep(c(0,0,1,0),76),
                      rep(c(0,0,0,1),38),
                      rep(c(0,0,0,0),471),
                      rep(c(0,1,-1,1),74),
                      rep(c(0,1,-1,0),278),
                      rep(c(0,1,1,1),56),
                      rep(c(0,1,1,0),45),
                      rep(c(0,1,0,1),83),
                      rep(c(0,1,0,0),372)), byrow = T,ncol = 4)
#rr = (is.na(data_ctpes[,3]))+0
ctpes = cbind(data_ctpes[,3],data_ctpes[,1:2],data_ctpes[,4],(data_ctpes[,3] != -1)+0)
yy = ctpes[,1]
xx1 = ctpes[,2]
xx2 = ctpes[,3]
xx3 = ctpes[,4]
rr = ctpes[,5]



###1. goodness of fit to find the suitable regression model
library(rms)
ryy = yy[rr==1]
rxx1= xx1[rr==1]
rxx2= xx2[rr==1]
rxx3= xx3[rr==1]
o1 = lrm(ryy ~rxx1+rxx2,x = T, y = T )
out= resid(o1,"gof")
o1
out

##Model:
## logit(p(Y=1|X1,X2))
## logit(P(R=1|Y,X2))

# o1
# Logistic Regression Model
# 
# lrm(formula = ryy ~ rxx1 + rxx2, x = T, y = T)
# 
# Model Likelihood      Discrimination    Rank Discrim.    
# Ratio Test             Indexes          Indexes    
# Obs          1425    LR chi2     14.36      R2       0.016    C       0.565    
# 0           1164    d.f.            2     R2(2,1425)0.009    Dxy     0.129    
# 1            261    Pr(> chi2) 0.0008    R2(2,639.6)0.019    gamma   0.192    
# max |deriv| 1e-05                           Brier    0.148    tau-a   0.039    
# 
# Coef    S.E.   Wald Z Pr(>|Z|)
# Intercept -1.7372 0.1070 -16.24 <0.0001 
# rxx1       0.5419 0.1607   3.37 0.0007  
# rxx2       0.2465 0.1380   1.79 0.0740  
# 
# > out
# Sum of squared errors     Expected value|H0                    SD                     Z                     P 
# 210.89765164          210.92807479            0.07857282           -0.38719688            0.69861046 





### 2. to estimate
startpar <- function(yy, xx1, xx2, xx3, rr, inip = rep(1, 3)) {
  
  ## 1) estimate xi: NOW use (1, xx1, xx2)
  glom <- glm(
    yy[rr == 1] ~ xx1[rr == 1] + xx2[rr == 1],
    family = binomial(link = "logit")
  )
  glom.coef <- coef(glom)
  
  ## hatxi: (intercept, xx1, xx2)
  hatxi <- as.numeric(glom.coef)
  
  ## 2) estimate theta: unchanged
  ell <- function(par) {
    
    ## xi linear predictor: xx1 + xx2 (NO xx3)
    xixx <- hatxi[1] + hatxi[2] * xx1 + hatxi[3] * xx2
    px   <- 1 / (1 + exp(-xixx))
    qx   <- 1 / (1 + exp(xixx))
    
    ## theta unchanged: alp + beta*xx2 + log(...)
    tt <- par[1] + par[2] * xx2 + log(qx + exp(par[3]) * px)
    
    sum(rr * tt + log(1 + exp(-tt)))
  }
  
  Grad <- function(par) {
    
    ## xi linear predictor: xx1 + xx2
    ss <- hatxi[1] + hatxi[2] * xx1 + hatxi[3] * xx2
    px <- 1 / (1 + exp(-ss))
    qx <- 1 / (1 + exp(ss))
    
    ## theta unchanged
    tt <- par[1] + par[2] * xx2 + log(qx + exp(par[3]) * px)
    
    tem0 <- 1 - rr - 1 / (1 + exp(-tt))
    
    g_alp  <- -sum(tem0)
    g_beta <- -sum(tem0 * xx2)
    g_gam  <- -sum(tem0 * px * exp(par[3]) /
                     (1 + (exp(par[3]) - 1) * px))
    
    c(g_alp, g_beta, g_gam)
  }
  
  hattheta <- optim(
    par    = inip,
    fn     = ell,
    gr     = Grad,
    method = "L-BFGS-B",
    lower  = rep(-5, 3),
    upper  = rep(5, 3)
  )$par
  
  ## return: theta(3) + xi(3) => length 6
  c(as.numeric(hattheta), as.numeric(hatxi))
}



estpar <- function(yy, xx1, xx2, xx3, rr, inip = rep(1, 3)) {
  
  start <- startpar(yy, xx1, xx2, xx3, rr, inip)  # theta(3)+xi(3) => length 6
  
  ell <- function(par){
    alp  <- par[1]
    beta <- par[2]   # coefficient for xx2 only
    gam  <- par[3]
    
    ## xi model: intercept + xx1 + xx2 (NO xx3)
    ss  <- par[4] + par[5] * xx1 + par[6] * xx2
    ppy <- 1 / (1 + exp(-ss))
    
    ## theta: alp + beta*xx2 + log(1 + (exp(gam)-1)*ppy)
    tt <- alp + beta * xx2 + log(1 + (exp(gam) - 1) * ppy)
    
    tem1 <- rr * ((1 - yy) * ss + log(1 + exp(-ss)))
    tem2 <- rr * tt + log(1 + exp(-tt))
    
    sum(tem1 + tem2)
  }
  
  Grad <- function(par){
    alp  <- par[1]
    beta <- par[2]
    gam  <- par[3]
    
    ## xi model: intercept + xx1 + xx2
    ss <- par[4] + par[5] * xx1 + par[6] * xx2
    px <- 1 / (1 + exp(-ss))
    qx <- 1 / (1 + exp(ss))
    
    tt <- alp + beta * xx2 + log(qx + exp(gam) * px)
    
    tem0 <- 1 - rr - 1 / (1 + exp(-tt))
    
    ## gradients for theta: (alp, beta_xx2, gam)
    g_alp  <- -sum(tem0)
    g_beta <- -sum(tem0 * xx2)
    g_gam  <- -sum(tem0 * px * exp(gam) / (1 + (exp(gam) - 1) * px))
    
    ## gradients for xi: (1, xx1, xx2)
    tt4 <- rr * (yy - px) + tem0 * (exp(gam) - 1) * px * qx / (1 + (exp(gam) - 1) * px)
    g_xi <- -as.numeric(colSums((matrix(tt4) %*% matrix(1, 1, 3)) * cbind(1, xx1, xx2)))
    
    c(g_alp, g_beta, g_gam, g_xi)
  }
  
  eres <- optim(
    par    = start,
    fn     = ell,
    gr     = Grad,
    method = "L-BFGS-B",
    lower  = rep(-5, 6),
    upper  = rep(5, 6),hessian=T
  )
  
  epar <- eres$par
  eval <- eres$value
  
  ## recompute ppy with fitted xi (uses xx1, xx2; NO xx3)
  ess <- epar[4] + epar[5] * xx1 + epar[6] * xx2
  ppy <- 1 / (1 + exp(-ess))
  
  ## theta: alp + beta*xx2 + log(1 - ppy + exp(gam)*ppy)
  nlp <- epar[1] + epar[2] * xx2 + log(1 - ppy + exp(epar[3]) * ppy)
  pix <- 1 / (1 + exp(nlp))
  
  HH <- (rr - pix)^2 - pix * (1 - pix)
  nn <- length(rr)
  tstat <- mean(HH) * sqrt(nn)
  
  return(list(
    hatxi    = as.numeric(epar[4:6]),  # (xi0, xi1, xi2) for (1, xx1, xx2)
    hattheta = as.numeric(epar[1:3]),  # (alp, beta_xx2, gam)
    tstat    = tstat,
    likelihood_value = eval,hes=eres$hessian
  ))
}


boot <- function(data, nboot = 100, inip = rep(1, 3)) {
  
  yy  <- data[,1]
  xx1 <- data[,2]
  xx2 <- data[,3]
  xx3 <- data[,4]
  rr  <- data[,5]
  nn  <- length(rr)
  
  out <- estpar(yy, xx1, xx2, xx3, rr, inip)
  hatxi    <- out$hatxi      # length 3: (1, xx1, xx2)
  hattheta <- out$hattheta   # length 3: (alp, beta_xx2, gam)
  tstat    <- out$tstat
  likval   <- out$likelihood_value
  
  bootstat <- c()
  
  for (i in 1:nboot) {
    set.seed(137*i + 837)
    
    indx <- sample(1:nn, nn, replace = TRUE)
    bxx1 <- xx1[indx]
    bxx2 <- xx2[indx]
    bxx3 <- xx3[indx]
    
    ## xi part: uses xx1 and xx2 (NO xx3)
    bhatolp <- hatxi[1] + hatxi[2] * bxx1 + hatxi[3] * bxx2
    bhatppy <- 1 / (1 + exp(-bhatolp))
    
    ## theta part: alp + beta*xx2 + log(...)
    bnlp <- hattheta[1] + hattheta[2] * bxx2 +
      log(1 + (exp(hattheta[3]) - 1) * bhatppy)
    
    bp  <- 1 / (1 + exp(bnlp))
    brr <- rbinom(nn, 1, bp)
    byy <- rbinom(nn, 1, bhatppy)
    
    bout <- estpar(byy, bxx1, bxx2, bxx3, brr, inip)
    bootstat <- c(bootstat, bout$tstat)
  }
  
  list(tstat = tstat, bootstat = bootstat, likelihood_value = likval)
}

### directly estimate asymptotic variance
asyVar_fun <- function(data, inip = rep(1, 3)) {
  
  yy  <- data[,1]
  xx1 <- data[,2]
  xx2 <- data[,3]
  xx3 <- data[,4]
  rr  <- data[,5]
  nn  <- length(rr)
  
  out <- estpar(yy, xx1, xx2, xx3, rr, inip)
  hatxi    <- out$hatxi       # length 3: (1, xx1, xx2)
  hattheta <- out$hattheta    # length 3: (alp, beta_xx2, gam)
  tstat    <- out$tstat
  likval   <- out$likelihood_value
  
  ## xi linear predictor: uses xx1 and xx2 (NO xx3)
  hatolp <- hatxi[1] + hatxi[2] * xx1 + hatxi[3] * xx2
  hatppy <- 1 / (1 + exp(-hatolp))
  
  ## theta: REMOVE xx1, keep xx2
  tt  <- hattheta[1] + hattheta[2] * xx2 + log(1 + (exp(hattheta[3]) - 1) * hatppy)
  pix <- 1 / (1 + exp(tt))
  
  ## theta covariates (only xx2 now)
  rx <- cbind(xx2)
  dim_rx <- ncol(rx)   # = 1
  
  ## covariates for xi score part: (1, xx1, xx2)
  cX_xi <- cbind(1, xx1, xx2)
  colnames(cX_xi) <- NULL
  
  alp  <- hattheta[1]
  beta <- hattheta[2]   # scalar
  gam  <- hattheta[3]
  
  Rat10 <- exp(gam + hatolp) / (1 + exp(hatolp))  ## A11/A01
  
  ## Delta: dim_theta x nn, dim_theta = (dim_rx+2) + 3 = 3 + 3 = 6
  Delta <- lapply(1:nn, function(i){
    term1 <- cX_xi[i,] * (Rat10[i] - hatppy[i])  # length 3
    term  <- c(1, rx[i,], Rat10[i], term1)       # length 1+1+1+3 = 6
    matrix(term)
  })
  Delta <- do.call("cbind", Delta)  ## 6 x nn
  
  Lambda <- lapply(1:nn, function(i){
    term1 <- cX_xi[i,] * (yy[i] - hatppy[i])     # length 3
    term  <- c(rep(0, (dim_rx + 2)), term1)      # first 3 zeros, then 3 -> 6
    matrix(term)
  })
  Lambda <- do.call("cbind", Lambda)  ## 6 x nn
  
  dim_theta <- 6
  
  Rmat   <- matrix(rep(1, dim_theta)) %*% matrix(rr,  nrow = 1)   ## 6 x nn
  Pimat  <- matrix(rep(1, dim_theta)) %*% matrix(pix, nrow = 1)   ## 6 x nn
  coEmat <- matrix(rep(1, dim_theta)) %*% matrix((1 + 2*rr - 4*pix) * pix * (1 - pix), nrow = 1)
  
  Smat <- Rmat * Lambda + (Pimat - Rmat) * Delta  ## 6 x nn
  
  HH <- (rr - pix)^2 - pix * (1 - pix)
  
  EH     <- matrix(rowMeans(coEmat * Delta))   ## 6 x 1
  JJ     <- Smat %*% t(Smat) / nn              ## 6 x 6
  inv_JJ <- solve(JJ)
  
  KK   <- HH + t(t(EH) %*% inv_JJ %*% Smat)
  asyV <- var(KK)
  
  ## Smat: dim_theta x nn, columns are samples
  SigmaS <- cov(t(Smat))               # 6 x 6
  asyVar_theta <- inv_JJ %*% SigmaS %*% inv_JJ
  
  ## CI for theta only: theta has dim_rx+2 = 3 components (alp, beta_xx2, gam)
  ci_theta <- matrix(0, dim_rx + 2, 6)
  for (i in 1:(dim_rx + 2)) {
    SE = sqrt(diag(solve(out$hes)))[i]
    waldZ = hattheta[i]/SE
    pvalue <- 2 * (1 - pnorm(abs(waldZ)))
    ci_theta[i,] <- c(
      hattheta[i],SE,waldZ,pvalue,
      hattheta[i] - 1.96 *SE,
      hattheta[i] + 1.96 *SE)
  }
  rownames(ci_theta) <- c("alp", "beta", "gam")
  colnames(ci_theta) <- c("est","SE","Wald Z","p-value","low", "upp")
  
  return(list(
    out = c(tstat, asyV, likval),
    ci_theta = ci_theta
  ))
}


######################################################################################
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



set.seed(10010111)
lenini =50 
inipmat = matrix(runif(lenini*3,min = -2,max = 2), ncol = 3)
######################################################
####### calculate P-value for different initial point
cll <- makeSOCKcluster(numCores)
registerDoSNOW(cll)
parallel::clusterExport(cll, varlist = c("startpar","estpar","boot","asyVar_fun"))
parallel::clusterExport(cl=cll,varlist = c('ctpes','nboot','inipmat'))

pb <- txtProgressBar(min=1, max=lenini, style=3) ### set progress bar
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

pval_out <- foreach(i=1:lenini,.options.snow=opts,.combine='rbind') %dopar% {
  bout = boot(ctpes, nboot, inipmat[i,])
  pval_boot=mean(abs(bout$bootstat) > abs(bout$tstat) )
  likval_boot =  bout$likelihood_value
  
  asyres=asyVar_fun(ctpes,inipmat[i,])
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
# print(minidx_boot)
# print(minidx_asym)
# print(minidx_boot)
# [1] 38
# > print(minidx_asym)
# [1] 38



optim_pval <- rbind(pval_out[minidx_boot,],
                    pval_out[minidx_asym,])
rownames(optim_pval) = c("minidx_boot","minidx_asym")
asyres=asyVar_fun(ctpes,inipmat[minidx_asym,])
tab_ci_theta = asyres$ci_theta


#### save resulta 
DGP_type_model = "results_pvalue"
###
path <- file.path("application", "ctpes")
sub_path <- file.path(path, DGP_type_model)

sub_path = paste(path,DGP_type_model,sep = "")
if (!dir.exists(sub_path)){
  dir.create(sub_path, recursive = TRUE)  #recursive = TRUE 确保所有父目录也被创建
}

filename_pvals <- paste("res_pvalues.csv", sep="")
full_res_pvals_path <- file.path(sub_path,  filename_pvals)
write.csv(pval_out,file = full_res_pvals_path, row.names = FALSE)
write.csv(optim_pval,file = file.path(sub_path, "optimal_pvalues.csv"), row.names = TRUE)
write.csv(tab_ci_theta,file = file.path(sub_path, "optimal_tab_ci_theta.csv"), row.names = TRUE)

