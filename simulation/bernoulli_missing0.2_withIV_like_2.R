###*##——————————————————————————————————————————————————————##
#  The data generation process: Bernoulli 
#  the missing probability is 0.2
remove(list = ls())
gendat=function(n,nullcase = 0, seed = 0)
{ set.seed(seed)
  xx1 <- rnorm(n, 0, 1)
  xx2 <- rbinom(n,1, 0.5)
  xx3 <- rnorm(n, 1, 1)
  tilr0 = 0
  tilr1 = xx1*xx1*xx2
  tilr2 = 0.5*xx1*xx1
  tilr3 = 0.5*xx1*xx1 + 0.5*xx1*xx1*xx2
  gx4 = xx1*xx1*xx2
  gx5 = xx1*xx1
  gx6 = 0.5*xx1*xx1 + xx1*xx1*xx2
  
  if(nullcase==0){
    p <- c(-1.1,-1.5,-1.5,-0.5); 
    p0 <- p[1]; p1 <- p[2]; p2 <- p[3];p3 <- p[4]
    g <- c(1,-1,-1,-2); 
    g0 <- g[1]; g1 <- g[2]; g2 <- g[3];g3 <- g[4]
    
    ypp =1/(1 + exp(-g0-g1*xx1-g2*xx2-g3*xx3))
    cx= log((1-ypp) + ypp*exp(p3))
    rr <- rbinom(n,1,1/(1+exp(p0+p1*xx1 + p2*xx2 + cx + tilr0)))
    #mean(rr)
  }else if(nullcase==1){
    p <- c(-1.6,-2.3,-2.5,-0.5); 
    p0 <- p[1]; p1 <- p[2]; p2 <- p[3];p3 <- p[4]
    g <- c(1,-1,-1,-2); 
    g0 <- g[1]; g1 <- g[2]; g2 <- g[3];g3 <- g[4]
    
    ypp =1/(1 + exp(-g0-g1*xx1-g2*xx2-g3*xx3))
    cx= log((1-ypp) + ypp*exp(p3))
    rr <- rbinom(n,1,1/(1+exp(p0+p1*xx1 + p2*xx2 + cx + tilr1)))
    #mean(rr)
  }else if(nullcase == 2){
    p <- c(-1.6,-2,-2,-0.5); 
    p0 <- p[1]; p1 <- p[2]; p2 <- p[3];p3 <- p[4]
    g <- c(1,-1,-1,-2); 
    g0 <- g[1]; g1 <- g[2]; g2 <- g[3];g3 <- g[4]
  
    ypp =1/(1 + exp(-g0-g1*xx1-g2*xx2-g3*xx3))
    cx= log((1-ypp) + ypp*exp(p3))
    rr <- rbinom(n,1,1/(1+exp(p0+p1*xx1 + p2*xx2 + cx + tilr2)))
    #mean(rr)
  }else if(nullcase == 3){
    p <- c(-1.6,-1.5,-2,-0.5); 
    p0 <- p[1]; p1 <- p[2]; p2 <- p[3];p3 <- p[4]
    g <- c(1,-1,-1,-2); 
    g0 <- g[1]; g1 <- g[2]; g2 <- g[3];g3 <- g[4]
    
    ypp =1/(1 + exp(-g0-g1*xx1-g2*xx2-g3*xx3))
    cx= log((1-ypp) + ypp*exp(p3))
    rr <- rbinom(n,1,1/(1+exp(p0+p1*xx1 + p2*xx2 + cx + tilr3)))
    #mean(rr)
  }else if(nullcase == 4){
    p <- c(-0.8,-0.5,-1.8,-0.5);
    p0 <- p[1]; p1 <- p[2]; p2 <- p[3];p3 <- p[4]
    g <- c(1,-1,-1,-2); 
    g0 <- g[1]; g1 <- g[2]; g2 <- g[3];g3 <- g[4]
    
    ypp =1/(1 + exp(-g0-g1*xx1-g2*xx2-g3*xx3))
    cx= log( (1-ypp) + ypp*exp(p3 + gx4) )
    rr <- rbinom(n,1,1/(1+exp(p0+p1*xx1 + p2*xx2 + cx)))
    mean(rr)
  }else if(nullcase == 5){
    p <- c(-1,1,-2.5,-0.5);
    p0 <- p[1]; p1 <- p[2]; p2 <- p[3];p3 <- p[4]
    g <- c(1,-1,-1,-2); 
    g0 <- g[1]; g1 <- g[2]; g2 <- g[3];g3 <- g[4]
    
    ypp =1/(1 + exp(-g0-g1*xx1-g2*xx2-g3*xx3))
    cx= log((1-ypp) + ypp*exp(p3 + gx5))
    rr <- rbinom(n,1,1/(1+exp(p0+p1*xx1 + p2*xx2 + cx )))
    mean(rr)
  }else if(nullcase == 6){
    p <- c(-1,-1,-2.5,-0.5) 
    p0 <- p[1]; p1 <- p[2]; p2 <- p[3];p3 <- p[4]
    g <- c(1,-1,-1,-2); 
    g0 <- g[1]; g1 <- g[2]; g2 <- g[3];g3 <- g[4]
    
    ypp =1/(1 + exp(-g0-g1*xx1-g2*xx2-g3*xx3))
    cx= log((1-ypp) + ypp*exp(p3 + gx6))
    rr <- rbinom(n,1,1/(1+exp(p0+p1*xx1 + p2*xx2 + cx)))
    mean(rr)
  }

  yy <- rbinom(n, 1, ypp )
  cbind(yy,xx1,xx2,xx3,rr)
}


startpar<-function(yy,xx1,xx2,xx3,rr){
  # test the significance of the coefficient of IV in f(y|x,R=1)
  ## to estimate xi
  glom = glm(yy[rr==1] ~ xx1[rr==1] + xx2[rr==1] + xx3[rr==1],family=binomial(link = "logit"))
  glom.coef = coef(glom) ##summary(glom)
  olp = glom.coef[1] + glom.coef[2]*xx1 + glom.coef[3]*xx2+ glom.coef[4]*xx3
  hatxi = as.numeric(glom.coef)
  
  ## to estimate theta
  ell <- function(par){
    xixx=hatxi[1]+hatxi[2]*xx1+hatxi[3]*xx2+hatxi[4]*xx3
    px=1/(1+exp(-xixx))
    qx=1/(1+exp(xixx))
    tt=par[1]+par[2]*xx1+par[3]*xx2+log(qx+exp(par[4])*px)
    res=sum(rr*tt + log(1+exp(-tt)))
    return(res)
  }
  
  Grad <- function(par){
    ss=hatxi[1]+hatxi[2]*xx1+hatxi[3]*xx2+hatxi[4]*xx3
    px=1/(1+exp(-ss))
    qx=1/(1+exp(ss))
    tt=par[1]+par[2]*xx1+par[3]*xx2+log(qx+exp(par[4])*px)
  
    tem0= 1-rr - 1/(1+exp(-tt))
    tem1 = sum(tem0)
    tem2 = colSums((matrix(tem0) %*% matrix(1,1,2))*cbind(xx1,xx2))
    tem2 = as.numeric(tem2)
    tem3 = sum(tem0*px*exp(par[4])/(1 + (exp(par[4]) -1)*px))
    c(-tem1,-tem2,-tem3)
  }
  # grad(ell,p)
  # Grad(p)
  hattheta = optim(par = rep(-1,4),fn = ell,gr = Grad,
                   method = "L-BFGS-B",
                   lower = rep(-3,4), upper = rep(3,4))$par
  return(c(as.numeric(hattheta),as.numeric(hatxi)))
}


estpar <- function(yy,xx1,xx2,xx3,rr)
{
  start = startpar(yy,xx1,xx2,xx3,rr)
  ell <- function(par){
    alp = par[1]
    beta= par[2:3]
    gam = par[4]
    ss = par[5] + par[6]*xx1 + par[7]*xx2 + par[8]*xx3
    ppy = 1/(1 + exp(-ss) )
    
    tt = alp + beta[1]*xx1 + beta[2]*xx2 + log(1 + (exp(gam)-1)*ppy )
    tem1 = rr*((1-yy)*ss + log(1 + exp(-ss)))
    tem2 = rr*tt + log(1 + exp(-tt))
    res = sum(tem1+tem2 ) 
    return(res)
  }
  
  Grad <- function(par){
    alp = par[1]
    beta= par[2:3]
    gam = par[4]
    ss = par[5] + par[6]*xx1 + par[7]*xx2 + par[8]*xx3
    px=1/(1+exp(-ss))
    qx=1/(1+exp(ss))
    tt=par[1]+par[2]*xx1+par[3]*xx2+log(qx+exp(par[4])*px)
    
    tem0= 1-rr - 1/(1+exp(-tt))
    tem1 = sum(tem0)
    tem2 = colSums((matrix(tem0) %*% matrix(1,1,2))*cbind(xx1,xx2))
    tem2 = as.numeric(tem2)
    tem3 = sum(tem0*px*exp(par[4])/(1 + (exp(par[4]) -1)*px))
    tt4 = rr*(yy - px) + tem0*(exp(gam)-1)*px*qx/(1+(exp(gam)-1)*px)
    tem4=colSums((matrix(tt4) %*% matrix(1,1,4))*cbind(1,xx1,xx2,xx3))
    tem4=as.numeric(tem4)
    c(-tem1,-tem2,-tem3,-tem4)
  }
  # ell(pp)
  # Grad(pp)
  # grad(ell,pp)
  epar = optim(par = start,fn = ell,gr = Grad,method = "L-BFGS-B",
               lower = rep(-3,8), upper = rep(3,8))$par
  ess = epar[5] + epar[6]*xx1+epar[7]*xx2 + epar[8]*xx3
  ppy = 1/(1 + exp(-ess) )
  
  nlp = epar[1] + epar[2]*xx1+  epar[3]*xx2 + log(1 -ppy  + exp(epar[4])*ppy)
  pix=1/(1+exp(nlp))
  
  HH = (rr-pix)^2 - pix*(1-pix)
  nn=length(rr)
  tstat=mean(HH)*sqrt(nn)
  return(list(hatxi = as.numeric(epar[5:8]),
              hattheta=as.numeric(epar[1:4]),tstat=tstat))
}


boot<- function(data,nboot = 100)
{
  yy=data[,1]
  xx1=data[,2]
  xx2=data[,3]
  xx3=data[,4]
  rr=data[,5]
  nn=length(rr)
  out=estpar(yy,xx1,xx2,xx3,rr)
  hatxi=out$hatxi
  hattheta=out$hattheta
  tstat=out$tstat
  
  bootstat=c()
  
  for(i in 1:nboot)
  { 
    set.seed(137*i+ 37)
    indx=sample(1:nn,nn,replace=T)
    bxx1=xx1[indx]
    bxx2=xx2[indx]
    bxx3=xx3[indx]
    bhatolp = hatxi[1]+hatxi[2]*bxx1 + hatxi[3]*bxx2 + hatxi[4]*bxx3
    bhatppy = 1/(1 + exp(-bhatolp) )
   
    bnlp=hattheta[1]+hattheta[2]*bxx1 +hattheta[3]*bxx2+log(1+(exp(hattheta[4])-1)*bhatppy)
    bp=1/(1+exp(bnlp))
    brr=rbinom(nn,1,bp)
    byy= rbinom(nn,1,bhatppy)
    bout=estpar(byy,bxx1,bxx2,bxx3,brr)
    bootstat=c(bootstat,bout$tstat)
    #print(i)
  }
  c(tstat,sd(bootstat),quantile(abs(bootstat),c(0.9,0.95,0.99)))
}



### directly estimate asymptotic variance
asyVar_fun <- function(data){
  yy=data[,1]
  xx1=data[,2]
  xx2=data[,3]
  xx3=data[,4]
  rr=data[,5]
  nn=length(rr)
  out=estpar(yy,xx1,xx2,xx3,rr)
  hatxi=out$hatxi
  hattheta=out$hattheta
  tstat=out$tstat
  
  hatolp = hatxi[1]+hatxi[2]*xx1 + hatxi[3]*xx2 + hatxi[4]*xx3
  hatppy = 1/(1 + exp(-hatolp) )
  tt=hattheta[1]+hattheta[2]*xx1 +hattheta[3]*xx2+log(1+(exp(hattheta[4])-1)*hatppy)
  pix=1/(1+exp(tt))
  
  rx = cbind(xx1,xx2)
  dim_rx = ncol(rx)
  cX = cbind(1,xx1,xx2,xx3)
  colnames(cX) = NULL
  
  alp =  hattheta[1]
  beta = hattheta[2:(1+dim_rx)]
  gam = hattheta[4]
  
  Rat10 = exp(gam + hatolp)/(1 + exp(hatolp) ) ## A11/A01
  Delta = lapply(1:nn, function(i){
    term1 = cX[i,]*(Rat10[i] - hatppy[i])
    term = c(1,rx[i,],Rat10[i], term1)
    return(matrix(term))
  })
  Delta =  do.call("cbind", Delta)##dim_theta*nn
  
  
  Lambda = lapply(1:nn, function(i){
    term1 = cX[i,]*(yy[i]- hatppy[i])
    term = c(rep(0, (dim_rx + 2)), term1)
    return(matrix(term))
  })
  Lambda  =  do.call("cbind", Lambda )##dim_theta*nn
  
  dim_theta = 4+4
  Rmat =  matrix(rep(1,dim_theta)) %*% matrix(rr, nrow  =1) ##dim_theta*nn
  Pimat = matrix(rep(1,dim_theta)) %*% matrix(pix, nrow  =1)
  coEmat = matrix(rep(1,dim_theta))  %*% matrix((1+2*rr - 4*pix)*pix*(1-pix), nrow  =1)
  Smat = Rmat*Lambda + (Pimat - Rmat)*Delta ## dim_theta*nn
  
  HH = (rr - pix)^2 - pix*(1 - pix) 
  EH = matrix(rowMeans(coEmat*Delta)) ##dim_theta*1
  JJ = Smat %*% t(Smat)/nn ##dim_theta*dim_theta
  inv_JJ = solve(JJ)
  KK = HH + t(t(EH)%*%inv_JJ%*%Smat)
  asyV = var(KK)
  return(c(tstat,asyV))
}

###******************************* adjusted bootstrap procedure ********************###
###*
adj_boot <- function(data,nboot = 100)
{
  yy=data[,1]
  xx1=data[,2]
  xx2=data[,3]
  xx3=data[,4]
  rr=data[,5]
  nn=length(rr)
  out=estpar(yy,xx1,xx2,xx3,rr)
  hatxi=out$hatxi
  hattheta=out$hattheta
  tstat=out$tstat
  
  bootstat=c()
  for(i in 1:nboot)
  { set.seed(137*i+ 37)
    indx=sample(1:nn,nn,replace=T)
    bxx1=xx1[indx]
    bxx2=xx2[indx]
    bxx3=xx3[indx]
    byy=yy[indx]
    brr=rr[indx]
    bout=estpar(byy,bxx1,bxx2,bxx3,brr)
    bootstat=c(bootstat,bout$tstat)
  }
  c(tstat,sd(bootstat),quantile(bootstat,c(0.025,0.975)))
}





### simulation process
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(pbapply))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(doSNOW))
suppressPackageStartupMessages(library(pracma))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))



## data generating random seed
path <- "simulation"
rs_file <- file.path(path, "gof_simu_rs.txt")
rs_read <- scan(rs_file, what = numeric())


for(n in c(500,1000,2000,4000)){
  print(n)
  NS = 1000
  nboot = 500
  numCores = 70
  ##********** parameters setting *************##
  DGP_type = "Bernoulli"
  DGP_type_model = "Bernoulli_boot_v2_2"
  miss_prob = 0.2
  
  
  ###
  #path <- "D:/20241031GOF-testing/boot-test/"
  
  sub_path = paste(path,DGP_type_model,sep = "")
  if (!dir.exists(sub_path)){
    dir.create(sub_path, recursive = TRUE)  #recursive = TRUE 确保所有父目录也被创建
  }
  
  subsub_path = paste(sub_path,"/missProb",miss_prob,"_n",n,"_NS",NS,sep = "")
  if (!dir.exists(subsub_path)){
    dir.create(subsub_path, recursive = TRUE)  
  }

  rej_tab <- matrix(0, nrow = 3, ncol = 5)
  rownames(rej_tab) = c("mixboot_variance","mixboot_percentile","closedasym_variance")
  colnames(rej_tab) = paste("nullcase=", c(0:4))
  

  for (tnullcase in c(0,2,3,5,6)) {##
    if(tnullcase == 0) {
      tNS = 10*NS
      trs = rs_read[1:tNS]
    }else{
      tNS = NS 
      trs = rs_read[1:tNS]
    }
    
    
    cll <- makeSOCKcluster(numCores)
    registerDoSNOW(cll)
    parallel::clusterExport(cll, varlist = c("gendat","startpar","estpar","boot"))
    parallel::clusterExport(cl=cll,varlist = c('n','tNS','nboot','trs','tnullcase'))


    pb <- txtProgressBar(min=1, max=tNS, style=3) ### set progress bar
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)

    ttout <- foreach(i=1:tNS,.options.snow=opts,
                     .combine='rbind') %dopar% {
                       data=gendat(n, nullcase = tnullcase,seed = trs[i])
                       out=boot(data, nboot)
                       return(out)
                     }
    close(pb)
    stopCluster(cll)

    rownames(ttout) = NULL
    assign(paste0("output_nullcase",tnullcase), ttout)


    tstat=ttout[,1]
    se=ttout[,2]
    c2=ttout[,4]

    rej_tab[1,tnullcase+1] = mean(abs(tstat/se)>1.96)
    rej_tab[2,tnullcase+1] = mean(abs(tstat)>c2)
    print(rej_tab)
    #
    
    #####################################################
    ## others: directly estimate asymptotic variance
    cll <- makeSOCKcluster(numCores)
    registerDoSNOW(cll)
    parallel::clusterExport(cll, varlist = c("gendat","estpar","asyVar_fun"))
    parallel::clusterExport(cl=cll,varlist = c('n','tNS','trs','tnullcase'))
    
    pb <- txtProgressBar(min=1, max=tNS, style=3) ### set progress bar
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
    
    ttres <- foreach(i=1:tNS,.options.snow=opts, 
                     .combine='rbind') %dopar% {
                       data=gendat(n, nullcase = tnullcase,seed = trs[i])
                       res=asyVar_fun(data)
                       return(res)
                     }
    close(pb)
    stopCluster(cll)
    
    rownames(ttres) = NULL
    assign(paste0("resput_nullcase",tnullcase), ttres)
    
    tstat=ttres[,1]
    asyV=ttres[,2]
    rej_tab[3,tnullcase+1] = mean(abs(tstat/sqrt(asyV))>1.96)
    print(rej_tab)

  }
  filename_Rej_Tab <- paste("rej_TAB_",DGP_type_model,"_missProb",miss_prob,"_n",n,".csv", sep="")
  full_Rej_Tab_path <- file.path(subsub_path,  filename_Rej_Tab)
  write.csv(rej_tab,file = full_Rej_Tab_path)
  
  ######
  filname_image <- paste(DGP_type_model,"_missProb",miss_prob,"_n",n,".RData", sep="")
  full_image_path <- file.path(subsub_path, filname_image)
  save.image(full_image_path)
}





