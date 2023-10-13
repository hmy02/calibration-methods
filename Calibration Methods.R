#' This is some description of this function.
#' @title Calibration Methods
#' 
#' @description There are six calibration methods in this package. L2, KO, LS, SGP, OGP and PK.
#' 
#' @details you can use this function to caculate estimation of calibration parameters, Standard deviation, posterior samples or bootstrap samples and confidence interval.
#' 
#' @param D D is the design points, which is a matrix.
#' @param zeta zeta is the ture process, which is a known function and related to design points.
#' @param y.s y.s is the computer model, which is a known function and related to design and calibration parameters.
#' @param tau2 tau2 is the variance of random error, which is a number.
#' @param lower.x lower.x is the lower bound of design points, which is an vector.
#' @param upper.x upper.x is the upper bound of design points, which is an vector.
#' @param lower.th lower.th is the lower bound of calibration parameters, which is an vector.
#' @param upper.th upper.th is the upper bound of calibration parameters, which is an vector.
#' @param m m is the number of bootstrap, the default is 1000. When the method is "L2" or "LS", it can be used.
#' @param B B is the number of burn-in samples, the default is 1000. When the method is "KO", "SGP" or "OGP", it can be used.
#' @param M M is the number of postrtior samples, the default is 1000. When the method is "KO", "SGP" or "OGP", it can be used.
#' @param thinning thinning is the interval times of sampling posterior samples, the default is 5. When the method is "KO", "SGP" or "OGP", it can be used.
#' @param jump jump is the step size, which is a number or an vector. When the method is "KO", "SGP" or "OGP", it can be used.
#' @param method method contains "TV", "L2", "LS", "KO", "SGP", "OGP". "TV" means the true value of calibration parameters.
#' @param confidence confidence is a logistic variable. "T" is the given confidence interval, "F" means not. When the method is "L2" or "LS", it can be used.
#' @return a list containing:
#' @return theta.hat: the estimation of calibration parameters.
#' @return 
#' @examples x=1;f(x)
#' lx = 0; ux = 1
#' lth = 0; uth = 3  
#' D =  as.matrix(seq(lx,ux,length=10))
#' zeta<-function(x) 5*x*cos(15*x/2)+5*x  
#' y.s<-function(x,th) sin(5*th*x)+5*x 
#' tau2 = 0.1^2; jump=rep(0.01,7)
#' y.p=NULL
#' f1 = fun(D,y.p,zeta,y.s,tau2,lx,ux,lth,uth,m=5,B=1000,M=1000,thinning=5,jump,method="TV",confidence=T)  
#' f2 = fun(D,zeta,y.s,tau2,lx,ux,lth,uth,m=500,B=1000,M=1000,thinning=5,jump,method="LS",confidence=T)  
#' f3 = fun(D,zeta,y.s,tau2,lx,ux,lth,uth,m=5,B=1000,M=1000,thinning=5,jump,method="OGP",confidence=T)  
#' @export

set.seed(1234)
library(cubature) ##Computing integral
library(lava) ##Computing trace
library(mcGlobaloptim) 
library(DiceKriging) ##Fitting Kringing Model
library(nloptr) ##Optimization
library(geoR)
library(lhs)
library(numDeriv)
library(RobustCalibration)
library(MASS)

fun<-function(D,y.p,zeta,y.s,tau2,lower.x,upper.x,lower.th,upper.th,m=1000,B=1000,M=1000,thinning=5,jump,method,confidence){
  p = length(lower.x); q = length(lower.th)
  n = nrow(D)
  if(is.null(y.p) == T){
    y.p = apply(D,1,zeta)+rnorm(n,0,sqrt(tau2))
  }else{
    y.p=y.p
  }
  ys.D<-function(th) apply(D,1,function(x) y.s(x,th))
  
  norm0<-function(x1,x2) (x1-x2)^2
  mat0<-function(x) matern(x, phi=1/2, kappa=5/2)
  R<-function(da1,da2){
    term = matrix(0,nrow(da1),nrow(da2))
    for(i in 1:length(lower.x)){
      term = term + outer(as.vector(da1[,i]),as.vector(da2[,i]),"norm0")
    }
    R0 = sqrt(term)
    R = matrix(sapply(R0,mat0),nrow(da1),nrow(da2))
    return(R)
  }
  
  if(method == "TV"){######True Value######
    model.tv<-function(da,th) (zeta(da)-y.s(da,th))^2
    tv.loss<-function(th) adaptIntegrate(model.tv,lower.x,upper.x,th=th)$integral
    theta.star = directL(tv.loss, lower.th, upper.th,
                         nl.info = TRUE, control=list(xtol_rel=1e-8, maxeval=1000))$par
    return(theta.star=theta.star)
  }
  else if(method == "L2"){######L2 Calibration######
    #Step1: Fitting an non-parametic model
    #L2.non<-function(x,data,y.data){
    #  E.p=as.matrix(dist(data, diag=T, upper=T)); theta.p=3
    #  R.p=exp(-theta.p*E.p^2)
    #  basis.p<-function(h) exp(-theta.p*sum(h^2))
    #  r.p<-function(x){
    #    A=t(t(data)-x)
    #    vec=apply(A,1,basis.p)
    #    return(vec)
    #  }
    #  coef.p=solve(R.p+1e-6*diag(n),y.data)
    #  return(t(r.p(x))%*%coef.p)
    #}
    model.L2.hat = rgasp(design=D,response=y.p)
    zeta.hat<-function(x){
      new = t(as.matrix(x))
      zeta.pre = predict(model.L2.hat,testing_input=new)$mean#newdata=new,type="UK",checkNames=F)$mean
      return(zeta.pre)
    }
    #zeta.hat<-function(x) L2.non(x,D,y.p)
    zeta.hat.D = as.matrix(apply(D,1,zeta.hat))
    
    #Step2: L2 Estimation
    model.L2<-function(da,th) (zeta.hat(da)-y.s(da,th))^2
    L2.loss<-function(th) adaptIntegrate(model.L2,lower.x,upper.x,th=th)$integral
    theta.L2 = directL(L2.loss, lower.th, upper.th,
                       nl.info = TRUE, control=list(xtol_rel=1e-8, maxeval=1000))$par
    
    #Step3: 95% Confidence Interval
    if(confidence == T){
      L2.boot = matrix(0,m,q)
      for(i in 1:m){
        ind.L2 = sample(n,n,replace=T); sam.L2 = as.matrix(D[ind.L2,])
        y.sam.L2 = as.matrix(y.p[ind.L2])#apply(sam.L2,1,zeta)+rnorm(n,0,sqrt(tau2)))
        model.sam = rgasp(design=sam.L2,response=y.sam.L2)
        zeta.sam<-function(x){
          new = t(as.matrix(x))
          zeta.pre = predict(model.sam,testing_input=new)$mean#newdata=new,type="UK",checkNames=F)$mean
          return(zeta.pre)
        }
        #zeta.sam<-function(x) L2.non(x,sam.L2,y.sam.L2)
        model.sam.L2<-function(da,th) (zeta.sam(da)-y.s(da,th))^2
        sam.L2.loss<-function(th) adaptIntegrate(model.sam.L2,lower.x,upper.x,th=th)$integral
        L2.boot[i,] = directL(sam.L2.loss,lower.th,upper.th,control=list(xtol_rel=1e-8, maxeval=1000))$par
      }
      sd.L2 = apply(L2.boot,2,sd)
      lwr.L2=theta.L2-sd.L2*qnorm(0.975); hwr.L2=theta.L2+sd.L2*qnorm(0.975)
      con.L2 = rbind(lwr.L2,hwr.L2)
      result.L2 = list(theta.hat=theta.L2, sd.hat=sd.L2, post.samples=L2.boot, con.hat=con.L2)  
      return(result.L2)
    }else{
      result.L2 = list(theta.hat=theta.L2)
      return(result.L2) 
    }}
  else if(method == "LS"){######LS Calibration######
    theta.LS = directL(function(th) sum((y.p-ys.D(th))^2)/n, lower.th, upper.th,
                       nl.info = TRUE, control=list(xtol_rel=1e-8, maxeval=1000))$par
    if(confidence == T){
      LS.boot = matrix(0,m,q)
      for(i in 1:m){
        ind.LS = sample(n,n,replace=T); sam.LS = as.matrix(D[ind.LS,])
        y.sam.LS = as.matrix(y.p[ind.LS])#apply(sam.LS,1,zeta)+rnorm(n,0,sqrt(tau2)))
        ys.sam.LS = function(th0) as.matrix(apply(sam.LS,1,function(x) y.s(x,th0)))
        model.sam.LS<-function(th) sum((y.sam.LS-ys.sam.LS(th))^2)/n
        LS.boot[i,] = directL(model.sam.LS,lower.th,upper.th,control=list(xtol_rel=1e-8, maxeval=1000))$par
      }
      sd.LS = apply(LS.boot,2,sd)
      lwr.LS=max(theta.LS-sd.LS*qnorm(0.975),lower.th);hwr.LS=min(theta.LS+sd.LS*qnorm(0.975),upper.th)
      con.LS = rbind(lwr.LS,hwr.LS)
      result.LS = list(theta.hat=theta.LS, sd.hat=sd.LS, post.samples=LS.boot, con.hat=con.LS) 
      return(result.LS)
    }else{
      result.LS = list(theta.hat=theta.LS)
      return(result.LS) 
    }}
  else if(method == "KO"){######KO Calibration######
    R.ko = R(D,D)+1e-6*diag(n)
    S.ko <-function(sigma2) sigma2*R.ko+tau2*diag(n)
    A.ko<-function(th,sigma2) sigma2*R.ko%*%solve(S.ko(sigma2))%*%(y.p-ys.D(th))
    b.ko<-function(sigma2) sigma2*R.ko-sigma2^2*R.ko%*%solve(S.ko(sigma2))%*%R.ko
    B.ko<-function(sigma2) 0.5*(b.ko(sigma2)+t(b.ko(sigma2)))
    
    KO.cali<-function(B,M,thinning,jump){
      ko.mc = matrix(NA, nrow = B+thinning*M, ncol=q)
      ko.mc[1,] = runif(q,lower.th,upper.th)
      sigma2  = vector(); sigma2[1] = tau2
      log.pi.ko<-function(th,sigma2) -
        0.5*t(y.p-ys.D(th))%*%solve(S.ko(sigma2))%*%(y.p-ys.D(th))
      for(iter in 2:(B+thinning*M)){
        for(j in 1:q){
          theta_new = runif(1, min=max(ko.mc[iter-1,j]-jump[j],lower.th[j]),
                            max=min(ko.mc[iter-1,j]+jump[j],upper.th[j]))
          theta.new = ko.mc[iter-1, ]
          theta.new[j] = theta_new
          loglik_old = log.pi.ko(ko.mc[iter-1,],sigma2[iter-1])
          loglik_new = log.pi.ko(theta.new,sigma2[iter-1])
          U = runif(1, min = 0, max = 1)
          if(log(U) < (loglik_new - loglik_old)){
            ko.mc[iter,j] = theta_new
          }else{
            ko.mc[iter,j] = ko.mc[iter-1,j]
          }
        }
        z = A.ko(ko.mc[iter,],sigma2[iter-1])#mvrnorm(1, A.ko(ko.mc[iter],sigma2[iter-1]), B.ko(sigma2[iter-1]))  
        sigma2.star = 1 / rgamma(1,n/2,t(z)%*%solve(R.ko)%*%z/2)
        loglik_old.s = log.pi.ko(ko.mc[iter,],sigma2[iter-1])
        loglik_new.s = log.pi.ko(ko.mc[iter,],sigma2.star)
        U = runif(1, min = 0, max = 1)
        if(log(U) < loglik_new.s - loglik_old.s){
          sigma2[iter] = sigma2.star
        }else{
          sigma2[iter] = sigma2[iter-1]
        }
      }
      ko.mc = ko.mc[seq(B+1, B+thinning*M, by=thinning),]
      result=list(ko.mc=ko.mc,sigma2 = sigma2[seq(B+1, B+thinning*M, by=thinning)])
      return(result)
    }
    result.ko = KO.cali(B,M,thinning,jump)
    ko.mh = as.matrix(result.ko$ko.mc) 
    theta.ko = apply(ko.mh,2,mean)
    sd.ko = apply(ko.mh,2,sd)
    lwr.ko=theta.ko-sd.ko*qnorm(0.975);hwr.ko=theta.ko+sd.ko*qnorm(0.975) ##Confidence Interval
    con.ko = as.matrix(rbind(lwr.ko,hwr.ko))
    result.KO = list(theta.hat=theta.ko, sd.hat=sd.ko, post.samples=ko.mh, con.hat=con.ko)
    return(result.KO)
  }
  else if(method == "SGP"){######SGP Calibration######
    #y.s1 <- function(x0,th){
    #  ys 
    #  return(ys)
    #}
    model_sgasp=rcalibration(design=D,observations=y.p,p_theta=q,simul_type=1,
                             math_model=y.s,theta_range=matrix(c(lower.th,upper.th),q,2)
                             ,S=(B-1)*thinning+M,S_0=M,thinning=thinning, discrepancy_type='S-GaSP')
    sgp.mcmc = as.matrix(model_sgasp@post_sample[,1:q])
    theta.sgp = apply(sgp.mcmc,2,mean)
    sd.sgp = apply(sgp.mcmc,2,sd)
    lwr.sgp=max(theta.sgp-sd.sgp*qnorm(0.975),lower.th);hwr.sgp=min(theta.sgp+sd.sgp*qnorm(0.975),upper.th)
    con.sgp = rbind(lwr.sgp,hwr.sgp)
    result.sgp = list(theta.hat=theta.sgp, sd.hat=sd.sgp, post.samples=sgp.mcmc, con.hat=con.sgp)
    return(result.sgp)
  }
  else if(method == "OGP"){######OGP Calibration######
    N=2000; x.test = as.matrix(maximinLHS(N,2)*(upper.th-lower.th)+lower.x)
    WW = R(x.test,x.test)
    R.co = R(D,D); w.ogp=R(D,x.test)
    der.ys<-function(x0,th0){
      d.new<-function(th) y.s(x0,th) 
      d.ys = jacobian(func=d.new, x=th0)  
      return(d.ys)
    }
    R.ogp<-function(th){
      F1 = t(as.matrix(apply(x.test,1,function(x) der.ys(x,th))))
      R.th = R.co-w.ogp%*%F1%*%solve(t(F1)%*%WW%*%F1)%*%t(F1)%*%t(w.ogp)
      return(R.th)
    }
    
    S.ogp <-function(th,sigma2) sigma2*R.ogp(th)+tau2*diag(n)
    A.ogp<-function(th,sigma2) sigma2*R.ogp(th)%*%solve(S.ogp(th,sigma2))%*%(y.p-ys.D(th))
    b.ogp<-function(th,sigma2) sigma2*R.ogp(th)-sigma2^2*R.ogp(th)%*%solve(S.ogp(th,sigma2))%*%R.ogp(th)
    B.ogp<-function(th,sigma2) 0.5*(b.ogp(th,sigma2)+t(b.ogp(th,sigma2)))
    
    OGP.cali<-function(B,M,thinning,jump){
      ogp.mc = matrix(NA, nrow = B+thinning*M, ncol=q);ogp.mc[1,] = runif(q,lower.th,upper.th)
      sigma2  = vector(); sigma2[1] = tau2
      log.pi.ogp<-function(th,sigma2) -0.5*log(det(S.ogp(th,sigma2)))-0.5*
        t(y.p-ys.D(th))%*%solve(S.ogp(th,sigma2))%*%(y.p-ys.D(th))
      for (iter in 2:(B+thinning*M)){
        theta.new = ogp.mc[iter-1,]
        for(j in 1:q){
          theta_new = runif(1, min=max(ogp.mc[iter-1,j]-jump[j],lower.th[j]),
                            max=min(ogp.mc[iter-1,j]+jump[j],upper.th[j]))
          theta.new[j] = theta_new
          loglik_old = log.pi.ogp(ogp.mc[iter-1,],sigma2[iter-1])
          loglik_new = log.pi.ogp(theta.new,sigma2[iter-1])
          U = runif(1, min = 0, max = 1)
          if(log(U) < (loglik_new - loglik_old)){
            ogp.mc[iter,j] = theta_new
          }else{
            ogp.mc[iter,j] = ogp.mc[iter-1,j]
          }
        }
        z = A.ogp(ogp.mc[iter,],sigma2[iter-1])#mvrnorm(1, A.ogp(ogp.mc[iter],sigma2[iter-1]), B.ogp(ogp.mc[iter],sigma2[iter-1]))  
        sigma2.star = 1 / rgamma(1,n/2,t(z)%*%solve(R.ogp(ogp.mc[iter,]))%*%z/2)
        loglik_old.s = log.pi.ogp(ogp.mc[iter,],sigma2[iter-1])
        loglik_new.s = log.pi.ogp(ogp.mc[iter,],sigma2.star)
        U = runif(1, min = 0, max = 1)
        if(log(U) < loglik_new.s - loglik_old.s){
          sigma2[iter] = sigma2.star
        }else{
          sigma2[iter] = sigma2[iter - 1]
        }  
      }
      ogp.mc = ogp.mc[seq(B+1, B+thinning*M, by=thinning),]
      result=list(ogp.mc=ogp.mc,sigma2 = sigma2[seq(B+1, B+thinning*M, by=thinning)])
      return(result)
    }
    result.ogp = OGP.cali(B,M,thinning,jump)
    ogp.mh = as.matrix(result.ogp$ogp.mc)
    theta.ogp = apply(ogp.mh,2,mean)
    sd.ogp = apply(ogp.mh,2,sd)
    lwr.ogp=max(theta.ogp-sd.ogp*qnorm(0.975),lower.th);hwr.ogp=min(theta.ogp+sd.ogp*qnorm(0.975),upper.th)
    con.ogp = as.matrix(rbind(lwr.ogp,hwr.ogp)) 
    result.OGP = list(theta.hat=theta.ogp, sd.hat=sd.ogp, post.samples=ogp.mh, con.hat=con.ogp)
    return(result.OGP)
  }
}

