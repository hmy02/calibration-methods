library(cubature) ##Computing integral
library(lava) ##Computing trace
library(DiceKriging) ##Fitting Kringing Model
library(nloptr) ##Optimization
library(geoR)
library(lhs)
library(numDeriv)
library(RobustCalibration)
library(RobustGaSP)
library(MASS)
library(lhs)

# eg
zeta<-function(x) 2*x[1]*exp(-2*x[2])
y.s<-function(x,th) th*x[1]*exp(-th*x[2])
lx = c(0,0); ux = c(1,1)
lth = 1.9; uth = 2.1
theta.GALS=rep(0,30)
theta.ALS=rep(0,30)
theta.GOLS=rep(0,30)
tau2=0.05^2
ne=20
ns=20

D.p= maximinLHS(ne,2)*(ux-lx)+lx #real design



for(k in 1:30){
  D = maximinLHS(ns,3) #simulation samples
  D[,2:3] = as.matrix(D[,2:3]*(ux-lx)+lx) #simulation design
  D[,1]= as.matrix(D[,1]*(uth-lth)+lth) #simulation calibration parameters
  
  ys3=apply(D,1,function(x0) y.s(x0[2:3],x0[1])) #computer model outputs
  model.yshat = rgasp(design=D,response=ys3,kernel_type='pow_exp',alpha=rep(2,dim(D)[2])) #simulator
  beta=model.yshat@beta_hat
  sigmaZ2=model.yshat@sigma2_hat
  plot(predict(model.yshat,testing_input=D)$mean,ys3) #rgasp fits ys well
  plot(1:ne,apply(D.p,1,zeta))
  ys.hat<-function(x,th){ #simulator model
    test.new=t(as.matrix(c(th,x)))
    pre.test=predict(model.yshat,testing_input=test.new)$mean
    pre.sd=predict(model.yshat,testing_input=test.new)$sd
    return(list(mean=pre.test,sd=pre.sd))
  }
  ys.D<-function(th) apply(D.p,1,function(x) y.s(x,th))
  y.p = apply(D.p,1,zeta)+rnorm(ne,0,sqrt(tau2)) #real experiments
  D_th<-function(th) cbind(matrix(rep(th,ne),ncol=1,byrow=T),D.p)
  Rij<-function(v1,v2){
    a<-matrix(nrow=nrow(v1),ncol=nrow(v2))
    for(i in 1:nrow(v1)){
      for(j in 1:nrow(v2)){
        a[i,j]<-exp(-(sum((v1[i,]-v2[j,])^2/beta)))
      }
    }
    return(a)
  }
  
  R<-function(th) Rij(D_th(th),D_th(th))
  U<-function(th) sigmaZ2*Rij(D_th(th),D)
  Vs<-sigmaZ2*Rij(D,D)
  Fe<-D.p
  Fs<-D[,2:3]
  O<-matrix(0,nrow=2,ncol=2)
  
  item1<-function(th) cbind(U(th),Fe)
  matr1 = cbind(Vs,Fs); matr2 = cbind(t(Fs),O); item2=rbind(matr1,matr2) 
  
  K<-function(th) sigmaZ2*(R(th)-as.matrix(item1(th))%*%solve(item2)%*%t(as.matrix(item1(th))))+tau2*diag(ne)
  
  GALS<-function(th){
    item11 = as.matrix(y.p-apply(D.p,1,function(x) ys.hat(x,th)$mean))
    LS = t(item11)%*%solve(K(th))%*%item11
    return(LS)
  }
  
  theta.GALS[k] = directL(GALS, lth, uth,
                           nl.info = TRUE, control=list(xtol_rel=1e-8, maxeval=1000))$par

}


model.tv<-function(da,th) (zeta(da)-y.s(da,th))^2
tv.loss<-function(th) adaptIntegrate(model.tv,lx,ux,th=th)$integral
theta.star = directL(tv.loss, lth, uth,
                     nl.info = TRUE, control=list(xtol_rel=1e-4, maxeval=1000))$par

sqrt(mean((theta.GALS[1:20]-2)^2)+var(theta.GALS[1:20]))

mean((theta.GALS-2)^2/4)

mean(theta.GALS)

