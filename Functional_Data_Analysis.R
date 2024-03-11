#################################################################
#################################################################
# R examples I

library(fda)
library(refund)


library(ggplot2)
library(dplyr)
library(reshape2)

# basis functions
spline.basis<-create.bspline.basis(rangeval=c(0,10), nbasis=7)
par(mfrow=c(1,1))
plot(spline.basis, lty=1, lwd=2)


fourier.basis<-create.fourier.basis(rangeval=c(0,10), nbasis=7)
plot(fourier.basis, lty=1, lwd=2)


# Example based on Brownian motion
Nsim<-50
Times<-10000
W.mat=matrix(0, ncol=Nsim, nrow=Times)
for(n in 1:Nsim){W.mat[, n]=cumsum(rnorm(Times))/100}

plot.new()
for (i in 1:Nsim) {
  plot(1:Times,W.mat[,i],type = 'l',
       col=sample(1:Times,size = 1),ylim = c(min(W.mat),max(W.mat)),xlab = '',ylab = '') 
  par(new=TRUE)
}

plot.new()

# smooth by B-spline
B25.basis=create.bspline.basis(rangeval=c(0,Times), nbasis=25)
W.fd=smooth.basis(y=W.mat, fdParobj=B25.basis)

plot(W.fd$fd, ylab="", xlab="",col='gray',lty=1)

W.mean=mean.fd(W.fd$fd)
lines(W.mean, lty=2, lwd=3)


#NB! Fourier system is usually only suitable for functions which have
# approximately the same values at the beginning and the end of the interval (page 5 B1) 

# smooth by Fourier 
fourier.basis25<-create.fourier.basis(rangeval=c(0,Times), nbasis=25)
Wf.fd=smooth.basis(y=W.mat, fdParobj=fourier.basis25)

plot(Wf.fd$fd, ylab="", xlab="",col="gray",lty=1)


Wf.mean=mean.fd(Wf.fd$fd)
lines(Wf.mean, lty=2, lwd=3)

# smooth by B-spline
B25.basis=create.bspline.basis(rangeval=c(0,Times), nbasis=25)
W.fd=smooth.basis(y=W.mat, fdParobj=B25.basis)

plot(W.fd$fd, ylab="", xlab="",col="gray",lty=1)

W.mean=mean.fd(W.fd$fd)
lines(W.mean, lty=2, lwd=3)


# covariance function
W.cov=var.fd(W.fd$fd) 
grid=(1:100)*100
W.cov.mat=eval.bifd(grid, grid, W.cov)
persp(grid, grid, W.cov.mat, xlab="s",
      ylab="t", zlab="c(s,t)")

# PCA in L2
W.pca = pca.fd(W.fd$fd, nharm=4)
W.pca$values # eigenvalues for 25 basis
plot(W.pca$harmonics, lwd=3)
W.pca$varprop


# First four principal functions describe 95% of the data
sum(W.pca$varprop)

W.pca = pca.fd(W.fd$fd, nharm=10)

dataFrame1<-as.data.frame(cbind(1:10,round(W.pca$varprop,3)))
colnames(dataFrame1)<-c("PC","Varprop")
ggplot(data=dataFrame1, aes(x=PC, y=Varprop)) +
  ylim(0, 1)+
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=Varprop), vjust=-0.3, size=3.5)+
  theme_minimal()


library(refund)
data(DTI)
Y<-DTI$cca
Y<-Y[-c(126,130,131,125,319,321),] # missing values
N<-dim(Y)[1]; M<-dim(Y)[2]
argvals<-seq(0,1,length=M)
data_basis<-create.bspline.basis(c(0,1),nbasis=10)
data_basis<-create.fourier.basis(c(0,1), nbasis=10)
Y.f<-Data2fd(argvals,t(Y),data_basis)
plot(Y.f,lty=1,col="gray",xlab="",ylab="",ylim=c(0.1,.9))
lines(mean.fd(Y.f),lwd=2)
lines(mean.fd(Y.f)+std.fd(Y.f),lwd=2,lty=2,col="green")
lines(mean.fd(Y.f)-std.fd(Y.f),lwd=2,lty=2,col="green")

# PCA in L2
W.pca = pca.fd(Y.f, nharm=4)
plot(W.pca$harmonics, lwd=3)
W.pca$varprop

dataFrame1<-as.data.frame(cbind(1:4,round(W.pca$varprop,3)))
colnames(dataFrame1)<-c("PC","Varprop")
ggplot(data=dataFrame1, aes(x=PC, y=Varprop)) +
  ylim(0, 1)+
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=Varprop), vjust=-0.3, size=3.5)+
  theme_minimal()

#################################################################
#################################################################
# R examples II

library(fdasrvf)

data("simu_data")
f1<-simu_data$f[,1]
f2<-0.7*simu_data$f[,20]
ylim<-c(0,max(c(f1,f2)))
time<-simu_data$time
plot(time,f1,type = 'l',col='blue',xlab = '',ylab = '',ylim = ylim)
lines(time,f2,type = 'l',col='red',xlab = '',ylab = '',ylim = ylim)

q1<-f_to_srvf(f1, time)
q2<-f_to_srvf(f2, time)

ylimTemp<-c(min(c(q1,q2)),max(c(q1,q2)))
plot(time,q1,type = 'l',col='blue',xlab = '',ylab = '',ylim = ylimTemp)
lines(time,q2,type = 'l',col='red',xlab = '',ylab = '',ylim = ylimTemp)

out = pair_align_functions(f1,f2,time)
plot(time,f1,type = 'l',col='blue',xlab = '',ylab = '',ylim = ylim)
lines(time,out$f2tilde,type = 'l',col='red',xlab = '',ylab = '',ylim = ylim)

out2 = pair_align_functions(f2,f1,time)
plot(time,f2,type = 'l',col='blue',xlab = '',ylab = '',ylim = ylim)
lines(time,out2$f2tilde,type = 'l',col='red',xlab = '',ylab = '',ylim = ylim)

#NB! the aligned functions are also registered for all t in [a,b]

# Comparing warping functions
# out$gam should be the inverse of out2$gam based on inverse symmetry 
# see slide 29
plot(1:length(out$gam),out$gam,type = 'l',col='blue',xlab = '',ylab = '')
lines(out2$gam,col='red',xlab = '',ylab = '')

#################################################################
#################################################################
# R examples III

library(refund)
library(ggplot2)
library(dplyr)
library(reshape2)


set.seed(9000)

n = 1000
grid = seq(0, 1, length = 101)

# regressor functions
X <- matrix(0, nrow=n, ncol=length(grid))
for(i2 in 1:n){
  X[i2,]=X[i2,]+rnorm(length(grid), 0, 1)
  X[i2,]=X[i2,]+runif(1, 0, 5)
  X[i2,]=X[i2,]+rnorm(1, 1, 0.2)*grid
  for(j2 in 1:10){
    e=rnorm(2, 0, 1/j2^(2))
    X[i2,]=X[i2,]+e[1]*sin((2*pi)*grid*j2)
    X[i2,]=X[i2,]+e[2]*cos((2*pi)*grid*j2)
  }
}

#plot one regressor curve
plot(1:NCOL(X),X[10,],type = 'l')

#Beta function
beta = -dnorm(grid, mean=.2, sd=.03) +3*dnorm(grid, mean=.5,
                                              sd=.04)+dnorm(grid, mean=.75, sd=.05)
# response variables based on integral in discrete format
Y = X %*% beta * .01 + rnorm(n, 0, .4)
plot(1:dim(Y)[1],Y,type = 'l')

# "ps" stands for "penalized splines", fx= TRUE means no penalty is used
fit.lin = pfr(Y ~ lf(X, bs = "ps", k = 15, fx = TRUE))
# if sp is not specified, data driven smoothing is used
fit.pfr = pfr(Y ~ lf(X, bs = "ps", k =50))
# FPC
fit.fpcr = pfr(Y ~ fpc(X))

#NB coefficients are discrete functions as 
# (beta1,grid[1]),...,(beta100,grid[100])

# plot
coefs = data.frame(grid = grid,
                   FPCR = coef(fit.fpcr)$value,
                   Basis = coef(fit.lin)$value,
                   Penalized = coef(fit.pfr)$value,
                   Truth = beta)
coefs.m = melt(coefs, id = "grid")
colnames(coefs.m) = c("grid", "Method", "Value")
ggplot(coefs.m, aes(x = grid, y = Value, color = Method, group
                    = Method),width=12,height=6) + geom_path() + theme_bw()

par(mfrow=c(1,1))
plot(1:1000,Y,type = 'l',col='black')
lines(1:1000,X %*% coef(fit.fpcr)$value *0.01 ,type = 'l',col='blue')
plot(1:1000,Y,type = 'l',col='black')
lines(1:1000,X %*% coef(fit.lin)$value *0.01,type = 'l',col='red')
plot(1:1000,Y,type = 'l',col='black')
lines(1:1000,X %*% coef(fit.pfr)$value *0.01,type = 'l',col='green')
par(mfrow=c(1,1))

#################################################################
#################################################################
# R examples IV

require(fda)
require(refund)

daybasis25 <- create.fourier.basis(rangeval=c(0,365),
                                   nbasis=25,axes=list('axesIntervals'))
Temp.fd <- with(CanadianWeather, smooth.basisPar(day.5,
                                                 dailyAv[,,'Temperature.C'], daybasis25)$fd)
par(mfrow=c(1,1))
plot(Temp.fd)

modmat<-cbind(1,model.matrix(~factor(CanadianWeather$region)-1))
constraints <- matrix(c(0,1,1,1,1),1)

#function-on-scalar

# Ordinary Least Square estimator (OLS)
olsmod <- fosr(fdobj = Temp.fd, X = modmat, con = constraints,
              method="OLS", lambda=100*10:30)
# Generalized Least Square (GLS)
glsmod <- fosr(fdobj = Temp.fd, X = modmat, con = constraints,
              method="GLS")

par(mfrow=c(2,5), mar=c(5,2,4,1))
plot(olsmod, split=1, set.mfrow=FALSE,titles=c("OLS: Intercept (mean)", levels(factor(CanadianWeather$region))),ylab="", xlab="Day")
plot(glsmod, split=1, set.mfrow=FALSE,titles=c("GLS: Intercept (mean)", levels(factor(CanadianWeather$region))),ylab="", xlab="Day")
par(mfrow=c(1,1))
# GLS is slightly smoother than OLS

#################################################################
#################################################################
# R examples V

library(refund)
fplot <- function(x, y, lty=1, col=rgb(0, 0, 0, max(.1,sqrt(1/
                                                              nrow(y)))),
                  lwd=.5, ...){
  if (missing(x)) {
    if (missing(y))
      stop("must specify at least one of 'x' and 'y'")
    else x <- 1L:NROW(y)
  } else {
    if (missing(y)) {
      y <- t(x)
      x <- 1L:NROW(y)
    } else {
      y <- t(y)
    }
  }
  matplot(x, y, type="l", lty=lty, col=col, lwd=lwd, ...)
}

#generate 100 simulated curve for model 5.10
set.seed(9312)
data_ff <- pffrSim(scenario=c("int", "ff"), n=100)

#actual kernel 
psi_st <- function(s, t) {s*cos(pi * abs(s - t)) - .19}
s <- seq(0, 1, length = 40)
t <- seq(0, 1, length = 60)
psi_true <- outer(s, t, psi_st)

par(mfrow=c(1,2))
#all regressor curves
fplot(s, data_ff$X1, xlab = "s", ylab = "", main="X(s)")
# highlight the first four regressor functions X_i
matlines(s, t(data_ff$X1[1:4,]), col=1, lwd=2)

#all transformed regressor curves
fplot(t, attr(data_ff, "truth")$etaTerms$X1, xlab = "t",
      main=expression(integral(psi(t,s)*X(s)* ds)))
# highlight the first four transformed curves
matlines(t, t(attr(data_ff, "truth")$etaTerms$X1[1:4,]), col
         =1, lwd=2)

# find the kernel function
m_ff <- pffr(Y ~ ff(X1), data = data_ff)
psi_plot <- plot(m_ff, select = 2, pers=TRUE)[[2]]
layout(t(1:2))
# true psi surface
par(mfrow=c(1,2))
persp(s, t, psi_true, xlab = "s", ylab = "t", main=expression(
  psi(t,s)), phi=40, theta = 30, ticktype="detailed", zlab =
    "", border= NA, col="grey", shade = .7, zlim = c(-1,1))
# estimated psi surface

persp(psi_plot$x, psi_plot$y, matrix(psi_plot$fit, 40, 40),
      xlab = "s", ylab = "t", phi=40, theta = 30, ticktype="detailed",
      main=expression(hat(psi)(t,s)), zlim = c(-1,1),
      zlab = "", border= NA, col="grey", shade = .7)

par(mfrow=c(1,3))
fplot(t, data_ff$Y, xlab = "t", ylab="", main="Observations",
      ylim=range(data_ff$Y))
# matlines(t, t(data_ff$Y[1:4,]), col =1, lty = 1:4, lwd=2)
fplot(t, attr(data_ff, "truth")$eta, xlab = "t", ylab="", main
      ="True predictor", ylim=range(data_ff$Y))
# matlines(t, t(attr(data_ff, "truth")$eta[1:4,]), col = 1, lty= 1:4, lwd=2)
fplot(t, fitted(m_ff), xlab = "t", ylab="", main="Estimates",
      ylim=range(data_ff$Y))
# matlines(t, t(fitted(m_ff)[1:4,]), col = 1, lty = 1:4, lwd=2)
par(mfrow=c(1,1))


#################################################################
#################################################################
# R examples VI 

# function-on-scalar

library(refund)
library(MASS)

N<-200
M<-50
time = seq(0,1,length=M)
# Mean Function
mu_f<-function(t){cos(pi*t + pi)}
# Coefficient Function
beta_f<-function(t){2*t}
# Matern Covariance Function for the error term epsilon(t)
C_f<-function(t,s){
  sig2<-1; rho<-0.5
  d<-abs(outer(t,s,"-"))
  tmp2<-sig2*(1+sqrt(3)*d/rho)*exp(-sqrt(3)*d/rho)
}

set.seed(2000)
Sigma<-C_f(time,time)
mu<-mu_f(time)
beta<-beta_f(time)
X<-rnorm(N,mean=0)
Xdata<- data.frame(X = X)
Z<-mvrnorm(N,mu,Sigma) + X%*%t(beta)
Y<-matrix(Z>0,nrow=N)
as.integer(Y)
dim(Y) #200 functions

#plot 2 samples
par(mfrow=c(1,2))
plot(time,Y[1,])
plot(time,Y[2,])

# function-on-scalar with probit model
pffr_fit<-pffr(Y~X ,family=binomial(link="probit"),yind=time, data=Xdata)

par(mfrow=c(1,2))
plot(pffr_fit,select=1,xlab="",ylab="Intercept",ylim=c(-1.25,1.5),cex.lab=1.25)
points(time,mu_f(time),typ="l",lty=4,lwd=4)
plot(pffr_fit,select=2,xlab="t",ylab="Slope",ylim=c(-0.25,2.75),cex.lab=1.25)
points(time,beta_f(time),typ="l",lty=4,lwd=4)

#################################################################
#################################################################
# R examples VII

# function-on-function
library(refund)
library(MASS)

N<-1000
M<-50
time = seq(0,1,length=M)
# Mean Function
mu_f<-function(t){cos(2*pi*t/2 + pi)}
# Coefficient Function
beta_f<-function(t,s){
  d<-abs(outer(t,s,"*"))
  return(4*d)}
# Matern Covariance Function
C_f<-function(t,s){
  sig2<-1; rho<-0.5
  d<-abs(outer(t,s,"-"))
  tmp2<-sig2*(1+sqrt(3)*d/rho)*exp(-sqrt(3)*d/rho)}

set.seed(2000)
Sigma<-C_f(time,time)
mu<-mu_f(time)
beta<-beta_f(time,time)
X=mvrnorm(N,mu,Sigma)
Xdata<- data.frame(X = X)
Z<-mvrnorm(N,mu,Sigma) + X%*%t(beta)/M
Y<-matrix(Z>0,nrow=N)

#Takes time!
pffr_fit<-pffr(Y~ff(X,basistype="te",xind=time),
               family=binomial(link="probit"),
               yind=time, data=Xdata)
par(mfrow=c(2,2))
plot(pffr_fit,select=1,xlab="",ylab='',main="Estimated intercept")
plot(time,mu_f(time),typ="l",xlab="",ylab='',main="True intercept")
plot(pffr_fit,select=2,pers=TRUE,xlab="",ylab="",main="Estimated Slope")
persp(time,time,beta_f(time,time),xlab="",ylab="",zlab="True Slope",theta=30,phi=30)
par(mfrow=c(1,1))


#################################################################
#################################################################
# R examples VIII

# library(ggplot2)
library(fdasrvf) #(elastic metric) contains 2D images

#2D images of SRVF
data("mpeg7")
# shapeGroup numbers 14->camel, 23->cup, 24->dear, 25->star,
# 35->dog, 36->elephant, 47-> c-shape
shapeGroup<-2
tempMatrix<-t(beta[,,shapeGroup,1])
mesh2D<-tempMatrix[nrow(tempMatrix):1,] #reverse the matrix
plotshapes(rbind(mesh2D,mesh2D[1,]),joinline = TRUE)

q = curve_to_q(beta[,,shapeGroup,1])$q

plotshapes(rbind(t(q),t(q)[1,]),joinline = TRUE)

f0<-q_to_curve(q)

#Note that the scale is removed
plotshapes(rbind(t(f0),t(f0)[1,]),joinline = TRUE)


#################################################################
#################################################################
# R examples IX

library(RiemBase)
library(matlib)
library(shapes)
library(pracma)
library(rgl)
library(mvShapiroTest)
library(mvtnorm)

# converted code of Sungkyu Jung MATLAB randS2.m
# generate random sample of small sphere on S2 of the second kind 
# mu0, mu1 are directions and kappa0>1, kappa1>1
randS2 <- function(mu0,mu1,kappa0,kappa1,n) {
  
  mu0<-mu0/norm(mu0,type = "2")
  mu1<-mu1/norm(mu1,type = "2")
  nu<-sum(mu0*mu1)
  
  #generate Bingham-Mardia random vectors by the north pole
  
  x<-rnorm(n = n,mean = nu,sd = 1/sqrt(2*kappa0))
  x<-x[x<1 & x>-1]
  nnow<-length(x)
  while(n>nnow){
    n_more<-ceiling(n/nnow*(n-nnow))
    xmore<-rnorm(n = n_more,mean = nu,sd = 1/sqrt(2*kappa0))
    xmore<-xmore[xmore<1 & xmore>-1]
    x<-c(x,xmore)
    nnow<-length(x)
  }
  z<-x[1:n]
  
  #generate von Mises for longitude that c=mu0-nu*mu1 is parallel to x-axis
  theta<-randVonMises(mean = 0, kappa = kappa1, n = n)
  X_axis_northpole<-cbind(sqrt(1-z^2)*cos(theta),sqrt(1-z^2)*sin(theta), z)
  
  cx<-(mu1-nu*mu0)/sqrt(1-nu^2)
  cy<-cross(mu0,cx)
  cz<-mu0
  
  #rotate
  X<-X_axis_northpole%*%rbind(cx,cy,cz)
  
  return(X)
}

# draw circle on unit sphere S2 by center of small circle and r
# converted code of Sungkyu Jung MATLAB drawCircleS2.m
drawCircleS2 <- function(center,theta) {
  # NB!!! theta is the angle from center
  if(theta==pi/2){
    t<-cbind(cos(seq(0,2*pi,length.out = 50)),sin(seq(0,2*pi,length.out = 50)),rep(0,50))
    sCirc<-t%*%rotMat(center,c(0,0,1))
  }else{
    t<-cbind(sin(theta)*cos(seq(0,2*pi,length.out = 50)),sin(theta)*sin(seq(0,2*pi,length.out = 50)),cos(theta)*rep(1,50))
    sCirc<-t%*%rotMat(center,c(0,0,1))
  }
  spheres3d(x = 0, y = 0, z = 0, radius = 1,col = "lightblue", alpha=0.1)
  plot3d(sCirc,type="l",col = "black",expand = 10,box=TRUE,add = TRUE)
}


# copyright belongs to Sungkyu Jung
# converted code from MATLAB
# output is the type of fitted circle 
# likelihood ratio test from shapes package
LRTpval2 <- function(resGreat,resSmall,n) {
  chi2 <- max(n*log(sum(resGreat^2)/sum(resSmall^2)))
  pval <- 1-pchisq(q = chi2, df = 1, lower.tail = T) # likelihood test p-value Also you can use chi2cdf(chi2,1) from library(PEIP) like matlab
}

kurtosisTestFunction <- function(sphericalData, alpha=0.1) {
  ndata<-dim(sphericalData)[2]
  
  subsphereSmall<-getSubSphere(sphericalData,geodesic = "small")
  subsphereGreat<-getSubSphere(sphericalData,geodesic = "great")
  
  currentSphere<-sphericalData
  
  rSmall<-subsphereSmall$r                     # rSmall is rs in matlab
  centerSmall<-subsphereSmall$center           # NB! center is the centerSmall is pnsSmall$PNS$orthaxis[[1]]
  # and centers in matlab
  resSmall <- acos(t(centerSmall)%*%currentSphere)-rSmall  # NB!!! resSmall==(pnsSmall$resmat)[2,] i.e., residuals are second coordinates of PNS
  
  rGreat<-subsphereGreat$r                     # rGreat is rg in matlab
  centerGreat<-subsphereGreat$center           # centerGreat is centers in matlab
  resGreat <- acos(t(centerGreat)%*%currentSphere)-rGreat  # NB!!! resGreat==(pnsGreat$resmat)[2,] i.e., residuals are second coordinates of PNS
  
  # LRTpval is the likelihood ratio test from 'shapes' package
  # Chi-squared statistic for a likelihood test
  pval1 <- LRTpval(resGreat,resSmall,n = ndata)
  pval1
  
  if(pval1>alpha){
    print('great by likelihood ratio test')
    return('great')
    break
  }
  
  # # equivalently we can find pval by pns function
  # pnsTest2<-pns(sphericalData)
  # pnsTest2$PNS$pvalues
  # sum(pnsTest2$resmat[2,]==resSmall)
  
  # kurtosis test routine
  X <- LogNPd(rotMat(centerSmall) %*% currentSphere)
  
  # plot3d(t(sphericalData),type="p",expand=10, add=TRUE)
  # plot3d(t(rbind(X,rep(1,dim(X)[2]))),type="p",col = "blue",expand=10, add=TRUE)
  
  # Note that the tangential point is the center of the small circle
  d<-dim(X)[1]
  n<-dim(X)[2]
  normX2 <- colSums(X^2)
  kurtosis <- sum( normX2^2 ) / n / ( sum( normX2 ) / (d * (n-1)) )^2
  M_kurt <- d * (d+2)^2 / (d+4)
  V_kurt <- (1/n) * (128*d*(d+2)^4) / ((d+4)^3*(d+6)*(d+8))
  pval2 <- pnorm((kurtosis - M_kurt) / sqrt(V_kurt))
  
  if(pval2>alpha){
    return('great')
  }else{
    # drawCircleS2(normalVec = centerSmall,radius = rSmall)
    return('small')
  }
}

# generate random von Mises distribution on circle in radian
# converted code of Sungkyu Jung 2013, and Byungwon Kim 2017, MATLAB randvonMises.m
# mean is in[0,2pi] and kappa>0
randVonMises <- function(mean, kappa, n) {
  tau<-1+sqrt(1+4*kappa^2)
  rho<-(tau-sqrt(2*tau))/(2*kappa)
  r<-(1+rho^2)/(2*rho)
  
  u1<-runif(n)
  z<-cos(pi*u1)
  f<-(1+r*z)/(r+z)
  c<-kappa*(r-f)
  u2<-runif(n)
  acceptanceid<-(c*(2-c)-u2>0) | (log(c/u2)+1-c>=0)
  u3<-runif(sum(acceptanceid))
  theta<-sign(u3-0.5)*acos(f[acceptanceid])
  nnow<-length(theta)
  
  while (n>nnow) {
    n_more<-ceiling(n/nnow*(n-nnow))
    u1<-runif(n_more)
    z<-cos(pi*u1)
    f<-(1+r*z)/(r+z)
    c<-kappa*(r-f)
    u2<-runif(n_more)
    acceptanceid<-(c*(2-c)-u2>0) | (log(c/u2)+1-c>=0)
    u3<-runif(sum(acceptanceid))
    thetamore<-sign(u3-0.5)*acos(f[acceptanceid])
    
    theta<-c(theta, thetamore)
    nnow<-length(theta)
  }
  
  theta<-theta[1:n] + mean
  
  return(theta)
}


randomData<-randS2(mu0 = c(0,0,1),mu1 = c(cos(pi/3),0,sin(pi/3)),kappa0 = 500,kappa1 = 2,n = 1000)
spheres3d(x = 0, y = 0, z = 0, radius = 1,col = "lightblue", alpha=0.1)
plot3d(randomData,type="p",col = "blue",expand = 10,box=TRUE,add = TRUE)

# PNS with small circle
smallCirclePns<-pns(t(randomData),sphere.type = "small")

# PNS with great circle
greatCirclePns<-pns(t(randomData),sphere.type = "great")

# Extrinsic mean
extrinsic_Mean<-colMeans(randomData)/sqrt(sum(colMeans(randomData)^2))

#plot
plot3d(randomData,type="p",col = "blue",expand = 10,box=TRUE,add = TRUE)
plot3d(rbind(smallCirclePns$PNS$mean,smallCirclePns$PNS$mean),type="s",radius = 0.015,col = "yellow",expand = 10,box=TRUE,add = TRUE)
plot3d(rbind(greatCirclePns$PNS$mean,greatCirclePns$PNS$mean),type="s",radius = 0.015,col = "orange",expand = 10,box=TRUE,add = TRUE)
drawCircleS2(center = greatCirclePns$PNS$orthaxis[[1]],theta = asin(greatCirclePns$PNS$radii[[2]]))
drawCircleS2(center = smallCirclePns$PNS$orthaxis[[1]],theta = asin(smallCirclePns$PNS$radii[[2]]))
plot3d(rbind(extrinsic_Mean,extrinsic_Mean),type="s",radius = 0.015,col = "red",expand = 10,box=TRUE,add = TRUE)



#plot residuals (Euclideanized data by PNS)
par(mfrow=c(1,2))
plotshapes(t(smallCirclePns$resmat))
title("PNS Residuals Small Circle")
plotshapes(t(greatCirclePns$resmat))
title("PNS Residuals Great Circle")
par(mfrow=c(1,1))


#################################################################
#################################################################
# R examples X

# Karcher mean for curves
library(fdasrvf)

data("simu_data")
plot(simu_data$time,simu_data$f[,1],type = 'l',ylim = c(0,1.5),xlab = '',ylab = '')
for (i in 1:ncol(simu_data$f)) {
  lines(simu_data$time,simu_data$f[,i])
}

#Group alignment 
# Takes time!
out <- function_group_warp_bayes(simu_data$f, simu_data$time, showplot = TRUE)

plot(simu_data$time,out$f_a[,1],type = 'l',ylim = c(0,1.5),xlab = '',ylab = '')
for (i in 1:ncol(simu_data$f)) {
  lines(simu_data$time,out$f_a[,i])
}

#################################################################
#################################################################
# R examples XI

# Karcher mean for planar closed curves
data("mpeg7")

plot(t(beta[,,1,1]),type = 'l',xlab = '',ylab = '')
plot(t(beta[,,1,2]),type = 'l',xlab = '',ylab = '')
plot(t(beta[,,1,3]),type = 'l',xlab = '',ylab = '')
plot(t(beta[,,1,4]),type = 'l',xlab = '',ylab = '')

data("mpeg7")
out = curve_srvf_align(beta[,,1,1:20],maxit=2)

shapes<-array(NA,dim = c(100,2,20))
for (i in 1:dim(out$qn)[3]) {
  shapes[,,i]<-t(q_to_curve(out$qn[,,i], scale = 1))
}

library(shapes)
plotshapes(shapes,joinline = TRUE)

meanShape<-t(q_to_curve(out$q_mu, scale = 1))
plotshapes(meanShape,joinline = TRUE)


#################################################################
#################################################################
# R examples XII

library(fda)
library(ggplot2)

m=100 # each function is observed at m+1 points, including 0 and 1
burnin=50 # the first 50 functions are a burn in period
N=200 # number of functions to simulate
N1=N+burnin
alpha=9/4
# Create 2 matrices, whose entries are 0s.
# Each column represents one function.
X<- matrix(rep(0, (m+1)*N1),m+1,N1)

epsilon<- matrix(rep(0, (m+1)*N1),m+1,N1)
epsilon[,1]<-rnorm(1)*sin(pi*(0:m/m))+0.5*rnorm(1)*cos(2*pi*(0:m/m))

# the following loop simulates FAR(1).
for(i in 2:N1){
  epsilon[,i]<-rnorm(1)*sin(pi*(0:m/m))+0.5*rnorm(1)*cos(2*pi*(0:m/m))
  X[,i]<-alpha*(1/m)^2*sum((1:m)*X[2:(m+1),i-1])*(0:m/m)+
    epsilon[,i]
}
X=X[,-(1:burnin)] # Remove the burn in period functions

last=10
plot.ts(c(X[,(N-last+1):N]),ylim=c(min(X[,(N-last):N])-0.5,0.5
                                   +max(X[,(N-last):N])),axes=F,xlab="",ylab="",lwd=2)
axis(2)
axis(1,tick=F,labels=F)
abline(h=0)
abline(v=seq(0,last*(m+1),by=m+1), lty=2)
box()

basisfd=10 # number of basis functions to represent each functional observation
basis=create.bspline.basis(c(0,1),nbasis=basisfd,norder=4)
fdX=Data2fd(argvals=0:m/m,X, basis)

p=4 # number of EFPC's
fdXpca=pca.fd(fdX, nharm=p)
eigenvalues=fdXpca$values
scoresX=fdXpca$scores
# jth column of scoresX contains scores of the jth EFPC
harmonicsX=fdXpca$harmonics # extract the EFPC's

varianceprop=fdXpca$varprop #proportion of variance explained by the EFP's

# we truncate the basis by j=4
round(varianceprop*100,0)
# plot 
dataFrame1<-as.data.frame(cbind(1:length(varianceprop),round(varianceprop,3)))
colnames(dataFrame1)<-c("PC","Varprop")
ggplot(data=dataFrame1, aes(x=PC, y=Varprop)) +
  ylim(0, 1)+
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=Varprop), vjust=-0.3, size=3.5)+
  theme_minimal()

phi=9/4*(0:m/m)%*%t(0:m/m)
# True surface evaluated on a discrete bivariate grid
par(mar=c(1.5, 1.5, 3.5, 0.2), mfrow=c(2,2), oma = c(4, 4, 0.2, 0.2))

# 4 panels 2 rows and 2 columns arrangement
axislabelsize=1.5 # controls the size of labels of axes
axisticksize=0.8 # controls the size of ticks of axes
persp((0:m/m),(0:m/m),z=phi,cex.axis = axisticksize,
      cex.lab=axislabelsize, xlab="t",ylab="s", zlab=" ", theta=30,phi=30,
      ticktype="detailed", main="True")

# Next we compute hat(phi)
# vivj is the matrix whose entries are products of v_j(s_k)*v_i(t_l).
# Blocks of vivj of m by m matrices represent products of
# v_j(s)v_i(t) evaluated on the (m+1) by (m+1) grid
for(npc in 1:3){
  vivj=matrix(0,p*(m+1),p*(m+1))
  for(j in 1:npc){
    for(i in 1:npc){
      vivj[1:(m+1)+(m+1)*(i-1),1:(m+1)+(m+1)*(j-1)]=
        eval.fd(evalarg=0:m/m, harmonicsX[i])%*%t(eval.fd(evalarg=0:m/m, harmonicsX[j]))
    }
  }
  # phip will be the estimated surface.
  phip=matrix(0,m+1,m+1)
  for(k in 1:(N-1)){
    temp=matrix(0,m+1,m+1)
    for(j in 1:npc){
      temp1=matrix(0,m+1,m+1)
      for(i in 1:npc){
        temp1=temp1+scoresX[k+1,i]*vivj[1:(m+1)+(m+1)*(i-1),1:(m+1)+(m+1)*(j-1)]
      }
      temp=temp+(eigenvalues[j])^(-1)*scoresX[k,j]*temp1
    }
    phip=phip+temp
  }
  phip=(1/(N-1))*phip
  if (npc==1)
    persp((0:m/m),(0:m/m), z=phip, cex.axis =axisticksize,cex.lab=axislabelsize,xlab="t",ylab="s", zlab=" ", theta=30, phi=30,
          ticktype="detailed", main="p=1")
  else if (npc==2)
    persp((0:m/m),(0:m/m), z=phip, cex.axis = axisticksize,cex.lab=axislabelsize,
          xlab="t", ylab="s",zlab=" ",theta=30,phi=30,ticktype="detailed",main="p=2")
  else if (npc==3)
    persp((0:m/m),(0:m/m), z=phip, cex.axis = axisticksize,cex.lab=axislabelsize,
          xlab="t", ylab="s", zlab=" ", theta=30, phi=30, ticktype="detailed",main="p=3")
}

