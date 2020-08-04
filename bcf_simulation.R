# data generating process
#p = 2 #two control variables
library(tidyverse)
library("bcf")
library(latex2exp)

source('~/Documents/GitHub/bcf_discussion/functions.R')
p=2
n = 250
alpha = 1
n_burn <- 500
n_sim <- 2000

#x = matrix(runif(n*p, min = 0, max = 1), nrow=n)

mu_function<- function(x1, x2){
  return(-3 + 6*pnorm(2*(x2 - x1)))  
}

pi_function<- function(x1, x2, alpha = 1){
  return( alpha*(0.8*pnorm(mu_function(x1, x2)/ (0.1*(2-x1- x2) +0.25 )) + 0.025*(x1 +x2) ) +0.05*(0.5*(19 - 17*alpha)) ) 
}

n_s = 30
x1<-seq(0,1,by=(2/(n_s - 1)))
x2<-seq(0,1, by=(2/(n_s - 1)))
mu_<- outer(x1,x2,Vectorize(mu_function))
pi_1 <-  outer(x1,x2,Vectorize(pi_function),alpha=1)
#(q- pnorm(-3))*(6/(pnorm(3)- pnorm(-3)) 

n_s = 30
x1<-seq(0,1,by=(2/(n_s - 1)))
x2<-seq(0,1, by=(2/(n_s - 1)))
pi_05 <-  outer(x1,x2,Vectorize(pi_function),alpha=0.5)


n_s = 30
x1<-seq(0,1,by=(2/(n_s - 1)))
x2<-seq(0,1, by=(2/(n_s - 1)))
mu_<- outer(x1,x2,Vectorize(mu_function))
pi_0 <-  outer(x1,x2,Vectorize(pi_function),alpha=0)

# applying smoothing to z matrix of integers

#pi_0 = lissage(pi_0) # do it about ten times :-)


pdf(file="mu.pdf",width=10,height=6)
build3ds1(x1,x2,mu_,par1= expression(paste(italic(mu))) )
dev.off()
pdf(file="pi_1.pdf",width=10,height=6)
build3ds1(x1,x2,pi_1,z_lim=c(0,1),par1= expression(paste(italic(pi))) )
dev.off()
pdf(file="pi_05.pdf",width=10,height=6)
build3ds1(x1,x2,pi_05,z_lim=c(0,1),par1= expression(paste(italic(pi))) )
dev.off()

pdf(file="pi_0.pdf",width=10,height=6)
build3ds1(x1,x2,pi_0,z_lim=c(0,1), par1= expression(paste(italic(pi))) )
dev.off()

b=ggplot()
x_tilde = x[,1]+  x[,2]
mu= seq(-3,3,by =0.05)
for(i in 1:20){
  b <- b + geom_line(aes(x=mu,y=pi), 
                     data = tibble(mu = mu, pi = 0.8* pnorm(mu/ (0.1*(2-x_tilde[i]) +0.25) ) + 0.025*(x_tilde[i]) +0.05 ))
}
print(b)

pdf(file="pi_vs_mu.pdf",width=5,height=3)
print(b+xlab(TeX(sprintf('$\\mu$')))+ylab(TeX(sprintf('$\\pi$'))) +theme_bw())
dev.off()

b=ggplot() + scale_color_gradient2(midpoint=mid, low="blue", mid="white",
                                   high="red", space ="Lab" )
#x_tilde = x[,1]+  x[,2]
x_tilde =1
mu= seq(-3,3,by =0.05)
alpha_seq=seq(0,0.5,by =0.05)
for(i in 1:length(alpha_seq)){
  b <- b + geom_line(aes(x=mu,y=pi), 
                     data = tibble(mu = mu, pi = alpha_seq[i]*(0.8*pnorm(mu/ (0.1*(2-x_tilde) +0.25) ) +0.025*(x_tilde)) +0.05*(0.5*(19 - 17*alpha_seq[i])) ), color = i)
}
print(b)


simulation <- function(i,alpha ){
  set.seed(i)
  #### control variables matrix x_i ~ Unif(0,1)
  x = matrix(runif(n*p, min = 0, max = 1), nrow=n)
  # create targeted selection
  q= -3 + 6*pnorm((2*(x[,2]- x[,1]))) 
  # generate treatment variable
  # probability of recievng the treatment
  pi = alpha*(0.8* pnorm(q/ (0.1*(2-x[,1]- x[,2]) +0.25 )) + 0.025*(x[,1] +x[,2])) +0.05*(0.5*(19 - 17*alpha))
  #pi = alpha*(0.8* pnorm(q/ (0.05*(4-x[,1]- x[,2]) +0.25 )) + 0.0125*(x[,1]+ x[,2]) ) +0.05*(0.5*(19 - 17*alpha))
  #treatment
  z = rbinom(n,1,pi)
  # tau is the true (homogeneous) treatment effect
  tau = 1
  # generate the response using q, tau and z
  mu = (q + tau*z)
  # set the noise level relative to the expected mean function of Y
  #sigma = diff(range(q + tau*pi))/8
  sigma=1
  # draw the response variable with additive error
  y = mu + sigma*rnorm(n)
  
  # If you didn't know pi, you would estimate it here
  pihat = pnorm(q)
  
  bcf_fit = bcf(y, z, x,x, pihat, nburn=n_burn, nsim=n_sim,include_pi= "control")
  ## convergence assessment for BCF
  #summary(bcf_fit)
  # Get posterior of treatment effects
  tau_ests <- data.frame(Mean  = colMeans(bcf_fit$tau),
                         Low95 = apply(bcf_fit$tau, 2, function(x) quantile(x, 0.025)),
                         Up95  = apply(bcf_fit$tau, 2, function(x) quantile(x, 0.975)))
  tau_post = bcf_fit$tau
  tauhat = colMeans(tau_post)
  ## rmse error
  RMSE= sqrt(sum((1 - tauhat)^2)/n)
  ## bias
  bias= mean(tauhat - 1)
  ## coverage
  isCovered <- function(i){
    ifelse(tau_ests$Low95[i] <= tau & tau <= tau_ests$Up95[i], 1, 0)
  }
  coverage <- lapply(1:length(tau_ests$Mean), isCovered)
  perCoverage <- sum(unlist(coverage))/length(tau_ests)
  bcf_ate = mean(tauhat)
  #Bart
  bartc1 <- bartc(response = y, treatment = z, confounders = x, method.rsp = "bart", method.trt = "none", n.samples = 3000)
  bart_ate = summary(bartc1)[9]$estimates[1]
  #ites <- extract(bartc1, type = "ite")
  #tau_x  = apply(ites, 2, mean)
  #RMSE_bart= sqrt(sum((1 - tau_x)^2)/n)
  #RMSE_bart
  return(list(bcf_ate =bcf_ate, bcf_cover =perCoverage, bcf_bias = bias, bcf_rmse= RMSE, bart_ate = bart_ate))
}

bart_ate_<-c()
bcf_ate_<- c()

for (i in 1:10){
  sim= simulation(i, alpha =0.5)
  bart_ate_ <- c(bart_ate_, sim$bart_ate)
  bcf_ate_ <- c(bcf_ate_, sim$bcf_ate)
}

#sqrt(mean((unlist(bart_ate_)  - 1)^2))

fin_rmse_bart = sqrt(mean((unlist(bart_ate_)  - 1)^2))
fin_rmse_bcf = sqrt(mean((unlist(bcf_ate_)  - 1)^2))

fin_rmse_bart
fin_rmse_bcf

bias_bart= mean((unlist(bart_ate_)  - 1))
bias_bcf =  mean((unlist(bcf_ate_)  - 1))
bias_bart
bias_bcf
# 
isCovered <- function(i){
  ifelse(tau_ests$Low95[i] <= tau & tau <= tau_ests$Up95[i], 1, 0)
}
coverage <- lapply(1:length(tau_ests$Mean), isCovered)
perCoverage <- sum(unlist(coverage))/length(tau_ests)
bcf_ate = mean(tauhat)
coverage
# 

bart_Low95 = quantile(unlist(bart_ate_), 0.025)
bart_Up95  = quantile(unlist(bart_ate_), 0.975)

bcf_Low95 = quantile(unlist(bcf_ate_), 0.025)
bcf_Up95  = quantile(unlist(bcf_ate_), 0.975)

# library(bartCause)
# bartc1 <- bartc(response = y, treatment = z, confounders = x, method.rsp = "bart", method.trt = "none", n.samples = 3000)
# summary(bartc1)
# 
# ites <- extract(bartc1, type = "ite")
# tau_x  = apply(ites, 2, mean)
# RMSE_bart= sqrt(sum((1 - tau_x)^2)/n)
# RMSE_bart
# 
# 
# 
# ####################################################################################################
# #bart_ <- bart(x, y, )
# 


