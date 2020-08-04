# data generating process
#p = 2 #two control variables
library(tidyverse)
library("bcf")
library(latex2exp)

source('~/Documents/GitHub/bcf_discussion/functions.R')
p=2
n = 250
alpha = 1
n_burn <- 1000
n_sim <- 3000

#
set.seed(1)
#### control variables matrix x_i ~ Unif(0,1)
x = matrix(runif(n*p, min = 0, max = 2), nrow=n)
# create targeted selection
q= -3 + 6*pnorm((x[,2]- x[,1])) 
#sq<- summary(q)
#sq

#l = (x[,2]- x[,1])
#plot(density(l))


mu_function<- function(x1, x2){
  return(-3 + 6*pnorm((x2 - x1)))  
}

pi_function<- function(x1, x2, alpha = 1){
    return( alpha*(0.8*pnorm(mu_function(x1, x2)/ (0.05*(4-x1- x2) +0.25 )) + 0.0125*(x1 +x2) ) +0.05*(0.5*(19 - 17*alpha)) ) 
  #return( alpha*(0.8*pnorm(mu_function(x1, x2)/ (0.1*(2-x1- x2) +0.25 )) + 0.025*(x1 +x2)) + 0.05*(0.5*(19 - 17*alpha)) ) 
}

n_s = 30
x1<-seq(0,2,by=(2/(n_s - 1)))
x2<-seq(0,2, by=(2/(n_s - 1)))
mu_<- outer(x1,x2,Vectorize(mu_function))
pi_ <-  outer(x1,x2,Vectorize(pi_function))
#(q- pnorm(-3))*(6/(pnorm(3)- pnorm(-3)) 




# applying smoothing to z matrix of integers

#pi_ = lissage(pi_) # do it about ten times :-)


pdf(file="mu.pdf",width=10,height=6)
build3ds1(x1,x2,mu_,par1= expression(paste(italic(mu))) )
dev.off()
pdf(file="pi.pdf",width=10,height=6)
build3ds1(x1,x2,pi_,par1= expression(paste(italic(pi))) )
dev.off()



# generate treatment variable
# probability of recievng the treatment
#pi = alpha*(0.8* pnorm(q/ (0.1*(2-x[,1]- x[,2]) +0.25 )) + 0.025*(x[,1] +x[,2])) +0.05*(0.5*(19 - 17*alpha))
pi = alpha*(0.8* pnorm(q/ (0.05*(4-x[,1]- x[,2]) +0.25 )) + 0.0125*(x[,1]+ x[,2]) ) +0.05*(0.5*(19 - 17*alpha))

###

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
summary(bcf_fit)




# Get posterior of treatment effects


tau_ests <- data.frame(Mean  = colMeans(bcf_fit$tau),
                       Low95 = apply(bcf_fit$tau, 2, function(x) quantile(x, 0.025)),
                       Up95  = apply(bcf_fit$tau, 2, function(x) quantile(x, 0.975)))



tau_post = bcf_fit$tau
tauhat = colMeans(tau_post)
#plot(tauhat,rep(tau,250));
#abline(0,1)

## rmse error
RMSE= sqrt(sum((1 - tauhat)^2)/n)
RMSE

## bias
bias= mean(tauhat - 1)
bias

## coverage
isCovered <- function(i){
  ifelse(tau_ests$Low95[i] <= tau & tau <= tau_ests$Up95[i], 1, 0)
}

coverage <- lapply(1:length(tau_ests$Mean), isCovered)
perCoverage <- sum(unlist(coverage))/length(tau_ests)
perCoverage

mean(tauhat)


#### BART 
library(bartCause)
bartc1 <- bartc(response = y, treatment = z, confounders = x, method.rsp = "bart", method.trt = "none", n.samples = 3000)
summary(bartc1)

ites <- extract(bartc1, type = "ite")
tau_x  = apply(ites, 2, mean)
RMSE_bart= sqrt(sum((1 - tau_x)^2)/n)
RMSE_bart



####################################################################################################
#bart_ <- bart(x, y, )



