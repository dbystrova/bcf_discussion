### Plots for the example 1

# data generating process
#p = 2 #two control variables
library(tidyverse)
library(bcf)
library(latex2exp)
library(viridis)
library(bartCause)
##change directory for functions.R file
source('~/Documents/GitHub/bcf_discussion/functions.R')
p=2
n = 250
alpha = 1
n_burn <- 500
n_sim <- 2000

### x_i ~ N(0,1)
x = matrix(runif(n*p, min = 0, max = 1), nrow=n)


## prognostic function \mu
mu_function<- function(x1, x2){
  return(-3 + 6*pnorm(2*(x1 - x2)))  
}


## propensity function \pi
pi_function<- function(x1, x2, alpha = 1){
  return( alpha*(0.8*pnorm(mu_function(x1, x2)/ (0.1*(2-x1- x2) +0.25 )) + 0.025*(x1 +x2) ) +0.05*(0.5*(19 - 17*alpha)) ) 
}

n_s = 70
x1<-seq(0,1,by=(2/(n_s - 1)))
x2<-seq(0,1, by=(2/(n_s - 1)))
mu_<- outer(x1,x2,Vectorize(mu_function))
pi_1 <-  outer(x1,x2,Vectorize(pi_function),alpha=1)
#(q- pnorm(-3))*(6/(pnorm(3)- pnorm(-3)) 
pi_05 <-  outer(x1,x2,Vectorize(pi_function),alpha=0.5)
pi_0 <-  outer(x1,x2,Vectorize(pi_function),alpha=0)

# applying smoothing to z matrix of integers
#pi_05 = lissage(pi_05) # do it about ten times :-)

## plot for mu(x_1, x_2) 
pdf(file="plots/mu.pdf",width=5,height=3)
build3ds1(x1,x2,mu_,par1= expression(paste(italic(mu))) )
dev.off()


## plot for pi(x_1, x_2) & alpha =1
pdf(file="plots/pi_1.pdf",width=10,height=9)
build3ds1(x1,x2,pi_1,z_lim=c(0,1),par1= expression(paste(italic(alpha)," = 1")) )
dev.off()

## plot for pi(x_1, x_2) & alpha =1/2
pdf(file="plots/pi_05.pdf",width=10,height=9)
build3ds1(x1,x2,pi_05,z_lim=c(0,1),par1= expression(paste(italic(alpha)," = 0.5")))
dev.off()

## plot for pi(x_1, x_2) & alpha =0
pdf(file="plots/pi_0.pdf",width=10,height=9)
build3ds1(x1,x2,pi_0,z_lim=c(0,1),par1= expression(paste(italic(alpha)," = 0")))
dev.off()


## plot for dependence between pi and mu (alpha =1)
b=ggplot()
x_tilde = x[,1]+  x[,2]
mu= seq(-3,3,by =0.05)
for(i in 1:20){
  b <- b + geom_line(aes(x=mu,y=pi), 
                     data = tibble(mu = mu, pi = 0.8* pnorm(mu/ (0.1*(2-x_tilde[i]) +0.25) ) + 0.025*(x_tilde[i]) +0.05 ))
}
print(b)
pdf(file="plots/pi_vs_mu.pdf",width=5,height=3)
print(b+xlab(TeX(sprintf('$\\mu$')))+ylab(TeX(sprintf('$\\pi$'))) +theme_bw())
dev.off()


## plot for dependence between pi and mu for alpha \in  [0,1]
alpha_colors= viridis(11)
b=ggplot() 
x_tilde =1
mu= seq(-3,3,by =0.05)
alpha_seq=seq(0,0.5,by =0.05)
for(i in 1:length(alpha_seq)){
  b <- b + geom_line(aes(x=mu,y=pi), 
                     data = tibble(mu = mu, pi = alpha_seq[i]*(0.8*pnorm(mu/ (0.1*(2-x_tilde) +0.25) ) +0.025*(x_tilde)) +0.05*(0.5*(19 - 17*alpha_seq[i]))), color=alpha_colors[i])
}

pdf(file="plots/pi_vs_mu_alpha.pdf",width=5,height=3)
print(b+xlab(TeX(sprintf('$\\mu$')))+ylab(TeX(sprintf('$\\pi$'))) +theme_bw())
dev.off()

