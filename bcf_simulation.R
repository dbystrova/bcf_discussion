# data generating process
#p = 2 #two control variables
library(tidyverse)
library("bcf")
library(latex2exp)
library(viridis)

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


pdf(file="mu.pdf",width=5,height=3)
build3ds1(x1,x2,mu_,par1= expression(paste(italic(mu))) )
dev.off()
pdf(file="pi_1.pdf",width=3,height=2)
build3ds1(x1,x2,pi_1,z_lim=c(0,1),par1= expression(paste(italic(pi))) )
dev.off()
pdf(file="pi_05.pdf",width=3,height=2)
build3ds1(x1,x2,pi_05,z_lim=c(0,1),par1= expression(paste(italic(pi))) )
dev.off()

pdf(file="pi_0.pdf",width=3,height=2)
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

b=ggplot() 

#x_tilde = x[,1]+  x[,2]
x_tilde =1
mu= seq(-3,3,by =0.05)
alpha_seq=seq(0,0.5,by =0.05)
for(i in 1:length(alpha_seq)){
  b <- b + geom_line(aes(x=mu,y=pi), 
                     data = tibble(mu = mu, pi = alpha_seq[i]*(0.8*pnorm(mu/ (0.1*(2-x_tilde) +0.25) ) +0.025*(x_tilde)) +0.05*(0.5*(19 - 17*alpha_seq[i]))))
}
pdf(file="pi_vs_mu_alpha.pdf",width=5,height=3)
print(b)
dev.off()


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
  tau_ate <- data.frame(Mean  = mean(colMeans(bcf_fit$tau)),
                         Low95 = quantile(rowMeans(bcf_fit$tau), 0.025),
                         Up95  = quantile(rowMeans(bcf_fit$tau), 0.975))
  
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
  coverage_bcf =  ifelse(tau_ate$Low95<= tau & tau <=tau_ate$Up95, 1, 0)
  #Bart
  bartc1 <- bartc(response = y, treatment = z, confounders = x, method.rsp = "bart", method.trt = "none", n.samples = 3000)
  bart_ate = unlist(summary(bartc1)[9]$estimates[1])
  #ites <- extract(bartc1, type = "ite")
  #tau_x  = apply(ites, 2, mean)
  #RMSE_bart= sqrt(sum((1 - tau_x)^2)/n)
  #RMSE_bart
  coverage_bart =  as.numeric(ifelse(summary(bartc1)[9]$estimates[3]<= tau & tau <=summary(bartc1)[9]$estimates[4], 1, 0))
  return( tibble(i=i, bcf_ate =bcf_ate,bcf_coverage =coverage_bcf,bcf_cover =perCoverage, bcf_bias = bias, bcf_rmse= RMSE,bart_ate = bart_ate, bart_coverage = coverage_bart))
 # return(list(bcf = list(bcf_ate =bcf_ate,bcf_coverage =coverage_bcf,bcf_cover =perCoverage, bcf_bias = bias, bcf_rmse= RMSE), bart= list(bart_ate = bart_ate, bart_coverage = coverage_bart))) 
}

alpha_seq <- seq(0, 1, by = 0.2)
data_alpha_list = list()
for (j in 1:length(alpha_seq)){
  datalist = list()
  for (i in 1:20){
   dat <- simulation(i, alpha =alpha_seq[j])
   dat$alpha = alpha_seq[j]
   #dat =  column_to_rownames(dat, var = "i")
   rownames(dat) <- NULL
   datalist[[i]] <- as.data.frame(dat)
  }
  df = do.call(rbind, datalist)
  dat_alpha = tibble(alpha =alpha_seq[j],
                    fin_rmse_bart = sqrt(mean((unlist(df$bart_ate ) - 1)^2)),
                    fin_rmse_bcf = sqrt(mean((df$bcf_ate  - 1)^2)),
                    bias_bart= mean(unlist(df$bart_ate ) - 1),
                    bias_bcf =  mean(df$bcf_ate  - 1),
                    bart_coverage = sum(df$bart_coverage)/length(df$bart_coverage),
                    bcf_coverage = sum(df$bcf_coverage)/length(df$bcf_coverage))
  data_alpha_list[[j]]<- dat_alpha
}

df_alpha_ <-  do.call(rbind, data_alpha_list)

df_alpha_[, c("alpha","bias_bart", "bias_bcf")]%>% gather(Model, bias, bias_bart:bias_bcf)%>%
  ggplot(aes(x=alpha,y=bias,col=Model))+geom_line(alpha=0.7)+ scale_color_viridis(discrete=TRUE)+
  #labs(title="Rmse")+
  xlab("iterations")+ylab("bias") +theme_bw()

  


df_alpha_[, c("alpha","fin_rmse_bart", "fin_rmse_bcf")]%>% gather(Model, Rmse, fin_rmse_bart:fin_rmse_bcf)%>%
  ggplot(aes(x=alpha,y=Rmse,col=Model))+geom_line(alpha=0.7)+ scale_color_viridis(discrete=TRUE)+
  #labs(title="Bias")+
  xlab("iterations")+ylab("RMSE") +theme_bw()


















# 
# 
# 
# bart_rmse<- c()
# bart_rmse<- c()
# 
# alpha_seq <- seq(0, 1, by = 0.1)
# for (j in 1:11){
#   bart_ate_<-c()
#   bcf_ate_<- c()
#   bcf_c <- c()
#   bart_c<- c()
#   for (i in 1:10){
#     sim= simulation(i, alpha =alpha_seq[j])
#     bart_ate_ <- c(bart_ate_, sim$bart$bart_ate)
#     bcf_ate_ <- c(bcf_ate_, sim$bcf$bcf_ate)
#     bart_c <-c(bart_c, sim$bart$bart_coverage) 
#     bcf_c <-c(bcf_c, sim$bcf$bcf_coverage) 
#   }
# fin_rmse_bart = sqrt(mean((unlist(bart_ate_)  - 1)^2))
# fin_rmse_bcf = sqrt(mean((unlist(bcf_ate_)  - 1)^2))
# 
# fin_rmse_bart
# fin_rmse_bcf
# 
# bias_bart= mean((  - 1))
# bias_bcf =  mean((unlist(bcf_ate_)  - 1))
# bias_bart
# bias_bcf
# 
# bart_coverage = sum(unlist(bart_c))/length(bart_c)
# bart_coverage
# 
# bcf_coverage = sum(unlist(bcf_c))/length(bcf_c)
# bcf_coverage
# 
# 
# }
# 
