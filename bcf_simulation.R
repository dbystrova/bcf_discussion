# data generating process
#p = 2 #two control variables
library(tidyverse)
library(bcf)
library(latex2exp)
library(viridis)
library(bartCause)
source('~/Documents/GitHub/bcf_discussion/functions.R')
p=2
n = 250
alpha = 1
n_burn <- 500
n_sim <- 2000

x = matrix(runif(n*p, min = 0, max = 1), nrow=n)

mu_function<- function(x1, x2){
  return(-3 + 6*pnorm(2*(x1 - x2)))  
}

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


pdf(file="mu.pdf",width=5,height=3)
build3ds1(x1,x2,mu_,par1= expression(paste(italic(mu))) )
dev.off()
pdf(file="pi_1.pdf",width=10,height=9)
build3ds1(x1,x2,pi_1,z_lim=c(0,1),par1= expression(paste(italic(alpha)," = 1")) )
dev.off()
pdf(file="pi_05.pdf",width=10,height=9)
build3ds1(x1,x2,pi_05,z_lim=c(0,1),par1= expression(paste(italic(alpha)," = 0.5")))
dev.off()

pdf(file="pi_0.pdf",width=10,height=9)
build3ds1(x1,x2,pi_0,z_lim=c(0,1),par1= expression(paste(italic(alpha)," = 0")))
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

alpha_colors= viridis(11)
b=ggplot() 
#x_tilde = x[,1]+  x[,2]
x_tilde =1
mu= seq(-3,3,by =0.05)
alpha_seq=seq(0,0.5,by =0.05)
for(i in 1:length(alpha_seq)){
  b <- b + geom_line(aes(x=mu,y=pi), 
                     data = tibble(mu = mu, pi = alpha_seq[i]*(0.8*pnorm(mu/ (0.1*(2-x_tilde) +0.25) ) +0.025*(x_tilde)) +0.05*(0.5*(19 - 17*alpha_seq[i]))), color=alpha_colors[i])
}

pdf(file="pi_vs_mu_alpha.pdf",width=5,height=3)
print(b+xlab(TeX(sprintf('$\\mu$')))+ylab(TeX(sprintf('$\\pi$'))) +theme_bw())
dev.off()


simulation <- function(i,alpha ){
  set.seed(i+1)
  #### control variables matrix x_i ~ Unif(0,1)
  x = matrix(runif(n*p, min = 0, max = 1), nrow=n)
  # create targeted selection
  q= -3 + 6*pnorm((2*(x[,2]- x[,1]))) 
  # generate treatment variable
  # probability of recievng the treatment
  pi = alpha*(0.8* pnorm(q/ (0.1*(2-x[,1]- x[,2]) +0.25 )) + 0.025*(x[,1] +x[,2])) +0.05*(0.5*(19 - 17*alpha))
  #treatment
  z = rbinom(n,1,pi)
  # tau is the true (homogeneous) treatment effect
  tau = 1
  # generate the response using q, tau and z
  mu = (q - tau*z)
  # set the noise level relative to the expected mean function of Y
  #sigma = diff(range(q + tau*pi))/8
  sigma=1
  # draw the response variable with additive error
  y = mu + sigma*rnorm(n)
  
  # If you didn't know pi, you would estimate it here
  #pihat = pnorm(q)
  pihat = pi

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
  RMSE= sqrt(sum((1+ tauhat)^2)/n)
  ## bias
  bias= mean(tauhat + 1)
  ## coverage
  isCovered <- function(i){
    ifelse(tau_ests$Low95[i] <= tau & tau <= tau_ests$Up95[i], 1, 0)
  }
  coverage <- lapply(1:length(tau_ests$Mean), isCovered)
  perCoverage <- sum(unlist(coverage))/length(tau_ests)
  bcf_ate = mean(tauhat)
  coverage_bcf =  ifelse(tau_ate$Low95<= -tau & -tau <=tau_ate$Up95, 1, 0)
  #Bart
  bartc1 <- bartc(response = y, treatment = z, confounders = x, method.rsp = "bart", method.trt = "none", p.scoreAsCovariate = FALSE,n.chains=2,n.samples = 2000)
  bart_ate = unlist(summary(bartc1)[9]$estimates[1])
  #ites <- extract(bartc1, type = "ite")
  #tau_x  = apply(ites, 2, mean)
  #RMSE_bart= sqrt(sum((1 - tau_x)^2)/n)
  #RMSE_bart
  coverage_bart =  as.numeric(ifelse(summary(bartc1)[9]$estimates[3]<= -tau & -tau <=summary(bartc1)[9]$estimates[4], 1, 0))
  return( tibble(i=i, bcf_ate =bcf_ate,bcf_coverage =coverage_bcf,bcf_cover =perCoverage, bcf_bias = bias, bcf_rmse= RMSE,bart_ate = bart_ate, bart_coverage = coverage_bart))
 # return(list(bcf = list(bcf_ate =bcf_ate,bcf_coverage =coverage_bcf,bcf_cover =perCoverage, bcf_bias = bias, bcf_rmse= RMSE), bart= list(bart_ate = bart_ate, bart_coverage = coverage_bart))) 
}

alpha_seq <- seq(0, 1, by = 0.2)
data_alpha_list = list()
for (j in 1:length(alpha_seq)){
  datalist = list()
  for (i in 1:100){
   dat <- simulation(i, alpha =alpha_seq[j])
   dat$alpha = alpha_seq[j]
   #dat =  column_to_rownames(dat, var = "i")
   rownames(dat) <- NULL
   datalist[[i]] <- as.data.frame(dat)
   print(c(i,j))
  }
  df = do.call(rbind, datalist)
  dat_alpha = tibble(alpha =alpha_seq[j],
                     fin_rmse_bart = sqrt(mean((unlist(df$bart_ate ) + 1)^2)),
                     fin_rmse_bcf = sqrt(mean((df$bcf_ate  + 1)^2)),
                     bias_bart= mean(unlist(df$bart_ate ) + 1),
                     bias_bcf =  mean(df$bcf_ate  + 1),
                     ate_bart_q10 = quantile(df$bart_ate, 0.1),
                     ate_bart_q90 = quantile(df$bart_ate, 0.9),
                     ate_bcf_q10 = quantile(df$bcf_ate, 0.1),
                     ate_bcf_q90 = quantile(df$bcf_ate, 0.9))
  data_alpha_list[[j]]<- dat_alpha
}

df_alpha_ <-  do.call(rbind, data_alpha_list)

pdf(file="bias.pdf",width=5,height=3)
df_alpha_[, c("alpha","bias_bart", "bias_bcf")]%>% gather(Model, bias, bias_bart:bias_bcf)%>%
  ggplot(aes(x=alpha,y=bias,col=Model))+geom_line(alpha=0.7)+ scale_color_viridis(discrete=TRUE)+
  #labs(title="Rmse")+
  xlab(TeX(sprintf('$\\alpha$')))+ylab("bias") +theme_bw() +theme(legend.position="none")
dev.off()

pdf(file="rmse.pdf",width=5,height=3)
df_alpha_[, c("alpha","fin_rmse_bart", "fin_rmse_bcf")]%>% gather(Model, Rmse, fin_rmse_bart:fin_rmse_bcf)%>%
  ggplot(aes(x=alpha,y=Rmse,col=Model))+geom_line(alpha=0.7)+ scale_color_viridis(discrete=TRUE)+
  #labs(title="Bias")+
  xlab(TeX(sprintf('$\\alpha$')))+ylab("RMSE") +theme_bw()+ theme(legend.position="none")
dev.off()



df_bias <- df_alpha_[, c("alpha","bias_bart", "bias_bcf")]%>% gather(Model, bias, bias_bart:bias_bcf)
dq10 <- df_alpha_[, c("alpha","ate_bart_q10", "ate_bcf_q10")]%>% gather(Model, q10, ate_bart_q10:ate_bcf_q10)
dq90 <- df_alpha_[, c("alpha","ate_bart_q90", "ate_bcf_q90")]%>% gather(Model, q90, ate_bart_q90:ate_bcf_q90)
for (i in 1:(12)){
  if (grepl("bart",dq90$Model[i], fixed = TRUE)) {dq90$Model[i] = "bart"}
  if (grepl("bcf",dq90$Model[i], fixed = TRUE)){dq90$Model[i] = "bcf"}
  if (grepl("bart",dq10$Model[i], fixed = TRUE)){dq10$Model[i] = "bart"}
  if (grepl("bcf",dq10$Model[i], fixed = TRUE)){dq10$Model[i] = "bcf"}
  if (grepl("bart",df_bias$Model[i], fixed = TRUE)){df_bias$Model[i] = "bart"}
  if (grepl("bcf",df_bias$Model[i], fixed = TRUE)){df_bias$Model[i] = "bcf"}
  
}

pdf(file="bias_ribbon.pdf",width=5,height=3)
bias_all <- cbind(df_bias, dq10[, "q10"], dq90[, "q90"])
bias_all%>% ggplot(aes(x=alpha,y=bias,col=Model))+geom_line(alpha=0.9, size=1)+ geom_ribbon(aes(ymin = q10 + 1, ymax = q90 + 1),linetype = 2,alpha=0.0) +
  xlab(TeX(sprintf('$\\alpha$')))+ylab("bias") +theme_bw() +theme(legend.position="none")
dev.off()





# 
# library(dbarts)
# bartFit<- bart(x, y, nskip = 350, ntree = 1000)
# plot(bartFit)
# 
# df2 <- data.frame(df, 
#                   ql = apply(bartFit$yhat.train, length(dim(bartFit$yhat.train)), quantile,probs=0.05),
#                   qm = apply(bartFit$yhat.train, length(dim(bartFit$yhat.train)), quantile,probs=.5),
#                   qu <- apply(bartFit$yhat.train, length(dim(bartFit$yhat.train)), quantile,probs=0.95)
# )
# 
# bartp <- ggplot(df2, aes(x= y, y = qm)) + geom_linerange(aes(ymin = ql, ymax = qu), col = "grey") +
#   geom_point() + geom_smooth() +
#   geom_abline(intercept = 0, slope = 1, col = "red", size = 1)
# 
# bartp
# 

