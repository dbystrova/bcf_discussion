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

x = matrix(runif(n*p, min = 0, max = 1), nrow=n)

mu_function<- function(x1, x2){
  return(-3 + 6*pnorm(2*(x1 - x2)))  
}

pi_function<- function(x1, x2, alpha = 1){
  return( alpha*(0.8*pnorm(mu_function(x1, x2)/ (0.1*(2-x1- x2) +0.25 )) + 0.025*(x1 +x2) ) +0.05*(0.5*(19 - 17*alpha)) ) 
}


##### simulation example #######################################################



simulation <- function(i,alpha ){
  set.seed(i+1)
  #### control variables matrix x_i ~ Unif(0,1)
  x = matrix(runif(n*p, min = 0, max = 1), nrow=n)
  # create targeted selection
  q= -3 + 6*pnorm((2*(x[,1]- x[,2]))) 
  # generate treatment variable
  # probability of recieving the treatment
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
  
  # estimate pihat here
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
  bcf_ate = mean(tauhat)
  #Bart
  bartc1 <- bartc(response = y, treatment = z, confounders = x, method.rsp = "bart", method.trt = "none", p.scoreAsCovariate = FALSE,n.chains=2,n.samples = 2000)
  bart_ate = unlist(summary(bartc1)[9]$estimates[1])
  return( tibble(i=i, bcf_ate = bcf_ate, bart_ate = bart_ate))
}

alpha_seq <- seq(0, 1, by = 0.1)
data_alpha_list = list()
for (j in 1:length(alpha_seq)){
  datalist = list()
  for (i in 1:200){
   dat <- simulation(j*50 + i, alpha =alpha_seq[j])
   dat$alpha = alpha_seq[j]
   datalist[[i]] <- as.data.frame(dat)
   print(c(i,j))
  }
  df = do.call(rbind, datalist)
  dat_alpha = tibble(alpha =alpha_seq[j],
                     fin_rmse_bart = sqrt(mean((df$bart_ate + 1)^2)),
                     fin_rmse_bcf = sqrt(mean((df$bcf_ate  + 1)^2)),
                     bias_bart= mean(df$bart_ate + 1),
                     bias_bcf =  mean(df$bcf_ate  + 1),
                     ate_bart_q10 = quantile(df$bart_ate, 0.1),
                     ate_bart_q90 = quantile(df$bart_ate, 0.9),
                     ate_bcf_q10 = quantile(df$bcf_ate, 0.1),
                     ate_bcf_q90 = quantile(df$bcf_ate, 0.9))
  data_alpha_list[[j]]<- dat_alpha
}

df_alpha_ <-  do.call(rbind, data_alpha_list)
#save(df_alpha_, file = "alpha_100.Rdata")
pdf(file="plots/bias.pdf",width=5,height=3)
df_alpha_[, c("alpha","bias_bart", "bias_bcf")]%>% gather(Model, bias, bias_bart:bias_bcf)%>%
  ggplot(aes(x=alpha,y=bias,col=Model))+geom_line(alpha=0.7)+ scale_color_viridis(discrete=TRUE)+
  xlab(TeX(sprintf('$\\alpha$')))+ylab("Bias") +theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 12),axis.text.y = element_text(angle = 0, hjust = 1,size = 12),legend.position = "none", plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(size = rel(1.5), angle = 90),
        axis.title.x = element_text(size = rel(1.5), angle = 00))

dev.off()

pdf(file="plots/rmse.pdf",width=5,height=3)
df_alpha_[, c("alpha","fin_rmse_bart", "fin_rmse_bcf")]%>% gather(Model, Rmse, fin_rmse_bart:fin_rmse_bcf)%>%
  ggplot(aes(x=alpha,y=Rmse,col=Model))+geom_line(alpha=0.7)+ scale_color_viridis(discrete=TRUE)+
  #labs(title="Bias")+
  xlab(TeX(sprintf('$\\alpha$')))+ylab("RMSE") +theme_bw()+ theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 12),axis.text.y = element_text(angle = 0, hjust = 1,size = 12),legend.position = "none", plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(size = rel(1.5), angle = 90),
        axis.title.x = element_text(size = rel(1.5), angle = 00))
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
# 
pdf(file="plots/bias_ribbon.pdf",width=5,height=3)
bias_all <- cbind(df_bias, dq10[, "q10"], dq90[, "q90"])
bias_all%>% ggplot(aes(x=alpha,y=bias,col=Model))+geom_line(alpha=0.9, size=1)+ geom_ribbon(aes(ymin = q10 + 1, ymax = q90 + 1),linetype = 2,alpha=0.0) +
  xlab(TeX(sprintf('$\\alpha$')))+ylab("bias") +theme_bw() +theme(legend.position="none")
dev.off()


