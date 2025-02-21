rm(list=ls())
library(R2OpenBUGS)
set.seed(1)
#stat_ed <- bugs.log("D:/My Document/Rcode/degradation-ED/ed_change_point/ed_dep/log.txt")[[1]]
stat_ed <- bugs.log("D:/My Document/Rcode/degradation-ED/ed_change_point/ed_dep/log.txt")[[1]]
stat_two <- read.table("D:/My Document/Rcode/degradation-GP/GPcode3/led_data/beta_two.txt")
stat_bi <- read.table("D:/My Document/Rcode/degradation-GP/GPcode3/led_data/stat_bi.txt")
stat_cpwp <- read.table("D:/My Document/Rcode/degradation-GP/GPcode3/led_data/stat.txt")

## ED
mu_h_ed <- -stat_ed["me[1]",1]
sigma_h_ed <- stat_ed["Sigma[1,1]",1]

mu_l_ed <- -stat_ed["me[2]",1]
sigma_l_ed <- stat_ed["Sigma[2,2]",1]

mu_tau_ed <- stat_ed["me[3]",1]
sigma_tau_ed <- stat_ed["Sigma[3,3]",1]


##Two
mu_h_two <- stat_two["me[1]",1]
sigma_h_two <- stat_two["Sigma[1,1]",1]

mu_l_two <- stat_two["me[2]",1]
sigma_l_two <- stat_two["Sigma[2,2]",1]

mu_tau_two <- stat_two["me[3]",1]
sigma_tau_two <- stat_two["Sigma[3,3]",1]

## BE
mu1 <- stat_bi["me[1]",1]
sigma1 <- stat_bi["Sigma[1,1]",1]

mu2 <- stat_bi["me[2]",1]
sigma2 <- stat_bi["Sigma[2,2]",1]

mu3 <- stat_bi["me[3]",1]
sigma3 <- stat_bi["Sigma[3,3]",1]

## CPWP
mu_h_cpwp <- stat_cpwp["me[1]",1]
sigma_h_cpwp <- stat_cpwp["Sigma[1,1]",1]

mu_l_cpwp <- stat_cpwp["me[2]",1]
sigma_l_cpwp <- stat_cpwp["Sigma[2,2]",1]

mu_tau_cpwp <- stat_cpwp["me[3]",1]
sigma_tau_cpwp <- stat_cpwp["Sigma[3,3]",1]


iter <- 10000
mspe_cpwp <- mspe_ed <- mspe_two <- mspe_bi <- numeric(iter)


old <- read.csv("D:\\My Document\\Rcode\\degradation-GP\\GPcode2\\Bae\\oledalt1.csv")
old <- old[-which(old[,1]==0),]
select <- c(18)
m <- length(select)
data <- old[(old[,3]%in%select),1:2]

iter <- 10000
y_fit_cpwp <- y_fit_ed <- y_fit_two <- y_fit_bi <- matrix(NA, 19, iter)

for(i in 1:iter){
  beta_h_ed <- rnorm(1, mu_h_ed, sqrt(sigma_h_ed))
  beta_l_ed <- rnorm(1,mu_l_ed, sqrt(sigma_l_ed))
  tau_ed <- rnorm(1,mu_tau_ed,sqrt(sigma_tau_ed))
  
  beta_h_two <- rnorm(1, mu_h_two, sqrt(sigma_h_two))
  beta_l_two <- rnorm(1,mu_l_two, sqrt(sigma_l_two))
  tau_two <- rnorm(1,mu_tau_two,sqrt(sigma_tau_two))
  
  gamma1 <- rnorm(1, mu1, sqrt(sigma1))
  gamma2 <- rnorm(1, mu2, sqrt(sigma2))
  gamma3 <- rnorm(1, mu3, sqrt(sigma3))
  
  beta_h_cpwp <- rnorm(1, mu_h_cpwp, sqrt(sigma_h_cpwp))
  beta_l_cpwp <- rnorm(1,mu_l_cpwp, sqrt(sigma_l_cpwp))
  tau_cpwp <- rnorm(1,mu_tau_cpwp,sqrt(sigma_tau_cpwp))
  
  tm_p <- log(data[,1])-log(data[,1])[1]
  tm_bi <- (data[,1]-data[1,1])/100
  for(j in 1:length(tm_p)){
    y_fit_cpwp[j,i] <- beta_h_cpwp*min(tm_p[j],tau_cpwp) + beta_l_cpwp*max(tm_p[j]-tau_cpwp,0)
    y_fit_ed[j,i] <- beta_h_ed*min(tm_p[j],tau_ed) + beta_l_ed*max(tm_p[j]-tau_ed,0)
    y_fit_two[j,i] <- beta_h_two*tm_p[j]-beta_l_two*min(tm_p[j],tau_two) 
    y_fit_bi[j,i] <- log(gamma1*exp(-(gamma2+gamma3)*tm_bi[j])+(1-gamma1)*exp(-gamma2*tm_bi[j]))*50
  }  
  
  mspe_cpwp[i] <- sum((y_fit_cpwp[,i]-data[1:19,2]+data[1,2])^2)
  mspe_ed[i] <- sum((y_fit_ed[,i]-data[1:19,2]+data[1,2])^2)
  mspe_two[i] <-sum((y_fit_two[,i]-data[1:19,2]+data[1,2])^2)
  mspe_bi[i] <- sum((y_fit_bi[,i]-data[1:19,2]+data[1,2])^2)
}

library(xtable)
xtable(cbind(mean(mspe_ed), mean(mspe_cpwp),mean(mspe_two), mean(mspe_bi)), digits=3)

