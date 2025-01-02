rm(list=ls())
library(R2OpenBUGS)
library(coda)
library(mvtnorm)
library(parallel)
library(statmod)
library(ggplot2)
set.seed(1)
sim <- 100

#set parameters' value
mu_tau <- 15
sigma2_tau <- 0.3

mu_h <- 7
sigma2_h <- 0.2

mu_l <- 2
sigma2_l <- 0.1

lambda  <- 10


I <- 10
n_i <- 21



##Data
#set.seed(1)
beta_h <- rnorm(I, mu_h, sqrt(sigma2_h))
#beta_h <- c(7.608,6.954,7.173,6.976,6.384,6.814,6.824,6.973,7.492,7.341)
beta_l <- rnorm(I, mu_l, sqrt(sigma2_l))
#beta_l <- c(2.126,1.806,2.108,1.643,2.453,2.626,1.884,1.67,2.18,1.957)
tau <- rnorm(I, mu_tau, sqrt(sigma2_tau))
#tau <- c(13.315,11.979,12.378,12.015,11.593,12.103,11.011,12.803,12.084,13.19)
#Sample data

time <- seq(0, 40,length=n_i)
t <- matrix(0,nrow=I, ncol=n_i)
for(i in 1:I){
  t[i,] <- time
}

delta_t <- matrix(0, nrow=I, ncol=(n_i-1))
for(i in 1:I){
  for(j in 1:(n_i-1)){
    delta_t[i,j] <- t[i,(j+1)]-t[i,j]
  }
}

t_up <- t[,-1]
t_low <- t[,-21]

Mean <- Shape <- matrix(nrow=I,ncol=(n_i-1))
for(i in 1:I){
  for(j in 1:(n_i-1)){
    Shape[i,j] <- lambda*(delta_t[i,j]^2)
  }
}

for(i in 1:I){
  for(j in 1:(n_i-1)){
    if(tau[i]>=t_up[i,j]){
      Mean[i,j] <- beta_h[i]*delta_t[i,j]
    } else{
      if(tau[i]>=t_low[i,j]){
        weight <- (tau[i]-t_low[i,j])/(t_up[i,j]-t_low[i,j])
        Mean[i,j] <- (weight*beta_h[i]+(1-weight)*beta_l[i])*delta_t[i,j]
      } else{
        Mean[i,j] <- beta_l[i]*delta_t[i,j]
      }
    }
  }
} 




for(g in 1:sim){
  set.seed(g)
  delta_y <- matrix(nrow=I,ncol=(n_i-1))
  for(i in 1:I){
    for(j in 1:(n_i-1)){
      delta_y[i,j] <- rinvgauss(1, mean=Mean[i,j], shape = Shape[i,j] )
    }
  }
  delta_y[delta_y<=0] <- abs(delta_y[delta_y<=0])
  
  #Y
  y <- matrix(nrow=I,ncol=n_i)
  y[, 1] <- rep(0, I)
  
  for(i in 1:I){
    y[i,2:n_i] <- y[i,1]+cumsum(delta_y[i,])
  }
  
  
  t_up <- t[,-1]
  t_low <- t[,-21]
  
  m <- I
  ni <- n_i-1
  R <- diag(3)
  mm <- rep(0,3)
  pm <- diag(rep(0.000001,3))
  L <- c(0,0, 0)
  U <- c(16,20,30)
  data <- list("m", "ni", "R", "mm", "pm", "L", "U","delta_t","delta_y","t_up","t_low")
  inits <- function (){
    list(lambda = 10,
         beta = matrix(rep(c(7,2, 15),m),m, 3, byrow = T),
         me = c(7,2, 15),
         Q=diag(3),
         xi=rep(1,3))
    list(lambda = 10,
         beta = matrix(rep(c(8,3, 14),m),m, 3, byrow = T),
         me = c(8,3, 14),
         Q=diag(3),
         xi=rep(1,3))
  }
  parameters <- c("beta","me","lambda")
  
  dir.create(paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/ig/bugs_wie",g, sep=""))
  m_wiener <- bugs(data,inits, parameters,
                   working.directory = paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/ig/","bugs_wie",g,sep=""),
                   model.file="D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/ig/wiener.odc",
                   n.chains=2, n.iter=20000, codaPkg=T,  n.burnin=12000,  debug=F, save.history = T)
  
  dir.create(paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/ig/bugs_gamma",g, sep=""))
  m_gamma <- bugs(data,inits, parameters,
                  working.directory = paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/ig/","bugs_gamma",g,sep=""),
                  model.file="D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/ig/gamma.odc",
                  n.chains=2, n.iter=20000, codaPkg=T,  n.burnin=15000,  debug=F, save.history = T)
  
  
  dir.create(paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/ig/bugs_ig",g, sep=""))
  m_ig <- bugs(data,inits, parameters,
                 working.directory = paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/ig/","bugs_ig",g,sep=""),
                 model.file="D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/ig/ig.odc",
                 n.chains=2, n.iter=20000, codaPkg=T,  n.burnin=15000,  debug=F, save.history = T)
  
  parameters <- c("beta","me","lambda","p")
  dir.create(paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/ig/bugs_ed",g, sep=""))
  m_ed_dep <- bugs(data, inits=inits, parameters,
                   working.directory = paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/ig/","bugs_ed",g,sep=""),
                   model.file="D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/ig/ED_two_dep.odc",
                   n.chains=2, n.iter=20000, n.burnin=5000,  
                   digits=5, codaPkg=T, debug=F, save.history = T, n.thin=2)
  #estimate_ed <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/ig/","bugs_ed",g, "/", "log.txt", sep=""))[[1]]
  
  cat("----------------g = ", g,"----------------\n")
}

degdata <- data.frame()
for(i in 1:I){
  degdata <- rbind(degdata, data.frame(time=t[i,], deg=y[i,], 
                                       Unit=paste("Unit #", i)))
}

deg <- ggplot(data=degdata, aes(x=time, y=deg, linetype = Unit, color=Unit))+
  geom_line(size=0.4)+geom_point(size=0.8)+xlab("Time")+ylab("Degradation Value")+
  theme(axis.text=element_text(size=8), axis.title=element_text(size = 8),
        legend.key.size=unit(0.3,'cm'),legend.text = element_text(size = 7),
        legend.title=element_blank())


