rm(list=ls())
library(R2OpenBUGS)
##BS distribtion
library(VGAM)
## IG distribution
library(statmod)
library(ggplot2)
library(plotly)
library(xtable)


beta_h_true <- c(6.719841, 7.082128, 6.626296, 7.713431, 7.147360, 6.633075, 7.217985, 7.330189, 7.257497, 6.863426)
beta_l_true <- c(2.478067, 2.123279, 1.803546, 1.299650, 2.355734, 1.985791, 1.994880, 2.298467, 2.259693, 2.187808)
tau_true <- c(15.50334, 15.42839, 15.04084, 13.91039, 15.33949, 14.96926, 14.91467, 14.19444, 14.73811, 15.22892)

#### Wiener process data

I=10
Df <- rep(200, I)

g=1
est_wie <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point_simulation/simulation/wiener/","bugs_wie",g, "/", "log.txt", sep=""))[[1]]
est_gamma <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point_simulation/simulation/wiener/","bugs_gamma",g, "/", "log.txt", sep=""))[[1]]

est_ig <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point_simulation/simulation/wiener/","bugs_ig",g, "/", "log.txt", sep=""))[[1]]
est_ed <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point_simulation/simulation/wiener/","bugs_ed",g, "/", "log.txt", sep=""))[[1]]

beta_h_wie <- est_wie[(0:9)*3+1,1] 
beta_l_wie <- est_wie[(0:9)*3+2,1] 
tau_wie <- est_wie[(0:9)*3+3,1] 
lambda_wie <- est_wie[31,1] 
p_wie=0

beta_h_gamma <- est_gamma[(0:9)*3+1,1] 
beta_l_gamma <- est_gamma[(0:9)*3+2,1] 
tau_gamma <- est_gamma[(0:9)*3+3,1] 
lambda_gamma <- est_gamma[31,1] 
p_gamma <- 2 

beta_h_ig <- est_ig[(0:9)*3+1,1] 
beta_l_ig <- est_ig[(0:9)*3+2,1] 
tau_ig <- est_ig[(0:9)*3+3,1] 
lambda_ig <- est_ig[31,1] 
p_ig=3

beta_h_ed <- est_ed[(0:9)*3+1,1] 
beta_l_ed <- est_ed[(0:9)*3+2,1] 
tau_ed <- est_ed[(0:9)*3+3,1] 
lambda_ed <- est_ed[31,1] 
p_ed <- est_ed[35,1] 

# ft1 <- function(t, ff, beta_h, beta_l, tau, lambda, p){
#   if(t<tau){
#     ft <- dbisa(t, scale=ff/beta_h, shape=beta_h^((p-1)/2)/sqrt(lambda*ff))
#   } else {
#     int_fun <- function(y){
#       temp1 <- dbisa(t-tau, scale=(ff-y)/beta_l, shape=beta_l^((p-1)/2)/sqrt(lambda*(ff-y)))
#       temp2 <- (y/tau)^(2-p)/(1-p)/(2-p)-y*beta_h^(1-p)/tau/(1-p)+beta_h^(2-p)/(2-p)
#       temp3 <- sqrt(lambda/2/pi/tau^(1-p)/y^p)*exp(-lambda*tau*temp2)
#       return(temp1*temp3)
#     }
#     ft <- integrate(int_fun, lower = 0.001, upper=ff-1)$value
#     
#   }
#   return(ft)
# }

ft_pdf_ig <- function(t, ff, beta_h, beta_l, tau, sigma2){
  if(t<tau){
    ft <- dinvgauss(t, mean=ff/beta_h, shape=ff^2/sigma2)
  } else ft <- dinvgauss(t, mean=(ff-(beta_h-beta_l)*tau)/beta_l, 
                         shape = (ff-(beta_h-beta_l)*tau)^2/sigma2)
  return(ft)
}

ft_cdf_ig <- function(t, ff, beta_h, beta_l, tau, sigma2){
  if(t<tau){
    ft <- pinvgauss(t, mean=ff/beta_h, shape=ff^2/sigma2)
  } else ft <- pinvgauss(t, mean=(ff-(beta_h-beta_l)*tau)/beta_l, 
                         shape = (ff-(beta_h-beta_l)*tau)^2/sigma2)
  return(ft)
}


ft_pdf_bs <- function(t, ff, beta_h, beta_l, tau, lambda, p){
  if(t<tau){
    ft <- dbisa(t, scale=ff/beta_h, shape=beta_h^((p-1)/2)/sqrt(lambda*ff))
  } else {
    ft <- dbisa(t, scale=(ff-(beta_h-beta_l)*tau)/beta_l, shape=beta_l^((p-1)/2)/sqrt(lambda*((ff-(beta_h-beta_l)*tau))))
  }
  return(ft)
}

ft_cdf_bs <- function(t, ff, beta_h, beta_l, tau, lambda, p){
  if(t<tau){
    ft <- pbisa(t, scale=ff/beta_h, shape=beta_h^((p-1)/2)/sqrt(lambda*ff))
  } else {
    ft <- pbisa(t, scale=(ff-(beta_h-beta_l)*tau)/beta_l, shape=beta_l^((p-1)/2)/sqrt(lambda*((ff-(beta_h-beta_l)*tau))))
  }
  return(ft)
}

## Fail-time distribution 
time <- seq(20, 120, length =1000)

fpdf_res_wie <-  fpdf_res_wie_ig <- fpdf_res_gamma <-fpdf_res_ig <-  fpdf_res_ed <- matrix(NA, I, length(time))
for(k in 1:I){
  for(i in 1:length(time)){
    fpdf_res_wie[k,i] <- ft_pdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_wie[k], beta_l = beta_l_wie[k], 
                               tau = tau_wie[k], lambda = lambda_wie, p=p_wie)
    fpdf_res_wie_ig[k,i] <- ft_pdf_ig(t=time[i], ff=Df[k], beta_h = beta_h_wie[k], beta_l = beta_l_wie[k], 
                                  tau = tau_wie[k], sigma2 = 1/lambda_wie)
    fpdf_res_gamma[k,i] <- ft_pdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_gamma[k], beta_l = beta_l_gamma[k], 
                                 tau = tau_gamma[k], lambda = lambda_gamma, p=p_gamma)
    fpdf_res_ig[k,i] <- ft_pdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_ig[k], beta_l = beta_l_ig[k], 
                              tau = tau_ig[k], lambda = lambda_ig, p=p_ig)
    fpdf_res_ed[k,i] <- ft_pdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_ed[k], beta_l = beta_l_ed[k], 
                              tau = tau_ed[k], lambda = lambda_ed, p=p_ed)
  }
}

fcdf_res_wie <-  fcdf_res_wie_ig <- fcdf_res_gamma <-fcdf_res_ig <-  fcdf_res_ed <- matrix(NA, I, length(time))
for(k in 1:I){
  for(i in 1:length(time)){
    fcdf_res_wie[k,i] <- ft_cdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_wie[k], beta_l = beta_l_wie[k], 
                                   tau = tau_wie[k], lambda = lambda_wie, p=p_wie)
    fcdf_res_wie_ig[k,i] <- ft_cdf_ig(t=time[i], ff=Df[k], beta_h = beta_h_wie[k], beta_l = beta_l_wie[k], 
                                      tau = tau_wie[k], sigma2 = 1/lambda_wie)
    fcdf_res_gamma[k,i] <- ft_cdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_gamma[k], beta_l = beta_l_gamma[k], 
                                     tau = tau_gamma[k], lambda = lambda_gamma, p=p_gamma)
    fcdf_res_ig[k,i] <- ft_cdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_ig[k], beta_l = beta_l_ig[k], 
                                  tau = tau_ig[k], lambda = lambda_ig, p=p_ig)
    fcdf_res_ed[k,i] <- ft_cdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_ed[k], beta_l = beta_l_ed[k], 
                                  tau = tau_ed[k], lambda = lambda_ed, p=p_ed)
  }
}

faildata <- data.frame()
for(k in 1:9){
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_wie[k,], ft_cdf=fcdf_res_wie[k,], GaAs=paste("Unit # 0", k, sep = ""), model="CPWP"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_wie_ig[k,], ft_cdf=fcdf_res_wie_ig[k,], GaAs=paste("Unit # 0", k, sep = ""), model="CPWP(IG)"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_gamma[k,], ft_cdf=fcdf_res_gamma[k,], GaAs=paste("Unit # 0", k, sep = ""), model="CPGP"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_ig[k,], ft_cdf=fcdf_res_ig[k,], GaAs=paste("Unit # 0", k, sep = ""), model="CPIGP"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_ed[k,], ft_cdf=fcdf_res_ed[k,], GaAs=paste("Unit # 0", k, sep = ""), model="CPTEDP"))
}

for(k in 10:I){
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_wie[k,], ft_cdf=fcdf_res_wie[k,], GaAs=paste("Unit #", k), model="CPWP"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_wie_ig[k,], ft_cdf=fcdf_res_wie_ig[k,], GaAs=paste("Unit #", k), model="CPWP(IG)"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_gamma[k,], ft_cdf=fcdf_res_gamma[k,], GaAs=paste("Unit #", k), model="CPGP"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_ig[k,], ft_cdf=fcdf_res_ig[k,], GaAs=paste("Unit #", k), model="CPIGP"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_ed[k,], ft_cdf=fcdf_res_ed[k,], GaAs=paste("Unit #", k), model="CPTEDP"))
}


ft_pdf_plot <- ggplot(data=faildata, aes(x=ft, y= ft_pdf, linetype=model, color=model))+facet_wrap(~GaAs, nrow=2)+
  xlim(40,100)+ylim(0,0.3)+
  geom_line(show.legend = T, size=1)+xlab("Failure-time")+ylab("PDF")+ theme_bw()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size = 8),
        legend.key.size=unit(0.3,'cm'),legend.text = element_text(size = 7),
        legend.title=element_blank())

ft_cdf_plot <- ggplot(data=faildata, aes(x=ft, y= ft_cdf, linetype=model, color=model))+facet_wrap(~GaAs, nrow=2)+xlim(40,100)+ylim(0,1)+
  geom_line(show.legend = T, size=1)+xlab("Failure-time")+ylab("CDF")+ theme_bw()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size = 8),
        legend.key.size=unit(0.3,'cm'),legend.text = element_text(size = 7),
        legend.title=element_blank())

postscript("D:/My Document/2023/paper/ED_change_point/figure/simulation/fail_pdf_wie.eps", width=10, height=5)
ft_pdf_plot
dev.off()

postscript("D:/My Document/2023/paper/ED_change_point/figure/simulation/fail_cdf_wie.eps", width=10, height=5)
ft_cdf_plot
dev.off()

## RUL
i=2
z <- seq(0, 40, by = 5)
rul <- seq(0, 70, length=300)
p_true=0
lambda_true=4

rul_pdf_bs <- function(x, t, ff, beta_h, beta_l, tau, lambda, p){
  if(t<tau){
    if(t+x<tau){
      ft <- dbisa(x, scale=(ff-beta_h*t)/beta_h, shape=beta_h^((p-1)/2)/sqrt(lambda*(ff-beta_h*t)))
    } else {
      ft <- dbisa(x, scale=(ff-(beta_h-beta_l)*tau-beta_l*t)/beta_l, 
                  shape=beta_l^((p-1)/2)/sqrt(lambda*((ff-(beta_h-beta_l)*tau-beta_l*t))))
    }
  }else{
    ft <- dbisa(x, scale=(ff-(beta_h-beta_l)*tau-beta_l*t)/beta_l, 
                shape=beta_l^((p-1)/2)/sqrt(lambda*((ff-(beta_h-beta_l)*tau-beta_l*t))))
  }
  return(ft)
}

rul_cdf_bs <- function(x, t, ff, beta_h, beta_l, tau, lambda, p){
  if(t<tau){
    if(t+x<tau){
      ft <- pbisa(x, scale=(ff-beta_h*t)/beta_h, shape=beta_h^((p-1)/2)/sqrt(lambda*(ff-beta_h*t)))
    } else {
      ft <- pbisa(x, scale=(ff-(beta_h-beta_l)*tau-beta_l*t)/beta_l, 
                  shape=beta_l^((p-1)/2)/sqrt(lambda*((ff-(beta_h-beta_l)*tau-beta_l*t))))
    }
  }else{
    ft <- pbisa(x, scale=(ff-(beta_h-beta_l)*tau-beta_l*t)/beta_l, 
                shape=beta_l^((p-1)/2)/sqrt(lambda*((ff-(beta_h-beta_l)*tau-beta_l*t))))
  }
  return(ft)
}

rul3d_df_ed <- rul3d_df_wie <- rul3d_df_ig <- rul3d_df_gamma <-  rul3d_df_true <- data.frame()

for(k in 1:length(z)){
  for(g in 1:length(rul)){
    dens <- rul_pdf_bs(x=rul[g], t=z[k], ff=Df[1], beta_h = beta_h_ed[1], beta_l = beta_l_ed[1], 
                       tau = tau_ed[1], lambda = lambda_ed, p=p_ed)
    dis <- rul_cdf_bs(x=rul[g], t=z[k],  ff=Df[1], beta_h = beta_h_ed[1], beta_l = beta_l_ed[1], 
                      tau = tau_ed[1], lambda = lambda_ed, p=p_ed)
    rul3d_df_ed <- rbind(rul3d_df_ed, data.frame(RUL=rul[g], t=as.character(z[k]), pdf=dens, cdf=dis))
    
    dens <- rul_pdf_bs(x=rul[g], t=z[k], ff=Df[1], beta_h = beta_h_wie[1], beta_l = beta_l_wie[1], 
                       tau = tau_wie[1], lambda = lambda_wie, p=p_wie)
    dis <- rul_cdf_bs(x=rul[g], t=z[k],  ff=Df[1], beta_h = beta_h_wie[1], beta_l = beta_l_wie[1], 
                      tau = tau_wie[1], lambda = lambda_wie, p=p_wie)
    rul3d_df_wie <- rbind(rul3d_df_wie, data.frame(RUL=rul[g], t=as.character(z[k]), pdf=dens, cdf=dis))
    
    dens <- rul_pdf_bs(x=rul[g], t=z[k], ff=Df[1], beta_h = beta_h_gamma[1], beta_l = beta_l_gamma[1], 
                       tau = tau_gamma[1], lambda = lambda_gamma, p=p_gamma)
    dis <- rul_cdf_bs(x=rul[g], t=z[k],  ff=Df[1], beta_h = beta_h_gamma[1], beta_l = beta_l_gamma[1], 
                      tau = tau_gamma[1], lambda = lambda_gamma, p=p_gamma)
    rul3d_df_gamma <- rbind(rul3d_df_gamma, data.frame(RUL=rul[g], t=as.character(z[k]), pdf=dens, cdf=dis))
    
    dens <- rul_pdf_bs(x=rul[g], t=z[k], ff=Df[1], beta_h = beta_h_ig[1], beta_l = beta_l_ig[1], 
                       tau = tau_ig[1], lambda = lambda_ig, p=p_ig)
    dis <- rul_cdf_bs(x=rul[g], t=z[k],  ff=Df[1], beta_h = beta_h_ig[1], beta_l = beta_l_ig[1], 
                      tau = tau_ig[1], lambda = lambda_ig, p=p_ig)
    rul3d_df_ig <- rbind(rul3d_df_ig, data.frame(RUL=rul[g], t=as.character(z[k]), pdf=dens, cdf=dis))
    
    dens <- rul_pdf_bs(x=rul[g], t=z[k], ff=Df[1], beta_h = beta_h_true[1], beta_l = beta_l_true[1], 
                       tau = tau_true[1], lambda = lambda_true, p=p_true)
    dis <- rul_cdf_bs(x=rul[g], t=z[k],  ff=Df[1], beta_h = beta_h_true[1], beta_l = beta_l_true[1], 
                      tau = tau_ig[1], lambda = lambda_true, p=p_true)
    rul3d_df_true <- rbind(rul3d_df_true, data.frame(RUL=rul[g], t=as.character(z[k]), pdf=dens, cdf=dis))
  }
}

## MRL
mrl_wie <- mrl_wie_ig <- mrl_gamma <- mrl_ig <- mrl_ed <- mrl_true <- matrix(NA, nrow=I, ncol=length(z))


for(k in 1:I){
  for(g in 1:length(z)){
    if(z[g]<=tau_wie[k]){
      temp1 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-beta_h_wie[k]*z[g])/beta_h_wie[k], 
                                             shape=beta_h_wie[k]^((p_wie-1)/2)/sqrt(lambda_wie*(Df[k]--beta_h_wie[k]*z[g]) )), 
                         lower = 0.01, upper = tau_wie[k]-z[g])$value
      temp2 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-(beta_h_wie[k]-beta_l_wie[k])*tau_wie[k]-beta_l_wie[k]*z[g])/beta_l_wie[k], 
                                             shape=beta_l_wie[k]^((p_wie-1)/2)/sqrt(lambda_wie*((Df[k]-(beta_h_wie[k]-beta_l_wie[k])*tau_wie[k]-beta_l_wie[k]*z[g])))),
                         lower = tau_wie[k]-z[g], upper = 1000)$value
      mrl_wie[k,g] <- temp1+temp2}else{
        mrl_wie[k,g] <- beta_l_wie[k]^{(p_wie-2)}/2/lambda_wie+(Df[k]-(beta_h_wie[k]-beta_l_wie[k])*tau_wie[k]-beta_l_wie[k]*z[g])/beta_l_wie[k]
      }
    
    if(z[g]<=tau_gamma[k]){
      temp1 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-beta_h_gamma[k]*z[g])/beta_h_gamma[k], 
                                             shape=beta_h_gamma[k]^((p_gamma-1)/2)/sqrt(lambda_gamma*(Df[k]--beta_h_gamma[k]*z[g]) )), 
                         lower = 0.01, upper = tau_gamma[k]-z[g])$value
      temp2 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-(beta_h_gamma[k]-beta_l_gamma[k])*tau_gamma[k]-beta_l_gamma[k]*z[g])/beta_l_gamma[k], 
                                             shape=beta_l_gamma[k]^((p_gamma-1)/2)/sqrt(lambda_gamma*((Df[k]-(beta_h_gamma[k]-beta_l_gamma[k])*tau_gamma[k]-beta_l_gamma[k]*z[g])))),
                         lower = tau_gamma[k]-z[g], upper = 1000)$value
      mrl_gamma[k,g] <- temp1+temp2}else{
        mrl_gamma[k,g] <- beta_l_gamma[k]^{(p_gamma-2)}/2/lambda_gamma+(Df[k]-(beta_h_gamma[k]-beta_l_gamma[k])*tau_gamma[k]-beta_l_gamma[k]*z[g])/beta_l_gamma[k]
      }
    
    if(z[g]<=tau_ig[k]){
      temp1 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-beta_h_ig[k]*z[g])/beta_h_ig[k], 
                                             shape=beta_h_ig[k]^((p_ig-1)/2)/sqrt(lambda_ig*(Df[k]--beta_h_ig[k]*z[g]) )), 
                         lower = 0.01, upper = tau_ig[k]-z[g])$value
      temp2 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-(beta_h_ig[k]-beta_l_ig[k])*tau_ig[k]-beta_l_ig[k]*z[g])/beta_l_ig[k], 
                                             shape=beta_l_ig[k]^((p_ig-1)/2)/sqrt(lambda_ig*((Df[k]-(beta_h_ig[k]-beta_l_ig[k])*tau_ig[k]-beta_l_ig[k]*z[g])))),
                         lower = tau_ig[k]-z[g], upper = 1000)$value
      mrl_ig[k,g] <- temp1+temp2}else{
        mrl_ig[k,g] <- beta_l_ig[k]^{(p_ig-2)}/2/lambda_ig+(Df[k]-(beta_h_ig[k]-beta_l_ig[k])*tau_ig[k]-beta_l_ig[k]*z[g])/beta_l_ig[k]
      }
    
    if(z[g]<=tau_ed[k]){
      temp1 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-beta_h_ed[k]*z[g])/beta_h_ed[k], 
                                             shape=beta_h_ed[k]^((p_ed-1)/2)/sqrt(lambda_ed*(Df[k]--beta_h_ed[k]*z[g]) )), 
                         lower = 0.01, upper = tau_ed[k]-z[g])$value
      temp2 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k]-beta_l_ed[k]*z[g])/beta_l_ed[k], 
                                             shape=beta_l_ed[k]^((p_ed-1)/2)/sqrt(lambda_ed*((Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k]-beta_l_ed[k]*z[g])))),
                         lower = tau_ed[k]-z[g], upper = 1000)$value
      mrl_ed[k,g] <- temp1+temp2}else{
        mrl_ed[k,g] <- beta_l_ed[k]^{(p_ed-2)}/2/lambda_ed+(Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k]-beta_l_ed[k]*z[g])/beta_l_ed[k]
      }
    
    if(z[g]<=tau_true[k]){
      temp1 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-beta_h_true[k]*z[g])/beta_h_true[k], 
                                             shape=beta_h_true[k]^((p_true-1)/2)/sqrt(lambda_true*(Df[k]--beta_h_true[k]*z[g]) )), 
                         lower = 0.01, upper = tau_true[k]-z[g])$value
      temp2 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-(beta_h_true[k]-beta_l_true[k])*tau_true[k]-beta_l_true[k]*z[g])/beta_l_true[k], 
                                             shape=beta_l_true[k]^((p_true-1)/2)/sqrt(lambda_true*((Df[k]-(beta_h_true[k]-beta_l_true[k])*tau_true[k]-beta_l_true[k]*z[g])))),
                         lower = tau_true[k]-z[g], upper = 1000)$value
      mrl_true[k,g] <- temp1+temp2}else{
        mrl_true[k,g] <- beta_l_true[k]^{(p_true-2)}/2/lambda_true+(Df[k]-(beta_h_true[k]-beta_l_true[k])*tau_true[k]-beta_l_true[k]*z[g])/beta_l_true[k]
      }
    
  }
}

mrl1 <- round(cbind(mrl_wie[,c(3,5,9)], mrl_gamma[,c(3,5,9)], mrl_ig[,c(3,5,9)], mrl_ed[,c(3,5,9)], mrl_true[,c(3,5,9)]),3)

data_rul <- data.frame(pdf=rep(0,length(z)), t=z, RUL= mrl_true[1,])

rul3d_df <- cbind(rul3d_df_wie, rul3d_df_gamma, rul3d_df_ig, rul3d_df_ed, rul3d_df_true)
names(rul3d_df) <- c("RUL1", "t1", "pdf1", "cdf1", "RUL2", "t2", "pdf2", "cdf2", 
                     "RUL3", "t3", "pdf3", "cdf3", "RUL4", "t4", "pdf4", "cdf4", 
                     "RUL5", "t5", "pdf5", "cdf5")

# plot1_pdf <- plot_ly(rul3d_df, x = ~RUL1, y = ~t1, z = ~pdf1,  type = 'scatter3d', mode = 'lines', 
#                  color=~t1, line = list(color="blue", dash="dot", width = 3), showlegend = FALSE) %>% 
#   layout(scene = list(yaxis = list( title = "Current time" ), xaxis = list( title = "RUL" ), zaxis = list( title = "PDF" )))
# plot2_pdf <- plot1_pdf %>%add_trace(rul3d_df, x = ~RUL4, y = ~t4, z = ~pdf4, color=~t4, 
#                             line = list(color="green", dash="dash", width =3))
# plot3_pdf <- plot2_pdf %>%add_trace(rul3d_df, x = ~RUL5, y = ~t5, z = ~pdf5, color=~t5, 
#                             line = list(color="pink", dash="solid", width =3))
# plot_pdf <- plot3_pdf%>% add_markers(data=data_rul, x = ~RUL, y = ~t, z = ~pdf, color = I("black"), symbol=I(8), size=1,
#                           line = list(color="black", dash="dash", width =3, showlegend = FALSE)) 

plot1_pdf <- plot_ly(rul3d_df, x = ~RUL1, y = ~t1, z = ~pdf1,  type = 'scatter3d', mode = 'lines', 
                     color=~t1, line = list(color="pink", dash="solid", width = 3), showlegend = FALSE) %>% 
  layout(scene = list(yaxis = list( title = "Current time" ), xaxis = list( title = "RUL" ), zaxis = list( title = "PDF" )))
plot2_pdf <- plot1_pdf %>%add_trace(rul3d_df, x = ~RUL4, y = ~t4, z = ~pdf4, color=~t4, 
                                    line = list(color="green", dash="dash", width =3))
plot_pdf <- plot2_pdf%>% add_markers(data=data_rul, x = ~RUL, y = ~t, z = ~pdf, color = I("black"), symbol=I(8), size=1,
                                     line = list(color="black", dash="dash", width =3, showlegend = FALSE)) 

## rul_3d_wie

plot1_cdf <- plot_ly(rul3d_df, x = ~RUL1, y = ~t1, z = ~cdf1,  type = 'scatter3d', mode = 'lines', 
                 color=~t1, line = list(color="pink", dash="solid", width = 4), showlegend = FALSE) %>% 
  layout(scene = list(yaxis = list( title = "Current time" ), xaxis = list( title = "RUL" ), zaxis = list( title = "CDF" )))
plot2_cdf <- plot1_cdf %>%add_trace(rul3d_df, x = ~RUL4, y = ~t4, z = ~cdf4, color=~t4, 
                            line = list(color="green", dash="dash", width =4))
plot_cdf <- plot2_cdf%>% add_markers(data=data_rul, x = ~RUL, y = ~t, z = ~pdf, color = I("black"), symbol=I(8), size=1,
                                     line = list(color="black", dash="dash", width =3, showlegend = FALSE)) 

data_plot=cbind(rul3d_df_ed[rul3d_df_ed$t%in%c(10,20,40),], model="CPTEDP")
data_plot=rbind(data_plot, cbind(rul3d_df_wie[rul3d_df_wie$t%in%c(10,20,40),], model="CPWP"))
data_plot=rbind(data_plot, cbind(rul3d_df_gamma[rul3d_df_gamma$t%in%c(10,20,40),], model="CPGP"))
data_plot=rbind(data_plot, cbind(rul3d_df_ig[rul3d_df_ig$t%in%c(10,20,40),], model="CPIGP"))



rul_pdf_ed <- ggplot(data=data_plot, aes(x=RUL, y= pdf, linetype = model, shape = model, color=t))+ xlim(0,80)+ylim(0,0.6)+
  geom_line( size=0.3)+geom_point(size=0.2)+xlab("RUL")+ylab("PDF")+ theme_bw()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size = 8),
        legend.key.size=unit(0.3,'cm'),legend.text = element_text(size = 7),
        legend.title=element_blank())+annotate("text", x=15, y=0.6, label="t=40")+
  annotate("text", x=35, y=0.38, label="t=20")+
  annotate("text", x=45, y=0.34, label="t=10")


rul_cdf_ed <- ggplot(data=data_plot, aes(x=RUL, y= cdf, linetype = model, shape = model, color=t))+ xlim(0,80)+ylim(0,1.05)+
  geom_line(size=0.3)+geom_point(size=0.2)+xlab("RUL")+ylab("CDF")+ theme_bw()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size = 8),
        legend.key.size=unit(0.3,'cm'),legend.text = element_text(size = 7),
        legend.title=element_blank())+annotate("text", x=20, y=1.02, label="t=40")+
  annotate("text", x=37, y=1.02, label="t=20")+
  annotate("text", x=50, y=1.02, label="t=10")

postscript("D:/My Document/2023/paper/ED_change_point/figure/simulation/rul_pdf_wie.eps", width=4, height=4)
rul_pdf_ed
dev.off()

postscript("D:/My Document/2023/paper/ED_change_point/figure/simulation/rul_cdf_wie.eps", width=4, height=4)
rul_cdf_ed
dev.off()

## mttf
mttf_wie <- mttf_wie_ig <-mttf_gamma <- mttf_ig <-mttf_ed <- numeric(I)

for(k in 1:I){
  temp1 <- integrate(function(t) t*dbisa(t, scale=Df[k]/beta_h_wie[k], shape=beta_h_wie[k]^((p_wie-1)/2)/sqrt(lambda_wie*Df[k])), 
                     lower = 0, upper = tau_wie[k])$value
  temp2 <- integrate(function(t) t*dbisa(t, scale=(Df[k]-(beta_h_wie[k]-beta_l_wie[k])*tau_wie[k])/beta_l_wie[k], 
                                         shape=beta_l_wie[k]^((p_wie-1)/2)/sqrt(lambda_wie*((Df[k]-(beta_h_wie[k]-beta_l_wie[k])*tau_wie[k])))),
                     lower = tau_wie[k], upper = 1000)$value
  mttf_wie[k] <- temp1+temp2
  
  temp1 <- integrate(function(t) t*dbisa(t, scale=Df[k]/beta_h_gamma[k], shape=beta_h_gamma[k]^((p_gamma-1)/2)/sqrt(lambda_gamma*Df[k])), 
                     lower = 0, upper = tau_gamma[k])$value
  temp2 <- integrate(function(t) t*dbisa(t, scale=(Df[k]-(beta_h_gamma[k]-beta_l_gamma[k])*tau_gamma[k])/beta_l_gamma[k], 
                                         shape=beta_l_gamma[k]^((p_gamma-1)/2)/sqrt(lambda_gamma*((Df[k]-(beta_h_gamma[k]-beta_l_gamma[k])*tau_gamma[k])))),
                     lower = tau_gamma[k], upper = 1000)$value
  mttf_gamma[k] <- temp1+temp2
  
  temp1 <- integrate(function(t) t*dbisa(t, scale=Df[k]/beta_h_ig[k], shape=beta_h_ig[k]^((p_ig-1)/2)/sqrt(lambda_ig*Df[k])), 
                     lower = 0, upper = tau_ig[k])$value
  temp2 <- integrate(function(t) t*dbisa(t, scale=(Df[k]-(beta_h_ig[k]-beta_l_ig[k])*tau_ig[k])/beta_l_ig[k], 
                                         shape=beta_l_ig[k]^((p_ig-1)/2)/sqrt(lambda_ig*((Df[k]-(beta_h_ig[k]-beta_l_ig[k])*tau_ig[k])))),
                     lower = tau_ig[k], upper = 1000)$value
  mttf_ig[k] <- temp1+temp2
  
  temp1 <- integrate(function(t) t*dbisa(t, scale=Df[k]/beta_h_ed[k], shape=beta_h_ed[k]^((p_ed-1)/2)/sqrt(lambda_ed*Df[k])), 
                     lower = 0, upper = tau_ed[k])$value
  temp2 <- integrate(function(t) t*dbisa(t, scale=(Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k])/beta_l_ed[k], 
                                         shape=beta_l_ed[k]^((p_ed-1)/2)/sqrt(lambda_ed*((Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k])))),
                     lower = tau_ed[k], upper = 1000)$value
  mttf_ed[k] <- temp1+temp2
}

for(i in 1:I){
  temp1 <- Df[i]/beta_h_wie[i]
  temp2 <- 1-pinvgauss(q=Df[i]^2/(tau_wie[i]*beta_h_wie[i]^2), mean=Df[i]/beta_h_wie[i], 
                       shape=Df[i]^2*lambda_wie)
  temp3 <- (Df[i]-(beta_h_wie[i]-beta_l_wie[i])*tau_wie[i])/beta_l_wie[i]
  temp4 <- pinvgauss(q=(Df[i]-(beta_h_wie[i]-beta_l_wie[i])*tau_wie[i])^2/(tau_wie[i]*beta_l_wie[i]^2), 
                     mean=(Df[i]-(beta_h_wie[i]-beta_l_wie[i])*tau_wie[i])/beta_l_wie[i], 
                     shape=(Df[i]-(beta_h_wie[i]-beta_l_wie[i])*tau_wie[i])^2*lambda_wie)
  mttf_wie_ig[i] <- temp1*temp2+temp3*temp4
}

mttf1 <- rbind(mttf_wie, mttf_wie_ig, mttf_gamma, mttf_ig, mttf_ed)






############################# Gamma process############################################################################ 
I=10
Df <- rep(200, I)
g=1
est_wie <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point_simulation/simulation/gamma/","bugs_wie", 1, "/", "log.txt", sep=""))[[1]]
est_gamma <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point_simulation/simulation/gamma/","bugs_gamma", 1, "/", "log.txt", sep=""))[[1]]

est_ig <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point_simulation/simulation/gamma/","bugs_ig", 1, "/", "log.txt", sep=""))[[1]]
est_ed <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point_simulation/simulation/gamma/","bugs_ed", 1, "/", "log.txt", sep=""))[[1]]

beta_h_wie <- est_wie[(0:9)*3+1,1] 
beta_l_wie <- est_wie[(0:9)*3+2,1] 
tau_wie <- est_wie[(0:9)*3+3,1] 
lambda_wie <- est_wie[31,1] 
p_wie=2

beta_h_gamma <- est_gamma[(0:9)*3+1,1] 
beta_l_gamma <- est_gamma[(0:9)*3+2,1] 
tau_gamma <- est_gamma[(0:9)*3+3,1] 
lambda_gamma <- est_gamma[31,1] 
p_gamma <- 2 

beta_h_ig <- est_ig[(0:9)*3+1,1] 
beta_l_ig <- est_ig[(0:9)*3+2,1] 
tau_ig <- est_ig[(0:9)*3+3,1] 
lambda_ig <- est_ig[31,1] 
p_ig=3

beta_h_ed <- est_ed[(0:9)*3+1,1] 
beta_l_ed <- est_ed[(0:9)*3+2,1] 
tau_ed <- est_ed[(0:9)*3+3,1] 
lambda_ed <- est_ed[31,1] 
p_ed <- est_ed[35,1] 



## Fail-time distribution 
time <- seq(20, 120, length =1000)

fpdf_res_wie <-  fpdf_res_gamma <-fpdf_res_ig <-  fpdf_res_ed <- matrix(NA, I, length(time))
for(k in 1:I){
  for(i in 1:length(time)){
    fpdf_res_wie[k,i] <- ft_pdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_wie[k], beta_l = beta_l_wie[k], 
                               tau = tau_wie[k], lambda = lambda_wie, p=p_wie)
    fpdf_res_gamma[k,i] <- ft_pdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_gamma[k], beta_l = beta_l_gamma[k], 
                                 tau = tau_gamma[k], lambda = lambda_gamma, p=p_gamma)
    fpdf_res_ig[k,i] <- ft_pdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_ig[k], beta_l = beta_l_ig[k], 
                              tau = tau_ig[k], lambda = lambda_ig, p=p_ig)
    fpdf_res_ed[k,i] <- ft_pdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_ed[k], beta_l = beta_l_ed[k], 
                              tau = tau_ed[k], lambda = lambda_ed, p=p_ed)
  }
}

fcdf_res_wie <-  fcdf_res_gamma <-fcdf_res_ig <-  fcdf_res_ed <- matrix(NA, I, length(time))
for(k in 1:I){
  for(i in 1:length(time)){
    fcdf_res_wie[k,i] <- ft_cdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_wie[k], beta_l = beta_l_wie[k], 
                                   tau = tau_wie[k], lambda = lambda_wie, p=p_wie)
    fcdf_res_gamma[k,i] <- ft_cdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_gamma[k], beta_l = beta_l_gamma[k], 
                                     tau = tau_gamma[k], lambda = lambda_gamma, p=p_gamma)
    fcdf_res_ig[k,i] <- ft_cdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_ig[k], beta_l = beta_l_ig[k], 
                                  tau = tau_ig[k], lambda = lambda_ig, p=p_ig)
    fcdf_res_ed[k,i] <- ft_cdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_ed[k], beta_l = beta_l_ed[k], 
                                  tau = tau_ed[k], lambda = lambda_ed, p=p_ed)
  }
}


faildata <- data.frame()
for(k in 1:9){
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_wie[k,], ft_cdf=fcdf_res_wie[k,], GaAs=paste("Unit # 0", k, sep = ""), model="CPWP"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_gamma[k,], ft_cdf=fcdf_res_gamma[k,], GaAs=paste("Unit # 0", k, sep = ""), model="CPGP"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_ig[k,], ft_cdf=fcdf_res_ig[k,], GaAs=paste("Unit # 0", k, sep = ""), model="CPIGP"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_ed[k,], ft_cdf=fcdf_res_ed[k,], GaAs=paste("Unit # 0", k, sep = ""), model="CPTEDP"))
}

for(k in 10:I){
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_wie[k,], ft_cdf=fcdf_res_wie[k,], GaAs=paste("Unit #", k), model="CPWP"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_gamma[k,], ft_cdf=fcdf_res_gamma[k,], GaAs=paste("Unit #", k), model="CPGP"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_ig[k,], ft_cdf=fcdf_res_ig[k,], GaAs=paste("Unit #", k), model="CPIGP"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_ed[k,], ft_cdf=fcdf_res_ed[k,], GaAs=paste("Unit #", k), model="CPTEDP"))
}




#fail_pdf_ed <- 


ft_pdf_plot <- ggplot(data=faildata, aes(x=ft, y= ft_pdf, linetype=model, color=model))+facet_wrap(~GaAs, nrow=2)+
  geom_line(show.legend = T, size=1)+xlab("Failure-time")+ylab("PDF")+ theme_bw()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size = 8),
        legend.key.size=unit(0.3,'cm'),legend.text = element_text(size = 7),
        legend.title=element_blank())

ft_cdf_plot <- ggplot(data=faildata, aes(x=ft, y= ft_cdf, linetype=model, color=model))+facet_wrap(~GaAs, nrow=2)+
  geom_line(show.legend = T, size=1)+xlab("Failure-time")+ylab("CDF")+ theme_bw()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size = 8),
        legend.key.size=unit(0.3,'cm'),legend.text = element_text(size = 7),
        legend.title=element_blank())

postscript("D:/My Document/2023/paper/ED_change_point/figure/simulation/fail_pdf_gamma.eps", width=10, height=5)
ft_pdf_plot
dev.off()

postscript("D:/My Document/2023/paper/ED_change_point/figure/simulation/fail_cdf_gamma.eps", width=10, height=5)
ft_cdf_plot
dev.off()

## RUL
i=1
z <- seq(0, 40, by = 5)
rul <- seq(0, 100, length=300)
p_true= 2
lambda_true=1

rul_pdf_bs <- function(x, t, ff, beta_h, beta_l, tau, lambda, p){
  if(t<tau){
    if(t+x<tau){
      ft <- dbisa(x, scale=(ff-beta_h*t)/beta_h, shape=beta_h^((p-1)/2)/sqrt(lambda*(ff-beta_h*t)))
    } else {
      ft <- dbisa(x, scale=(ff-(beta_h-beta_l)*tau-beta_l*t)/beta_l, 
                  shape=beta_l^((p-1)/2)/sqrt(lambda*((ff-(beta_h-beta_l)*tau-beta_l*t))))
    }
  }else{
    ft <- dbisa(x, scale=(ff-(beta_h-beta_l)*tau-beta_l*t)/beta_l, 
                shape=beta_l^((p-1)/2)/sqrt(lambda*((ff-(beta_h-beta_l)*tau-beta_l*t))))
  }
  return(ft)
}

rul_cdf_bs <- function(x, t, ff, beta_h, beta_l, tau, lambda, p){
  if(t<tau){
    if(t+x<tau){
      ft <- pbisa(x, scale=(ff-beta_h*t)/beta_h, shape=beta_h^((p-1)/2)/sqrt(lambda*(ff-beta_h*t)))
    } else {
      ft <- pbisa(x, scale=(ff-(beta_h-beta_l)*tau-beta_l*t)/beta_l, 
                  shape=beta_l^((p-1)/2)/sqrt(lambda*((ff-(beta_h-beta_l)*tau-beta_l*t))))
    }
  }else{
    ft <- pbisa(x, scale=(ff-(beta_h-beta_l)*tau-beta_l*t)/beta_l, 
                shape=beta_l^((p-1)/2)/sqrt(lambda*((ff-(beta_h-beta_l)*tau-beta_l*t))))
  }
  return(ft)
}

rul3d_df_ed <- rul3d_df_wie <- rul3d_df_ig <- rul3d_df_gamma <-  rul3d_df_true <- data.frame()

for(k in 1:length(z)){
  for(g in 1:length(rul)){
    dens <- rul_pdf_bs(x=rul[g], t=z[k], ff=Df[1], beta_h = beta_h_ed[1], beta_l = beta_l_ed[1], 
                       tau = tau_ed[1], lambda = lambda_ed, p=p_ed)
    dis <- rul_cdf_bs(x=rul[g], t=z[k],  ff=Df[1], beta_h = beta_h_ed[1], beta_l = beta_l_ed[1], 
                      tau = tau_ed[1], lambda = lambda_ed, p=p_ed)
    rul3d_df_ed <- rbind(rul3d_df_ed, data.frame(RUL=rul[g], t=as.character(z[k]), pdf=dens, cdf=dis))
    
    dens <- rul_pdf_bs(x=rul[g], t=z[k], ff=Df[1], beta_h = beta_h_wie[1], beta_l = beta_l_wie[1], 
                       tau = tau_wie[1], lambda = lambda_wie, p=p_wie)
    dis <- rul_cdf_bs(x=rul[g], t=z[k],  ff=Df[1], beta_h = beta_h_wie[1], beta_l = beta_l_wie[1], 
                      tau = tau_wie[1], lambda = lambda_wie, p=p_wie)
    rul3d_df_wie <- rbind(rul3d_df_wie, data.frame(RUL=rul[g], t=as.character(z[k]), pdf=dens, cdf=dis))
    
    dens <- rul_pdf_bs(x=rul[g], t=z[k], ff=Df[1], beta_h = beta_h_gamma[1], beta_l = beta_l_gamma[1], 
                       tau = tau_gamma[1], lambda = lambda_gamma, p=p_gamma)
    dis <- rul_cdf_bs(x=rul[g], t=z[k],  ff=Df[1], beta_h = beta_h_gamma[1], beta_l = beta_l_gamma[1], 
                      tau = tau_gamma[1], lambda = lambda_gamma, p=p_gamma)
    rul3d_df_gamma <- rbind(rul3d_df_gamma, data.frame(RUL=rul[g], t=as.character(z[k]), pdf=dens, cdf=dis))
    
    dens <- rul_pdf_bs(x=rul[g], t=z[k], ff=Df[1], beta_h = beta_h_ig[1], beta_l = beta_l_ig[1], 
                       tau = tau_ig[1], lambda = lambda_ig, p=p_ig)
    dis <- rul_cdf_bs(x=rul[g], t=z[k],  ff=Df[1], beta_h = beta_h_ig[1], beta_l = beta_l_ig[1], 
                      tau = tau_ig[1], lambda = lambda_ig, p=p_ig)
    rul3d_df_ig <- rbind(rul3d_df_ig, data.frame(RUL=rul[g], t=as.character(z[k]), pdf=dens, cdf=dis))
    
    dens <- rul_pdf_bs(x=rul[g], t=z[k], ff=Df[1], beta_h = beta_h_true[1], beta_l = beta_l_true[1], 
                       tau = tau_true[1], lambda = lambda_true, p=p_true)
    dis <- rul_cdf_bs(x=rul[g], t=z[k],  ff=Df[1], beta_h = beta_h_true[1], beta_l = beta_l_true[1], 
                      tau = tau_ig[1], lambda = lambda_true, p=p_true)
    rul3d_df_true <- rbind(rul3d_df_true, data.frame(RUL=rul[g], t=as.character(z[k]), pdf=dens, cdf=dis))
  }
}

## MRL

mrl_wie <- mrl_wie_ig <- mrl_gamma <- mrl_ig <- mrl_ed <- mrl_true <- matrix(NA, nrow=I, ncol=length(z))


for(k in 1:I){
  for(g in 1:length(z)){
    if(z[g]<=tau_wie[k]){
      temp1 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-beta_h_wie[k]*z[g])/beta_h_wie[k], 
                                             shape=beta_h_wie[k]^((p_wie-1)/2)/sqrt(lambda_wie*(Df[k]--beta_h_wie[k]*z[g]) )), 
                         lower = 0, upper = tau_wie[k]-z[g])$value
      temp2 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-(beta_h_wie[k]-beta_l_wie[k])*tau_wie[k]-beta_l_wie[k]*z[g])/beta_l_wie[k], 
                                             shape=beta_l_wie[k]^((p_wie-1)/2)/sqrt(lambda_wie*((Df[k]-(beta_h_wie[k]-beta_l_wie[k])*tau_wie[k]-beta_l_wie[k]*z[g])))),
                         lower = tau_wie[k]-z[g], upper = 1000)$value
      mrl_wie[k,g] <- temp1+temp2}else{
        mrl_wie[k,g] <- beta_l_wie[k]^{(p_wie-2)}/2/lambda_wie+(Df[k]-(beta_h_wie[k]-beta_l_wie[k])*tau_wie[k]-beta_l_wie[k]*z[g])/beta_l_wie[k]
      }
    
    if(z[g]<=tau_gamma[k]){
      temp1 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-beta_h_gamma[k]*z[g])/beta_h_gamma[k], 
                                             shape=beta_h_gamma[k]^((p_gamma-1)/2)/sqrt(lambda_gamma*(Df[k]--beta_h_gamma[k]*z[g]) )), 
                         lower = 0, upper = tau_gamma[k]-z[g])$value
      temp2 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-(beta_h_gamma[k]-beta_l_gamma[k])*tau_gamma[k]-beta_l_gamma[k]*z[g])/beta_l_gamma[k], 
                                             shape=beta_l_gamma[k]^((p_gamma-1)/2)/sqrt(lambda_gamma*((Df[k]-(beta_h_gamma[k]-beta_l_gamma[k])*tau_gamma[k]-beta_l_gamma[k]*z[g])))),
                         lower = tau_gamma[k]-z[g], upper = 1000)$value
      mrl_gamma[k,g] <- temp1+temp2}else{
        mrl_gamma[k,g] <- beta_l_gamma[k]^{(p_gamma-2)}/2/lambda_gamma+(Df[k]-(beta_h_gamma[k]-beta_l_gamma[k])*tau_gamma[k]-beta_l_gamma[k]*z[g])/beta_l_gamma[k]
      }
    
    if(z[g]<=tau_ig[k]){
      temp1 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-beta_h_ig[k]*z[g])/beta_h_ig[k], 
                                             shape=beta_h_ig[k]^((p_ig-1)/2)/sqrt(lambda_ig*(Df[k]--beta_h_ig[k]*z[g]) )), 
                         lower = 0, upper = tau_ig[k]-z[g])$value
      temp2 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-(beta_h_ig[k]-beta_l_ig[k])*tau_ig[k]-beta_l_ig[k]*z[g])/beta_l_ig[k], 
                                             shape=beta_l_ig[k]^((p_ig-1)/2)/sqrt(lambda_ig*((Df[k]-(beta_h_ig[k]-beta_l_ig[k])*tau_ig[k]-beta_l_ig[k]*z[g])))),
                         lower = tau_ig[k]-z[g], upper = 1000)$value
      mrl_ig[k,g] <- temp1+temp2}else{
        mrl_ig[k,g] <- beta_l_ig[k]^{(p_ig-2)}/2/lambda_ig+(Df[k]-(beta_h_ig[k]-beta_l_ig[k])*tau_ig[k]-beta_l_ig[k]*z[g])/beta_l_ig[k]
      }
    
    if(z[g]<=tau_ed[k]){
      temp1 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-beta_h_ed[k]*z[g])/beta_h_ed[k], 
                                             shape=beta_h_ed[k]^((p_ed-1)/2)/sqrt(lambda_ed*(Df[k]--beta_h_ed[k]*z[g]) )), 
                         lower = 0, upper = tau_ed[k]-z[g])$value
      temp2 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k]-beta_l_ed[k]*z[g])/beta_l_ed[k], 
                                             shape=beta_l_ed[k]^((p_ed-1)/2)/sqrt(lambda_ed*((Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k]-beta_l_ed[k]*z[g])))),
                         lower = tau_ed[k]-z[g], upper = 1000)$value
      mrl_ed[k,g] <- temp1+temp2}else{
        mrl_ed[k,g] <- beta_l_ed[k]^{(p_ed-2)}/2/lambda_ed+(Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k]-beta_l_ed[k]*z[g])/beta_l_ed[k]
      }
    
    if(z[g]<=tau_true[k]){
      temp1 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-beta_h_true[k]*z[g])/beta_h_true[k], 
                                             shape=beta_h_true[k]^((p_true-1)/2)/sqrt(lambda_true*(Df[k]--beta_h_true[k]*z[g]) )), 
                         lower = 0, upper = tau_true[k]-z[g])$value
      temp2 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-(beta_h_true[k]-beta_l_true[k])*tau_true[k]-beta_l_true[k]*z[g])/beta_l_true[k], 
                                             shape=beta_l_true[k]^((p_true-1)/2)/sqrt(lambda_true*((Df[k]-(beta_h_true[k]-beta_l_true[k])*tau_true[k]-beta_l_true[k]*z[g])))),
                         lower = tau_true[k]-z[g], upper = 1000)$value
      mrl_true[k,g] <- temp1+temp2}else{
        mrl_true[k,g] <- beta_l_true[k]^{(p_true-2)}/2/lambda_true+(Df[k]-(beta_h_true[k]-beta_l_true[k])*tau_true[k]-beta_l_true[k]*z[g])/beta_l_true[k]
      }
    
  }
}

mrl2 <- round(cbind(mrl_wie[,c(3,5,9)], mrl_gamma[,c(3,5,9)], mrl_ig[,c(3,5,9)], mrl_ed[,c(3,5,9)], mrl_true[,c(3,5,9)]),3)

data_rul <- data.frame(pdf=rep(0,length(z)), t=z, RUL= mrl_true[1,])


rul3d_df <- cbind(rul3d_df_wie, rul3d_df_gamma, rul3d_df_ig, rul3d_df_ed, rul3d_df_true)
names(rul3d_df) <- c("RUL1", "t1", "pdf1", "cdf1", "RUL2", "t2", "pdf2", "cdf2", 
                     "RUL3", "t3", "pdf3", "cdf3", "RUL4", "t4", "pdf4", "cdf4", 
                     "RUL5", "t5", "pdf5", "cdf5")

plot1_pdf <- plot_ly(rul3d_df, x = ~RUL2, y = ~t2, z = ~pdf2,  type = 'scatter3d', mode = 'lines', 
                     color=~t2, line = list(color="pink", dash="solid", width = 3), showlegend = FALSE) %>% 
  layout(scene = list(yaxis = list( title = "Current time" ), xaxis = list( title = "RUL" ), zaxis = list( title = "PDF" )))
plot2_pdf <- plot1_pdf %>%add_trace(rul3d_df, x = ~RUL4, y = ~t4, z = ~pdf4, color=~t4, 
                                line = list(color="green", dash="dash", width =3))
plot_pdf <- plot2_pdf%>% add_markers(data=data_rul, x = ~RUL, y = ~t, z = ~pdf, color = I("black"), symbol=I(8), size=1,
                                     line = list(color="black", dash="dash", width =3, showlegend = FALSE)) 

plot1_cdf <- plot_ly(rul3d_df, x = ~RUL2, y = ~t2, z = ~cdf2,  type = 'scatter3d', mode = 'lines', 
                     color=~t2, line = list(color="pink", dash="solid", width = 4), showlegend = FALSE) %>% 
  layout(scene = list(yaxis = list( title = "Current time" ), xaxis = list( title = "RUL" ), zaxis = list( title = "CDF" )))
plot2_cdf <- plot1_cdf %>%add_trace(rul3d_df, x = ~RUL4, y = ~t4, z = ~cdf4, color=~t4, 
                                    line = list(color="green", dash="dash", width =4))
plot_cdf <- plot2_cdf%>% add_markers(data=data_rul, x = ~RUL, y = ~t, z = ~pdf, color = I("black"), symbol=I(8), size=1,
                                     line = list(color="black", dash="dash", width =3, showlegend = FALSE)) 

data_plot=cbind(rul3d_df_ed[rul3d_df_ed$t%in%c(10,20,40),], model="CPTEDP")
data_plot=rbind(data_plot, cbind(rul3d_df_wie[rul3d_df_wie$t%in%c(10,20,40),], model="CPWP"))
data_plot=rbind(data_plot, cbind(rul3d_df_gamma[rul3d_df_gamma$t%in%c(10,20,40),], model="CPGP"))
data_plot=rbind(data_plot, cbind(rul3d_df_ig[rul3d_df_ig$t%in%c(10,20,40),], model="CPIGP"))
data_plot=data_plot[!(data_plot$t==10&data_plot$model=="CPWP"&data_plot$RUL<4.69),]

rul_pdf_ed <- ggplot(data=data_plot, aes(x=RUL, y= pdf, linetype = model, shape = model, color=t))+ xlim(0,100)+ylim(0,0.09)+
  geom_line( size=0.3)+geom_point(size=0.2)+xlab("RUL")+ylab("PDF")+ theme_bw()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size = 8),
        legend.key.size=unit(0.3,'cm'),legend.text = element_text(size = 7),
        legend.title=element_blank())+annotate("text", x=20, y=0.09, label="t=40")+
  annotate("text", x=40, y=0.065, label="t=20")+
  annotate("text", x=55, y=0.055, label="t=10")

rul_cdf_ed <- ggplot(data=data_plot, aes(x=RUL, y= cdf, linetype = model, shape = model, color=t))+ xlim(0,100)+ylim(0,1.05)+
  geom_line(size=0.3)+geom_point(size=0.2)+xlab("RUL")+ylab("CDF")+ theme_bw()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size = 8),
        legend.key.size=unit(0.3,'cm'),legend.text = element_text(size = 7),
        legend.title=element_blank())+annotate("text", x=20, y=0.85, label="t=40")+
  annotate("text", x=45, y=0.85, label="t=20")+
  annotate("text", x=55, y=0.85, label="t=10")

postscript("D:/My Document/2023/paper/ED_change_point/figure/simulation/rul_pdf_gamma.eps", width=4, height=4)
rul_pdf_ed
dev.off()

postscript("D:/My Document/2023/paper/ED_change_point/figure/simulation/rul_cdf_gamma.eps", width=4, height=4)
rul_cdf_ed
dev.off()


## mttf
mttf_wie <- mttf_wie_ig <-mttf_gamma <- mttf_ig <-mttf_ed <- numeric(I)

for(k in 1:I){
  temp1 <- integrate(function(t) t*dbisa(t, scale=Df[k]/beta_h_wie[k], shape=beta_h_wie[k]^((p_wie-1)/2)/sqrt(lambda_wie*Df[k])), 
                     lower = 0, upper = tau_wie[k])$value
  temp2 <- integrate(function(t) t*dbisa(t, scale=(Df[k]-(beta_h_wie[k]-beta_l_wie[k])*tau_wie[k])/beta_l_wie[k], 
                                         shape=beta_l_wie[k]^((p_wie-1)/2)/sqrt(lambda_wie*((Df[k]-(beta_h_wie[k]-beta_l_wie[k])*tau_wie[k])))),
                     lower = tau_wie[k], upper = 1000)$value
  mttf_wie[k] <- temp1+temp2
  
  temp1 <- integrate(function(t) t*dbisa(t, scale=Df[k]/beta_h_gamma[k], shape=beta_h_gamma[k]^((p_gamma-1)/2)/sqrt(lambda_gamma*Df[k])), 
                     lower = 0, upper = tau_gamma[k])$value
  temp2 <- integrate(function(t) t*dbisa(t, scale=(Df[k]-(beta_h_gamma[k]-beta_l_gamma[k])*tau_gamma[k])/beta_l_gamma[k], 
                                         shape=beta_l_gamma[k]^((p_gamma-1)/2)/sqrt(lambda_gamma*((Df[k]-(beta_h_gamma[k]-beta_l_gamma[k])*tau_gamma[k])))),
                     lower = tau_gamma[k], upper = 1000)$value
  mttf_gamma[k] <- temp1+temp2
  
  temp1 <- integrate(function(t) t*dbisa(t, scale=Df[k]/beta_h_ig[k], shape=beta_h_ig[k]^((p_ig-1)/2)/sqrt(lambda_ig*Df[k])), 
                     lower = 0, upper = tau_ig[k])$value
  temp2 <- integrate(function(t) t*dbisa(t, scale=(Df[k]-(beta_h_ig[k]-beta_l_ig[k])*tau_ig[k])/beta_l_ig[k], 
                                         shape=beta_l_ig[k]^((p_ig-1)/2)/sqrt(lambda_ig*((Df[k]-(beta_h_ig[k]-beta_l_ig[k])*tau_ig[k])))),
                     lower = tau_ig[k], upper = 1000)$value
  mttf_ig[k] <- temp1+temp2
  
  temp1 <- integrate(function(t) t*dbisa(t, scale=Df[k]/beta_h_ed[k], shape=beta_h_ed[k]^((p_ed-1)/2)/sqrt(lambda_ed*Df[k])), 
                     lower = 0, upper = tau_ed[k])$value
  temp2 <- integrate(function(t) t*dbisa(t, scale=(Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k])/beta_l_ed[k], 
                                         shape=beta_l_ed[k]^((p_ed-1)/2)/sqrt(lambda_ed*((Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k])))),
                     lower = tau_ed[k], upper = 1000)$value
  mttf_ed[k] <- temp1+temp2
}

for(i in 1:I){
  temp1 <- Df[i]/beta_h_wie[i]
  temp2 <- 1-pinvgauss(q=Df[i]^2/(tau_wie[i]*beta_h_wie[i]^2), mean=Df[i]/beta_h_wie[i], 
                       shape=Df[i]^2*lambda_wie)
  temp3 <- (Df[i]-(beta_h_wie[i]-beta_l_wie[i])*tau_wie[i])/beta_l_wie[i]
  temp4 <- pinvgauss(q=(Df[i]-(beta_h_wie[i]-beta_l_wie[i])*tau_wie[i])^2/(tau_wie[i]*beta_l_wie[i]^2), 
                     mean=(Df[i]-(beta_h_wie[i]-beta_l_wie[i])*tau_wie[i])/beta_l_wie[i], 
                     shape=(Df[i]-(beta_h_wie[i]-beta_l_wie[i])*tau_wie[i])^2*lambda_wie)
  mttf_wie_ig[i] <- temp1*temp2+temp3*temp4
}

mttf2 <- rbind(mttf_wie, mttf_wie_ig, mttf_gamma, mttf_ig, mttf_ed)





######################################### IG process #######################################
## Failuretime and RUL
I=10
Df <- rep(200, I)
g=1
est_wie <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point_simulation/simulation/ig/","bugs_wie",g, "/", "log.txt", sep=""))[[1]]
est_gamma <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point_simulation/simulation/ig/","bugs_gamma",g, "/", "log.txt", sep=""))[[1]]

est_ig <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point_simulation/simulation/ig/","bugs_ig",g, "/", "log.txt", sep=""))[[1]]
est_ed <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point_simulation/simulation/ig/","bugs_ed",g, "/", "log.txt", sep=""))[[1]]

beta_h_wie <- est_wie[(0:9)*3+1,1] 
beta_l_wie <- est_wie[(0:9)*3+2,1] 
tau_wie <- est_wie[(0:9)*3+3,1] 
lambda_wie <- est_wie[31,1] 
p_wie=0

beta_h_gamma <- est_gamma[(0:9)*3+1,1] 
beta_l_gamma <- est_gamma[(0:9)*3+2,1] 
tau_gamma <- est_gamma[(0:9)*3+3,1] 
lambda_gamma <- est_gamma[31,1] 
p_gamma <- 2 

beta_h_ig <- est_ig[(0:9)*3+1,1] 
beta_l_ig <- est_ig[(0:9)*3+2,1] 
tau_ig <- est_ig[(0:9)*3+3,1] 
lambda_ig <- est_ig[31,1] 
p_ig=3

beta_h_ed <- est_ed[(0:9)*3+1,1] 
beta_l_ed <- est_ed[(0:9)*3+2,1] 
tau_ed <- est_ed[(0:9)*3+3,1] 
lambda_ed <- est_ed[31,1] 
p_ed <- est_ed[35,1] 





## Fail-time distribution 
time <- seq(20, 120, length =1000)

fpdf_res_wie <-  fpdf_res_gamma <-fpdf_res_ig <-  fpdf_res_ed <- matrix(NA, I, length(time))
for(k in 1:I){
  for(i in 1:length(time)){
    fpdf_res_wie[k,i] <- ft_pdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_wie[k], beta_l = beta_l_wie[k], 
                               tau = tau_wie[k], lambda = lambda_wie, p=p_wie)
    fpdf_res_gamma[k,i] <- ft_pdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_gamma[k], beta_l = beta_l_gamma[k], 
                                 tau = tau_gamma[k], lambda = lambda_gamma, p=p_gamma)
    fpdf_res_ig[k,i] <- ft_pdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_ig[k], beta_l = beta_l_ig[k], 
                              tau = tau_ig[k], lambda = lambda_ig, p=p_ig)
    fpdf_res_ed[k,i] <- ft_pdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_ed[k], beta_l = beta_l_ed[k], 
                              tau = tau_ed[k], lambda = lambda_ed, p=p_ed)
  }
}


fcdf_res_wie <-  fcdf_res_gamma <-fcdf_res_ig <-  fcdf_res_ed <- matrix(NA, I, length(time))
for(k in 1:I){
  for(i in 1:length(time)){
    fcdf_res_wie[k,i] <- ft_cdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_wie[k], beta_l = beta_l_wie[k], 
                               tau = tau_wie[k], lambda = lambda_wie, p=p_wie)
    fcdf_res_gamma[k,i] <- ft_cdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_gamma[k], beta_l = beta_l_gamma[k], 
                                 tau = tau_gamma[k], lambda = lambda_gamma, p=p_gamma)
    fcdf_res_ig[k,i] <- ft_cdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_ig[k], beta_l = beta_l_ig[k], 
                              tau = tau_ig[k], lambda = lambda_ig, p=p_ig)
    fcdf_res_ed[k,i] <- ft_cdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_ed[k], beta_l = beta_l_ed[k], 
                              tau = tau_ed[k], lambda = lambda_ed, p=p_ed)
  }
}


faildata <- data.frame()
for(k in 1:9){
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_wie[k,], ft_cdf=fcdf_res_wie[k,], GaAs=paste("Unit # 0", k, sep = ""), model="CPWP"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_gamma[k,], ft_cdf=fcdf_res_gamma[k,], GaAs=paste("Unit # 0", k, sep = ""), model="CPGP"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_ig[k,], ft_cdf=fcdf_res_ig[k,], GaAs=paste("Unit # 0", k, sep = ""), model="CPIGP"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_ed[k,], ft_cdf=fcdf_res_ed[k,], GaAs=paste("Unit # 0", k, sep = ""), model="CPTEDP"))
}

for(k in 10:I){
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_wie[k,], ft_cdf=fcdf_res_wie[k,], GaAs=paste("Unit #", k), model="CPWP"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_gamma[k,], ft_cdf=fcdf_res_gamma[k,], GaAs=paste("Unit #", k), model="CPGP"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_ig[k,], ft_cdf=fcdf_res_ig[k,], GaAs=paste("Unit #", k), model="CPIGP"))
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_ed[k,], ft_cdf=fcdf_res_ed[k,], GaAs=paste("Unit #", k), model="CPTEDP"))
}




#fail_pdf_ed <- 


ft_pdf_plot <- ggplot(data=faildata, aes(x=ft, y= ft_pdf, linetype=model, color=model))+facet_wrap(~GaAs, nrow=2)+
  geom_line(show.legend = T, size=1)+xlab("Failure-time")+ylab("PDF")+ theme_bw()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size = 8),
        legend.key.size=unit(0.3,'cm'),legend.text = element_text(size = 7),
        legend.title=element_blank())
ft_cdf_plot <- ggplot(data=faildata, aes(x=ft, y= ft_cdf, linetype=model, color=model))+facet_wrap(~GaAs, nrow=2)+
  geom_line(show.legend = T, size=1)+xlab("Failure-time")+ylab("CDF")+ theme_bw()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size = 8),
        legend.key.size=unit(0.3,'cm'),legend.text = element_text(size = 7),
        legend.title=element_blank())

postscript("D:/My Document/2023/paper/ED_change_point/figure/simulation/fail_pdf_ig.eps", width=10, height=5)
ft_pdf_plot
dev.off()

postscript("D:/My Document/2023/paper/ED_change_point/figure/simulation/fail_cdf_ig.eps", width=10, height=5)
ft_cdf_plot
dev.off()



## RUL
i=1
z <- seq(0, 40, by = 5)
rul <- seq(0, 80, length=300)
p_true=3
lambda_true=10

rul_pdf_bs <- function(x, t, ff, beta_h, beta_l, tau, lambda, p){
  if(t<tau){
    if(t+x<tau){
      ft <- dbisa(x, scale=(ff-beta_h*t)/beta_h, shape=beta_h^((p-1)/2)/sqrt(lambda*(ff-beta_h*t)))
    } else {
      ft <- dbisa(x, scale=(ff-(beta_h-beta_l)*tau-beta_l*t)/beta_l, 
                  shape=beta_l^((p-1)/2)/sqrt(lambda*((ff-(beta_h-beta_l)*tau-beta_l*t))))
    }
  }else{
    ft <- dbisa(x, scale=(ff-(beta_h-beta_l)*tau-beta_l*t)/beta_l, 
                shape=beta_l^((p-1)/2)/sqrt(lambda*((ff-(beta_h-beta_l)*tau-beta_l*t))))
  }
  return(ft)
}

rul_cdf_bs <- function(x, t, ff, beta_h, beta_l, tau, lambda, p){
  if(t<tau){
    if(t+x<tau){
      ft <- pbisa(x, scale=(ff-beta_h*t)/beta_h, shape=beta_h^((p-1)/2)/sqrt(lambda*(ff-beta_h*t)))
    } else {
      ft <- pbisa(x, scale=(ff-(beta_h-beta_l)*tau-beta_l*t)/beta_l, 
                  shape=beta_l^((p-1)/2)/sqrt(lambda*((ff-(beta_h-beta_l)*tau-beta_l*t))))
    }
  }else{
    ft <- pbisa(x, scale=(ff-(beta_h-beta_l)*tau-beta_l*t)/beta_l, 
                shape=beta_l^((p-1)/2)/sqrt(lambda*((ff-(beta_h-beta_l)*tau-beta_l*t))))
  }
  return(ft)
}

rul3d_df_ed <- rul3d_df_wie <- rul3d_df_ig <- rul3d_df_gamma <-  rul3d_df_true <- data.frame()

for(k in 1:length(z)){
  for(g in 1:length(rul)){
    dens <- rul_pdf_bs(x=rul[g], t=z[k], ff=Df[1], beta_h = beta_h_ed[1], beta_l = beta_l_ed[1], 
                       tau = tau_ed[1], lambda = lambda_ed, p=p_ed)
    dis <- rul_cdf_bs(x=rul[g], t=z[k],  ff=Df[1], beta_h = beta_h_ed[1], beta_l = beta_l_ed[1], 
                      tau = tau_ed[1], lambda = lambda_ed, p=p_ed)
    rul3d_df_ed <- rbind(rul3d_df_ed, data.frame(RUL=rul[g], t=as.character(z[k]), pdf=dens, cdf=dis))
    
    dens <- rul_pdf_bs(x=rul[g], t=z[k], ff=Df[1], beta_h = beta_h_wie[1], beta_l = beta_l_wie[1], 
                       tau = tau_wie[1], lambda = lambda_wie, p=p_wie)
    dis <- rul_cdf_bs(x=rul[g], t=z[k],  ff=Df[1], beta_h = beta_h_wie[1], beta_l = beta_l_wie[1], 
                      tau = tau_wie[1], lambda = lambda_wie, p=p_wie)
    rul3d_df_wie <- rbind(rul3d_df_wie, data.frame(RUL=rul[g], t=as.character(z[k]), pdf=dens, cdf=dis))
    
    dens <- rul_pdf_bs(x=rul[g], t=z[k], ff=Df[1], beta_h = beta_h_gamma[1], beta_l = beta_l_gamma[1], 
                       tau = tau_gamma[1], lambda = lambda_gamma, p=p_gamma)
    dis <- rul_cdf_bs(x=rul[g], t=z[k],  ff=Df[1], beta_h = beta_h_gamma[1], beta_l = beta_l_gamma[1], 
                      tau = tau_gamma[1], lambda = lambda_gamma, p=p_gamma)
    rul3d_df_gamma <- rbind(rul3d_df_gamma, data.frame(RUL=rul[g], t=as.character(z[k]), pdf=dens, cdf=dis))
    
    dens <- rul_pdf_bs(x=rul[g], t=z[k], ff=Df[1], beta_h = beta_h_ig[1], beta_l = beta_l_ig[1], 
                       tau = tau_ig[1], lambda = lambda_ig, p=p_ig)
    dis <- rul_cdf_bs(x=rul[g], t=z[k],  ff=Df[1], beta_h = beta_h_ig[1], beta_l = beta_l_ig[1], 
                      tau = tau_ig[1], lambda = lambda_ig, p=p_ig)
    rul3d_df_ig <- rbind(rul3d_df_ig, data.frame(RUL=rul[g], t=as.character(z[k]), pdf=dens, cdf=dis))
    
    dens <- rul_pdf_bs(x=rul[g], t=z[k], ff=Df[1], beta_h = beta_h_true[1], beta_l = beta_l_true[1], 
                       tau = tau_true[1], lambda = lambda_true, p=p_true)
    dis <- rul_cdf_bs(x=rul[g], t=z[k],  ff=Df[1], beta_h = beta_h_true[1], beta_l = beta_l_true[1], 
                      tau = tau_ig[1], lambda = lambda_true, p=p_true)
    rul3d_df_true <- rbind(rul3d_df_true, data.frame(RUL=rul[g], t=as.character(z[k]), pdf=dens, cdf=dis))
  }
}

## MRL
mrl_wie <- mrl_wie_ig <- mrl_gamma <- mrl_ig <- mrl_ed <- mrl_true <- matrix(NA, nrow=I, ncol=length(z))

for(k in 1:I){
  for(g in 1:length(z)){
    if(z[g]<=tau_wie[k]){
      temp1 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-beta_h_wie[k]*z[g])/beta_h_wie[k], 
                                             shape=beta_h_wie[k]^((p_wie-1)/2)/sqrt(lambda_wie*(Df[k]--beta_h_wie[k]*z[g]) )), 
                         lower = 0, upper = tau_wie[k]-z[g])$value
      temp2 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-(beta_h_wie[k]-beta_l_wie[k])*tau_wie[k]-beta_l_wie[k]*z[g])/beta_l_wie[k], 
                                             shape=beta_l_wie[k]^((p_wie-1)/2)/sqrt(lambda_wie*((Df[k]-(beta_h_wie[k]-beta_l_wie[k])*tau_wie[k]-beta_l_wie[k]*z[g])))),
                         lower = tau_wie[k]-z[g], upper = 1000)$value
      mrl_wie[k,g] <- temp1+temp2}else{
        mrl_wie[k,g] <- beta_l_wie[k]^{(p_wie-2)}/2/lambda_wie+(Df[k]-(beta_h_wie[k]-beta_l_wie[k])*tau_wie[k]-beta_l_wie[k]*z[g])/beta_l_wie[k]
      }
    
    if(z[g]<=tau_gamma[k]){
      temp1 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-beta_h_gamma[k]*z[g])/beta_h_gamma[k], 
                                             shape=beta_h_gamma[k]^((p_gamma-1)/2)/sqrt(lambda_gamma*(Df[k]--beta_h_gamma[k]*z[g]) )), 
                         lower = 0, upper = tau_gamma[k]-z[g])$value
      temp2 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-(beta_h_gamma[k]-beta_l_gamma[k])*tau_gamma[k]-beta_l_gamma[k]*z[g])/beta_l_gamma[k], 
                                             shape=beta_l_gamma[k]^((p_gamma-1)/2)/sqrt(lambda_gamma*((Df[k]-(beta_h_gamma[k]-beta_l_gamma[k])*tau_gamma[k]-beta_l_gamma[k]*z[g])))),
                         lower = tau_gamma[k]-z[g], upper = 1000)$value
      mrl_gamma[k,g] <- temp1+temp2}else{
        mrl_gamma[k,g] <- beta_l_gamma[k]^{(p_gamma-2)}/2/lambda_gamma+(Df[k]-(beta_h_gamma[k]-beta_l_gamma[k])*tau_gamma[k]-beta_l_gamma[k]*z[g])/beta_l_gamma[k]
      }
    
    if(z[g]<=tau_ig[k]){
      temp1 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-beta_h_ig[k]*z[g])/beta_h_ig[k], 
                                             shape=beta_h_ig[k]^((p_ig-1)/2)/sqrt(lambda_ig*(Df[k]--beta_h_ig[k]*z[g]) )), 
                         lower = 0, upper = tau_ig[k]-z[g])$value
      temp2 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-(beta_h_ig[k]-beta_l_ig[k])*tau_ig[k]-beta_l_ig[k]*z[g])/beta_l_ig[k], 
                                             shape=beta_l_ig[k]^((p_ig-1)/2)/sqrt(lambda_ig*((Df[k]-(beta_h_ig[k]-beta_l_ig[k])*tau_ig[k]-beta_l_ig[k]*z[g])))),
                         lower = tau_ig[k]-z[g], upper = 1000)$value
      mrl_ig[k,g] <- temp1+temp2}else{
        mrl_ig[k,g] <- beta_l_ig[k]^{(p_ig-2)}/2/lambda_ig+(Df[k]-(beta_h_ig[k]-beta_l_ig[k])*tau_ig[k]-beta_l_ig[k]*z[g])/beta_l_ig[k]
      }
    
    if(z[g]<=tau_ed[k]){
      temp1 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-beta_h_ed[k]*z[g])/beta_h_ed[k], 
                                             shape=beta_h_ed[k]^((p_ed-1)/2)/sqrt(lambda_ed*(Df[k]--beta_h_ed[k]*z[g]) )), 
                         lower = 0, upper = tau_ed[k]-z[g])$value
      temp2 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k]-beta_l_ed[k]*z[g])/beta_l_ed[k], 
                                             shape=beta_l_ed[k]^((p_ed-1)/2)/sqrt(lambda_ed*((Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k]-beta_l_ed[k]*z[g])))),
                         lower = tau_ed[k]-z[g], upper = 1000)$value
      mrl_ed[k,g] <- temp1+temp2}else{
        mrl_ed[k,g] <- beta_l_ed[k]^{(p_ed-2)}/2/lambda_ed+(Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k]-beta_l_ed[k]*z[g])/beta_l_ed[k]
      }
    
    if(z[g]<=tau_true[k]){
      temp1 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-beta_h_true[k]*z[g])/beta_h_true[k], 
                                             shape=beta_h_true[k]^((p_true-1)/2)/sqrt(lambda_true*(Df[k]--beta_h_true[k]*z[g]) )), 
                         lower = 0, upper = tau_true[k]-z[g])$value
      temp2 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-(beta_h_true[k]-beta_l_true[k])*tau_true[k]-beta_l_true[k]*z[g])/beta_l_true[k], 
                                             shape=beta_l_true[k]^((p_true-1)/2)/sqrt(lambda_true*((Df[k]-(beta_h_true[k]-beta_l_true[k])*tau_true[k]-beta_l_true[k]*z[g])))),
                         lower = tau_true[k]-z[g], upper = 1000)$value
      mrl_true[k,g] <- temp1+temp2}else{
        mrl_true[k,g] <- beta_l_true[k]^{(p_true-2)}/2/lambda_true+(Df[k]-(beta_h_true[k]-beta_l_true[k])*tau_true[k]-beta_l_true[k]*z[g])/beta_l_true[k]
      }
    
  }
}

mrl3 <- round(cbind(mrl_wie[,c(3,5,9)], mrl_gamma[,c(3,5,9)], mrl_ig[,c(3,5,9)], mrl_ed[,c(3,5,9)], mrl_true[,c(3,5,9)]),3)

data_rul <- data.frame(pdf=rep(0,length(z)), t=z, RUL= mrl_true[1,])

rul3d_df <- cbind(rul3d_df_wie, rul3d_df_gamma, rul3d_df_ig, rul3d_df_ed, rul3d_df_true)
names(rul3d_df) <- c("RUL1", "t1", "pdf1", "cdf1", "RUL2", "t2", "pdf2", "cdf2", 
                     "RUL3", "t3", "pdf3", "cdf3", "RUL4", "t4", "pdf4", "cdf4", 
                     "RUL5", "t5", "pdf5", "cdf5")

plot1_pdf <- plot_ly(rul3d_df, x = ~RUL3, y = ~t3, z = ~pdf3,  type = 'scatter3d', mode = 'lines', 
                     color=~t3, line = list(color="pink", dash="solid", width = 3), showlegend = FALSE) %>% 
  layout(scene = list(yaxis = list( title = "Current time" ), xaxis = list( title = "RUL" ), zaxis = list( title = "PDF" )))
plot2_pdf <- plot1_pdf %>%add_trace(rul3d_df, x = ~RUL4, y = ~t4, z = ~pdf4, color=~t4, 
                                line = list(color="green", dash="dash", width =3))
plot_pdf <- plot2_pdf%>% add_markers(data=data_rul, x = ~RUL, y = ~t, z = ~pdf, color = I("black"), symbol=I(8), size=1,
                                     line = list(color="black", dash="dash", width =3, showlegend = FALSE)) 

plot1_cdf <- plot_ly(rul3d_df, x = ~RUL1, y = ~t1, z = ~cdf1,  type = 'scatter3d', mode = 'lines', 
                     color=~t1, line = list(color="pink", dash="solid", width = 4), showlegend = FALSE) %>% 
  layout(scene = list(yaxis = list( title = "Current time" ), xaxis = list( title = "RUL" ), zaxis = list( title = "CDF" )))
plot2_cdf <- plot1_cdf %>%add_trace(rul3d_df, x = ~RUL4, y = ~t4, z = ~cdf4, color=~t4, 
                                    line = list(color="green", dash="dash", width =4))
plot_cdf <- plot2_cdf%>% add_markers(data=data_rul, x = ~RUL, y = ~t, z = ~pdf, color = I("black"), symbol=I(8), size=1,
                                     line = list(color="black", dash="dash", width =3, showlegend = FALSE)) 


data_plot=cbind(rul3d_df_ed[rul3d_df_ed$t%in%c(10,20,40),], model="CPTEDP")
data_plot=rbind(data_plot, cbind(rul3d_df_wie[rul3d_df_wie$t%in%c(10,20,40),], model="CPWP"))
data_plot=rbind(data_plot, cbind(rul3d_df_gamma[rul3d_df_gamma$t%in%c(10,20,40),], model="CPGP"))
data_plot=rbind(data_plot, cbind(rul3d_df_ig[rul3d_df_ig$t%in%c(10,20,40),], model="CPIGP"))

rul_pdf_ed <- ggplot(data=data_plot, aes(x=RUL, y= pdf, linetype = model, shape = model, color=t))+ xlim(0,80)+ylim(0,0.2)+
  geom_line(show.legend = F, size=0.3)+geom_point(size=0.2)+xlab("RUL")+ylab("PDF")+ theme_bw()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size = 8),
        legend.key.size=unit(0.3,'cm'),legend.text = element_text(size = 7),
        legend.title=element_blank())+annotate("text", x=20, y=0.2, label="t=40")+
  annotate("text", x=40, y=0.14, label="t=20")+
  annotate("text", x=50, y=0.125, label="t=10")

rul_cdf_ed <- ggplot(data=data_plot, aes(x=RUL, y= cdf, linetype = model, shape = model, color=t))+ xlim(0,80)+ylim(0,1.05)+
  geom_line(show.legend = F, size=0.3)+geom_point(size=0.2)+xlab("RUL")+ylab("CDF")+ theme_bw()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size = 8),
        legend.key.size=unit(0.3,'cm'),legend.text = element_text(size = 7),
        legend.title=element_blank())+annotate("text", x=20, y=0.8, label="t=40")+
  annotate("text", x=40, y=0.8, label="t=20")+
  annotate("text", x=60, y=0.8, label="t=10")

postscript("D:/My Document/2023/paper/ED_change_point/figure/simulation/rul_pdf_ig.eps", width=4, height=4)
rul_pdf_ed
dev.off()

postscript("D:/My Document/2023/paper/ED_change_point/figure/simulation/rul_cdf_ig.eps", width=4, height=4)
rul_cdf_ed
dev.off()


## mttf
mttf_wie <- mttf_wie_ig <-mttf_gamma <- mttf_ig <-mttf_ed <- numeric(I)

for(k in 1:I){
  temp1 <- integrate(function(t) t*dbisa(t, scale=Df[k]/beta_h_wie[k], shape=beta_h_wie[k]^((p_wie-1)/2)/sqrt(lambda_wie*Df[k])), 
                     lower = 0, upper = tau_wie[k])$value
  temp2 <- integrate(function(t) t*dbisa(t, scale=(Df[k]-(beta_h_wie[k]-beta_l_wie[k])*tau_wie[k])/beta_l_wie[k], 
                                         shape=beta_l_wie[k]^((p_wie-1)/2)/sqrt(lambda_wie*((Df[k]-(beta_h_wie[k]-beta_l_wie[k])*tau_wie[k])))),
                     lower = tau_wie[k], upper = 1000)$value
  mttf_wie[k] <- temp1+temp2
  
  temp1 <- integrate(function(t) t*dbisa(t, scale=Df[k]/beta_h_gamma[k], shape=beta_h_gamma[k]^((p_gamma-1)/2)/sqrt(lambda_gamma*Df[k])), 
                     lower = 0, upper = tau_gamma[k])$value
  temp2 <- integrate(function(t) t*dbisa(t, scale=(Df[k]-(beta_h_gamma[k]-beta_l_gamma[k])*tau_gamma[k])/beta_l_gamma[k], 
                                         shape=beta_l_gamma[k]^((p_gamma-1)/2)/sqrt(lambda_gamma*((Df[k]-(beta_h_gamma[k]-beta_l_gamma[k])*tau_gamma[k])))),
                     lower = tau_gamma[k], upper = 1000)$value
  mttf_gamma[k] <- temp1+temp2
  
  temp1 <- integrate(function(t) t*dbisa(t, scale=Df[k]/beta_h_ig[k], shape=beta_h_ig[k]^((p_ig-1)/2)/sqrt(lambda_ig*Df[k])), 
                     lower = 0, upper = tau_ig[k])$value
  temp2 <- integrate(function(t) t*dbisa(t, scale=(Df[k]-(beta_h_ig[k]-beta_l_ig[k])*tau_ig[k])/beta_l_ig[k], 
                                         shape=beta_l_ig[k]^((p_ig-1)/2)/sqrt(lambda_ig*((Df[k]-(beta_h_ig[k]-beta_l_ig[k])*tau_ig[k])))),
                     lower = tau_ig[k], upper = 1000)$value
  mttf_ig[k] <- temp1+temp2
  
  temp1 <- integrate(function(t) t*dbisa(t, scale=Df[k]/beta_h_ed[k], shape=beta_h_ed[k]^((p_ed-1)/2)/sqrt(lambda_ed*Df[k])), 
                     lower = 0, upper = tau_ed[k])$value
  temp2 <- integrate(function(t) t*dbisa(t, scale=(Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k])/beta_l_ed[k], 
                                         shape=beta_l_ed[k]^((p_ed-1)/2)/sqrt(lambda_ed*((Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k])))),
                     lower = tau_ed[k], upper = 1000)$value
  mttf_ed[k] <- temp1+temp2
}

for(i in 1:I){
  temp1 <- Df[i]/beta_h_wie[i]
  temp2 <- 1-pinvgauss(q=Df[i]^2/(tau_wie[i]*beta_h_wie[i]^2), mean=Df[i]/beta_h_wie[i], 
                       shape=Df[i]^2*lambda_wie)
  temp3 <- (Df[i]-(beta_h_wie[i]-beta_l_wie[i])*tau_wie[i])/beta_l_wie[i]
  temp4 <- pinvgauss(q=(Df[i]-(beta_h_wie[i]-beta_l_wie[i])*tau_wie[i])^2/(tau_wie[i]*beta_l_wie[i]^2), 
                     mean=(Df[i]-(beta_h_wie[i]-beta_l_wie[i])*tau_wie[i])/beta_l_wie[i], 
                     shape=(Df[i]-(beta_h_wie[i]-beta_l_wie[i])*tau_wie[i])^2*lambda_wie)
  mttf_wie_ig[i] <- temp1*temp2+temp3*temp4
}

mttf3 <- rbind(mttf_wie, mttf_wie_ig, mttf_gamma, mttf_ig, mttf_ed)

xtable(rbind(mttf1, mttf2, mttf3), digits = 2)

xtable(rbind(mrl1[,1:12], mrl2[,1:12], mrl3[,1:12]), digits = 2)
