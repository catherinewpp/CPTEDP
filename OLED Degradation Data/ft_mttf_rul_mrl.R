rm(list=ls())
library(R2OpenBUGS)
##BS distribtion
library(VGAM)
## IG distribution
library(statmod)
library(ggplot2)
library(plotly)
library(xtable)

stat_ed <- bugs.log("D:/My Document/Rcode/degradation-ED/ed_change_point/ed_dep/log.txt")[[1]]

## ED
beta_h_ed <- stat_ed[(0:5)*3+1+9,1]
beta_l_ed <- stat_ed[(0:5)*3+2+9,1]
tau_ed <- stat_ed[(0:5)*3+3+9,1]
lambda_ed <- stat_ed["lambda",1]
p_ed <- stat_ed["p",1]

load("D:/My Document/Rcode/degradation_Wie_err/wie_error_real/real.RData")
Df <- sapply(1:6, function(x) y1[n[x]+1])/2
I=6

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
time <- seq(0, 30, length =200)

fpdf_res_ed <- matrix(NA, I, length(time))
for(k in 1:I){
  for(i in 1:length(time)){
    fpdf_res_ed[k,i] <- ft_pdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_ed[k], beta_l = beta_l_ed[k], 
                                  tau = tau_ed[k], lambda = lambda_ed, p=p_ed)
  }
}

fcdf_res_ed <- matrix(NA, I, length(time))
for(k in 1:I){
  for(i in 1:length(time)){
    fcdf_res_ed[k,i] <- ft_cdf_bs(t=time[i], ff=Df[k], beta_h = beta_h_ed[k], beta_l = beta_l_ed[k], 
                                  tau = tau_ed[k], lambda = lambda_ed, p=p_ed)
  }
}


faildata <- data.frame()
for(k in 1:I){
  faildata <- rbind(faildata, data.frame(ft=time, ft_pdf=fpdf_res_ed[k,], ft_cdf=fcdf_res_ed[k,],
                                         GaAs=paste("OLED # 0", k, sep = ""), model="ED"))
}

ft_pdf_plot <- ggplot(data=faildata, aes(x=ft, y= ft_pdf, linetype=GaAs, color=GaAs))+ 
  xlim(4,10)+ylim(0,0.8)+
  geom_line(show.legend = T, size=0.6)+xlab("Failure-time")+ylab("PDF")+ theme_bw()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size = 8),
        legend.key.size=unit(0.3,'cm'),legend.text = element_text(size = 7),
        legend.title=element_blank())

ft_cdf_plot <- ggplot(data=faildata, aes(x=ft, y= ft_cdf, linetype=GaAs, color=GaAs))+xlim(4,10)+ylim(0,1)+
  geom_line(show.legend = T, size=0.6)+xlab("Failure-time")+ylab("CDF")+ theme_bw()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size = 8),
        legend.key.size=unit(0.3,'cm'),legend.text = element_text(size = 7),
        legend.title=element_blank())

postscript("D:/My Document/2023/paper/ED_change_point/figure/case/fail_pdf.eps", width=3, height=3)
ft_pdf_plot
dev.off()

postscript("D:/My Document/2023/paper/ED_change_point/figure/case/fail_cdf.eps", width=3, height=3)
ft_cdf_plot
dev.off()


## RUL
i=1
z <- seq(0, 6, by = 1)
rul <- seq(0, 9, length=300)

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

rul3d_df_ed <- data.frame()

for(k in 1:length(z)){
  for(g in 1:length(rul)){
    dens <- rul_pdf_bs(x=rul[g], t=z[k], ff=Df[1], beta_h = beta_h_ed[1], beta_l = beta_l_ed[1], 
                       tau = tau_ed[1], lambda = lambda_ed, p=p_ed)
    dis <- rul_cdf_bs(x=rul[g], t=z[k],  ff=Df[1], beta_h = beta_h_ed[1], beta_l = beta_l_ed[1], 
                      tau = tau_ed[1], lambda = lambda_ed, p=p_ed)
    rul3d_df_ed <- rbind(rul3d_df_ed, data.frame(RUL=rul[g], t=as.character(z[k]), pdf=dens, cdf=dis))
    }
}

## MRL
mrl_ed <- matrix(NA, nrow=I, ncol=length(z))

for(k in 1:I){
  for(g in 1:length(z)){
    if(z[g]<=tau_ed[k]){
      temp1 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-beta_h_ed[k]*z[g])/beta_h_ed[k], 
                                             shape=beta_h_ed[k]^((p_ed-1)/2)/sqrt(lambda_ed*(Df[k]-beta_h_ed[k]*z[g]) )), 
                         lower = 0.01, upper = tau_ed[k]-z[g])$value
      temp2 <- integrate(function(x) x*dbisa(x, scale=(Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k]-beta_l_ed[k]*z[g])/beta_l_ed[k], 
                                             shape=beta_l_ed[k]^((p_ed-1)/2)/sqrt(lambda_ed*((Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k]-beta_l_ed[k]*z[g])))),
                         lower = tau_ed[k]-z[g], upper = 1000)$value
      mrl_ed[k,g] <- temp1+temp2}else{
        mrl_ed[k,g] <- beta_l_ed[k]^{(p_ed-2)}/2/lambda_ed+(Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k]-beta_l_ed[k]*z[g])/beta_l_ed[k]
      }
    
  }
}


data_rul <- data.frame(pdf=rep(0,length(z)), cdf=rep(0,length(z)), t=z, RUL= mrl_ed[1,])

plot1_pdf <- plot_ly(rul3d_df_ed, x = ~RUL, y = ~t, z = ~pdf,  type = 'scatter3d', mode = 'lines', 
                     color=~t, line = list(color="green", dash="solid", width = 3), showlegend = FALSE) %>% 
  layout(scene = list(yaxis = list( title = "Current time" ), xaxis = list( title = "RUL" ), zaxis = list( title = "PDF" )))
plot_pdf <- plot1_pdf%>% add_markers(data=data_rul, x = ~RUL, y = ~t, z = ~pdf, color = I("black"), symbol=I(8), size=1,
                                     line = list(color="black", dash="dash", width =3, showlegend = FALSE)) 

plot1_cdf <- plot_ly(rul3d_df_ed, x = ~RUL, y = ~t, z = ~cdf,  type = 'scatter3d', mode = 'lines', 
                     color=~t, line = list(color="green", dash="solid", width = 3), showlegend = FALSE) %>% 
  layout(scene = list(yaxis = list( title = "Current time" ), xaxis = list( title = "RUL" ), zaxis = list( title = "CDF" )))
plot_cdf <- plot1_cdf%>% add_markers(data=data_rul, x = ~RUL, y = ~t, z = ~cdf, color = I("black"), symbol=I(8), size=1,
                                     line = list(color="black", dash="dash", width =3, showlegend = FALSE)) 


data_plot=cbind(rul3d_df_ed[rul3d_df_ed$t%in%c(0,3,6),], model="ED")

rul_pdf_ed <- ggplot(data=data_plot, aes(x=RUL, y= pdf,  color=t))+ xlim(0,9)+
  geom_line( size=0.3)+geom_point(size=0.2)+xlab("RUL")+ylab("PDF")+ theme_bw()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size = 8),
        legend.key.size=unit(0.3,'cm'),legend.text = element_text(size = 7),
        legend.title=element_blank())+annotate("text", x=0.6, y=3, label="t=6")+
  annotate("text", x=2.7, y=1.2, label="t=3")+
  annotate("text", x=6.5, y=0.8, label="t=0")

rul_cdf_ed <- ggplot(data=data_plot, aes(x=RUL, y= cdf, color=t))+ xlim(0,9)+ylim(0,1.05)+
  geom_line(show.legend = F, size=0.3)+geom_point(size=0.2)+xlab("RUL")+ylab("PDF")+ theme_bw()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size = 8),
        legend.key.size=unit(0.3,'cm'),legend.text = element_text(size = 7),
        legend.title=element_blank())+annotate("text", x=0.6, y=0.8, label="t=6")+
  annotate("text", x=4, y=0.8, label="t=3")+
  annotate("text", x=6.5, y=0.8, label="t=0")

postscript("D:/My Document/2023/paper/ED_change_point/figure/case/rul_pdf.eps", width=3, height=3)
rul_pdf_ed
dev.off()

postscript("D:/My Document/2023/paper/ED_change_point/figure/case/rul_cdf.eps", width=3, height=3)
rul_cdf_ed
dev.off()

## mttf
mttf_ed <- numeric(I)

for(k in 1:I){
  temp1 <- integrate(function(t) t*dbisa(t, scale=Df[k]/beta_h_ed[k], shape=beta_h_ed[k]^((p_ed-1)/2)/sqrt(lambda_ed*Df[k])), 
                     lower = 0.01, upper = tau_ed[k])$value
  temp2 <- integrate(function(t) t*dbisa(t, scale=(Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k])/beta_l_ed[k], 
                                         shape=beta_l_ed[k]^((p_ed-1)/2)/sqrt(lambda_ed*((Df[k]-(beta_h_ed[k]-beta_l_ed[k])*tau_ed[k])))),
                     lower = tau_ed[k], upper = 100)$value
  mttf_ed[k] <- temp1+temp2
}

# mttf <- numeric(6)
# for(i in 1:6){
#   temp1 <- Df[i]/beta_h_ed[i]
#   temp2 <- 1-pinvgauss(q=Df[i]^2/(tau_ed[i]*beta_h_ed[i]^2), mean=Df[i]/beta_h_ed[i], 
#                        shape=Df[i]^2*lambda_ed)
#   temp3 <- (Df[i]-(beta_h_ed[i]-beta_l_ed[i])*tau_ed[i])/beta_l_ed[i]
#   temp4 <- pinvgauss(q=(Df[i]-(beta_h_ed[i]-beta_l_ed[i])*tau_ed[i])^2/(tau_ed[i]*beta_l_ed[i]^2), 
#                      mean=(Df[i]-(beta_h_ed[i]-beta_l_ed[i])*tau_ed[i])/beta_l_ed[i], 
#                      shape=(Df[i]-(beta_h_ed[i]-beta_l_ed[i])*tau_ed[i])^2*lambda_ed)
#   mttf[i] <- temp1*temp2+temp3*temp4
# }

