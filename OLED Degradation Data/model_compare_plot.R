rm(list=ls())
library(ggplot2)
library(scales)
library(grid)
library(ggpubr)## ggarange
library(R2OpenBUGS)

source("D:/My Document/bugscode/Change-point-Codes/R2OpenBUGS/led/read_led.R")

est_ed <- bugs.log("D:/My Document/Rcode/degradation-ED/ed_change_point/ed_dep/log.txt")$stats
est_gp_d <- read.table("D:/My Document/bugscode/Change-point-Codes/Codes/led_data/gp_con_dep.txt")
est_two <- read.table("D:/My Document/bugscode/Change-point-Codes/Codes/led_data/beta_two.txt")
est_bi <- read.table("D:/My Document/bugscode/Change-point-Codes/Codes/led_data/gamma_bi.txt")
tm_p <- seq(0,8,by=0.02)
y_gp_d <- y_ed <- y_bi <- y_two <- numeric(length(tm_p))

#######  Winener Process Gaussian Process with Dependent Prior Assumption
i <- 1
for(j in 1:length(tm_p)){
  y_gp_d[j] <- est_gp_d[(i-1)*3+1,1]*min(tm_p[j],est_gp_d[(i-1)*3+3,1]) + est_gp_d[(i-1)*3+2,1]*max(tm_p[j]-est_gp_d[(i-1)*3+3,1],0)
  y_ed[j] <- -est_ed[(i-1)*3+1+9,1]*min(tm_p[j],est_ed[(i-1)*3+3+9,1]) - est_ed[(i-1)*3+2+9,1]*max(tm_p[j]-est_ed[(i-1)*3+3+9,1],0)
}

data_plot1 <- data.frame(xx=exp(tm1[(n[i]+1):(n[i+1])]), yy=y1[(n[i]+1):(n[i+1])], OLED=paste("OLED #", i))
data_plot_a_d <- data.frame(xx=exp(tm_p+tm1[n[i]+1]), yy=y_gp_d+y1[n[i]+1], OLED=paste("OLED #", i))
data_plot_a_ed <- data.frame(xx=exp(tm_p+tm1[n[i]+1]), yy=y_ed+y1[n[i]+1], OLED=paste("OLED #", i))


for(i in 2:m){
  for(j in 1:length(tm_p)){
    y_gp_d[j] <- est_gp_d[(i-1)*3+1,1]*min(tm_p[j],est_gp_d[(i-1)*3+3,1]) + est_gp_d[(i-1)*3+2,1]*max(tm_p[j]-est_gp_d[(i-1)*3+3,1],0)
    y_ed[j] <- -est_ed[(i-1)*3+1+9,1]*min(tm_p[j],est_ed[(i-1)*3+3+9,1]) - est_ed[(i-1)*3+2+9,1]*max(tm_p[j]-est_ed[(i-1)*3+3+9,1],0)
  }
  # data_plot1 <- rbind(data_plot1, data.frame(xx=tm1[(n[i]+1):(n[i+1])]-tm1[n[i]+1],yy=y1[(n[i]+1):(n[i+1])],OLED=paste("OLED #", i)))
  # data_plot_a_d <- rbind(data_plot_a_d, data.frame(xx=tm_p, yy=y_gp_d+y1[n[i]+1], OLED=paste("OLED #", i)))
  # data_plot_a_ed <- rbind(data_plot_a_ed, data.frame(xx=tm_p, yy=y_ed+y1[n[i]+1], OLED=paste("OLED #", i)))
   data_plot1 <- rbind(data_plot1, data.frame(xx=exp(tm1[(n[i]+1):(n[i+1])]),yy=y1[(n[i]+1):(n[i+1])],OLED=paste("OLED #", i)))
   data_plot_a_d <- rbind(data_plot_a_d, data.frame(xx=exp(tm_p+tm1[n[i]+1]), yy=y_gp_d+y1[n[i]+1], OLED=paste("OLED #", i)))
   data_plot_a_ed <- rbind(data_plot_a_ed, data.frame(xx=exp(tm_p+tm1[n[i]+1]), yy=y_ed+y1[n[i]+1], OLED=paste("OLED #", i)))
  
}

p_d <- ggplot(data_plot1, aes(x=xx, y=yy, color=OLED, shape=OLED,linetype=OLED))+
  geom_point(size=0.5)+ coord_trans(x="log")+
  geom_line(data=data_plot_a_d, aes(x=xx, y=yy, color=OLED),size=0.2, show.legend = FALSE)+
  xlab("Time")+ylab("Luminorsity") +ggtitle("CPWP")+scale_x_continuous(breaks=c(0,500,1000,2000))+
  theme(title= element_text(size=8), axis.text=element_text(hjust=2,size=6), 
        axis.title=element_text(size = 8),legend.key.size=unit(0.3,'cm'),
        legend.text = element_text(size = 7),legend.title=element_blank())


p_ed <- ggplot(data_plot1, aes(x=xx, y=yy, color=OLED, shape=OLED,linetype=OLED))+
  geom_point(size=0.5,show.legend = FALSE)+ coord_trans(x="log")+
  geom_line(data=data_plot_a_ed, aes(x=xx, y=yy, color=OLED),size=0.2,show.legend = FALSE)+
  xlab("Time")+ylab("Luminorsity") +ggtitle("CPTEDP")+scale_x_continuous(breaks=c(0,500,1000,2000))+
  theme(title= element_text(size=8), axis.text=element_text(hjust=2,size=6), 
        axis.title=element_text(size = 8),legend.key.size=unit(0.3,'cm'),
        legend.text = element_text(size = 7),legend.title=element_blank())


#######  Two Phase 
i <- 1
for(j in 1:length(tm_p)){
  y_two[j] <- est_two[(i+2)*3+1,1]*tm_p[j]-est_two[(i+2)*3+2,1]*min(tm_p[j],est_two[(i+2)*3+3,1]) 
}

data_plot2 <- data.frame(xx=exp(tm1[(n[i]+1):(n[i+1])]),yy=y1[(n[i]+1):(n[i+1])],OLED=paste("OLED #", i))
data_plot_a_two <- data.frame(xx=exp(tm_p+tm1[n[i]+1]), yy=y_two+y1[n[i]+1], OLED=paste("OLED #", i))

for(i in 2:m){
  for(j in 1:length(tm_p)){
    y_two[j] <- est_two[(i+2)*3+1,1]*tm_p[j]-est_two[(i+2)*3+2,1]*min(tm_p[j],est_two[(i+2)*3+3,1]) 
  }
  data_plot2 <- rbind(data_plot2, data.frame(xx=exp(tm1[(n[i]+1):(n[i+1])]),yy=y1[(n[i]+1):(n[i+1])], OLED=paste("OLED #", i)))
  data_plot_a_two <- rbind(data_plot_a_two, data.frame(xx=exp(tm_p+tm1[n[i]+1]), yy=y_two+y1[n[i]+1], OLED=paste("OLED #", i)))  
}


p_two <- ggplot(data_plot2,aes(x=xx, y=yy, color=OLED, shape=OLED))+
  geom_point(size=0.5,show.legend = FALSE)+ coord_trans(x="log")+
  geom_line(data=data_plot_a_two, aes(x=xx, y=yy, color=OLED),size=0.2,show.legend = FALSE)+
  xlab("Time")+ylab("Luminorsity") +ggtitle("TPLCP")+scale_x_continuous(breaks=c(0,500,1000,2000))+
  theme(title= element_text(size=8), axis.text=element_text(hjust=2,size=6), 
        axis.title=element_text(size = 8),legend.key.size=unit(0.3,'cm'),
        legend.text = element_text(size = 7),legend.title=element_blank())
## CPWP
est_gp_in <- read.table("D:/My Document/bugscode/Change-point-Codes/Codes/led_data/gp_con_inde.txt")
y_gp_in <-  real_in <-numeric(length(tm_p))
yy_in<- numeric(length(tm11))
i <- 1
for(j in 1:length(tm_p)){
  y_gp_in[j] <- est_gp_in[(i-1)*3+1,1]*min(tm_p[j],est_gp_in[(i-1)*3+3,1]) + est_gp_in[(i-1)*3+2,1]*max(tm_p[j]-est_gp_in[(i-1)*3+3,1],0)
  real_in[j] <- est_gp_in[10,1]*min(tm_p[j],est_gp_in[19,1]) + est_gp_in[21,1]*max(tm_p[j]-est_gp_in[19,1],0)
}

data_plot_a_in <- data.frame(xx=exp(tm_p+tm1[n[i]+1]), yy=y_gp_in+y1[n[i]+1],OLED=paste("OLED #", i))
data_real_in <- data.frame(xx=exp(tm_p+tm1[n[i]+1]), yy=real_in+y1[n[i]+1])


for(i in 2:m){
  for(j in 1:length(tm_p)){
    y_gp_in[j] <- est_gp_in[(i-1)*3+1,1]*min(tm_p[j],est_gp_in[(i-1)*3+3,1]) + est_gp_in[(i-1)*3+2,1]*max(tm_p[j]-est_gp_in[(i-1)*3+3,1],0)
  }
  data_plot_a_in <- rbind(data_plot_a_in, data.frame(xx=exp(tm_p+tm1[n[i]+1]), yy=y_gp_in+y1[n[i]+1],OLED=paste("OLED #", i)))
}

p_in <- ggplot(data_plot2,aes(x=xx, y=yy, color=OLED, shape=OLED))+
  geom_point(size=0.5,show.legend = FALSE)+ coord_trans(x="log")+
  geom_line(data=data_plot_a_in, aes(x=xx, y=yy, color=OLED),size=0.2,show.legend = FALSE)+
  xlab("Time")+ylab("Luminorsity") +ggtitle("CPWP")+scale_x_continuous(breaks=c(0,500,1000,2000))+
  theme(title= element_text(size=8), axis.text=element_text(hjust=2,size=6), 
        axis.title=element_text(size = 8),legend.key.size=unit(0.3,'cm'),
        legend.text = element_text(size = 7),legend.title=element_blank())


#######  Bi-exponential Phase 
y_bi <- real_bi <- numeric(length(tm_p))
i <- 1
for(j in 1:length(tm_p)){
  y_bi[j] <- est_bi[(i-1)*3+1,5]*exp(-(est_bi[(i-1)*3+2,5]+est_bi[(i-1)*3+3,5])*tm_p[j])+(1-est_bi[(i-1)*3+1,5])*exp(-est_bi[(i-1)*3+2,5]*tm_p[j])
  real_bi[j]<- est_bi[(6)*3+1,5]*exp(-(est_bi[(6)*3+2,5]+est_bi[(6)*3+3,5])*tm_p[j])+(1-est_bi[(6)*3+1,5])*exp(-est_bi[(6)*3+2,5]*tm_p[j])
}

data_plot3 <- data.frame(xx=tm2[(n[i]+1):(n[i+1])],yy=y1[(n[i]+1):(n[i+1])],OLED=paste("OLED #", i))
data_plot_a_bi <- data.frame(xx=tm_p*100+tm2[n[i]+1], yy=log(y_bi)*50+y1[n[i]+1],OLED=paste("OLED #", i))

for(i in 2:m){
  for(j in 1:length(tm_p)){
    y_bi[j] <- est_bi[(i-1)*3+1,5]*exp(-(est_bi[(i-1)*3+2,5]+est_bi[(i-1)*3+3,5])*tm_p[j])+(1-est_bi[(i-1)*3+1,5])*exp(-est_bi[(i-1)*3+2,5]*tm_p[j])
  }
  
  data_plot3 <- rbind(data_plot3, data.frame(xx=tm2[(n[i]+1):(n[i+1])],yy=y1[(n[i]+1):(n[i+1])],OLED=paste("OLED #", i)))
  data_plot_a_bi <- rbind(data_plot_a_bi, data.frame(xx=tm_p*100+tm2[n[i]+1], yy=log(y_bi)*50+y1[n[i]+1],OLED=paste("OLED #", i)))
  #data_log_bi <- rbind(data_log_bi, data.frame(xx=log(tm_p*100+tm2[n[i]+1]), yy=log(y_bi)*50+y1[n[i]+1],OLED=paste("OLED #", i)))
}

data_plot3$log_xx <- rep(c( 0, cumsum(diff(data_plot3$xx[1:(nrow(data_plot3)/6)]))),6)
data_plot_a_bi$log_xx <- rep(c( 0, cumsum(diff(data_plot_a_bi$xx[1:(nrow(data_plot_a_bi)/6)]))),6)


p_bi <- ggplot(data_plot3, aes(x= xx, y=yy, shape=OLED, color=OLED, linetype=OLED))+
  geom_point(size=0.5,show.legend = FALSE)+
  geom_line(data=data_plot_a_bi, aes(x=xx, y=yy),size=0.2,show.legend = FALSE)+ coord_trans(x="log")+
  xlab("Time")+ylab("Luminorsity")+ggtitle("BE")+scale_x_continuous(breaks=c(0,500,1000,2000))+
  theme(title= element_text(size=8), axis.text=element_text(hjust=2,size=6), 
        axis.title=element_text(size = 8),legend.key.size=unit(0.3,'cm'),
        legend.text = element_text(size = 7),legend.title=element_blank())

##### CPWPME
# load("D:/My Document/Rcode/degradation_Wie_err/wie_error_real/real.RData")
# 
# tm_p <- seq(0,8,by=0.02)
# fit <- matrix(NA, m, length(tm_p))
# for(i in 1:m){
#   for(j in 1:length(tm_p)){
#     fit[i,j] <-  -summary[(i-1)*3+10, 1]*min(summary[(i-1)*3+12, 1],tm_p[j]) -(summary[(i-1)*3+11, 1])*max(tm_p[j]-summary[(i-1)*3+12, 1],0)
#   }
#   #lines(tm11[i,], fit[i,]+y1[n[i]+1], col=i)
# }
# 
# i <- 1
# data_origin <- data.frame(xx=exp(tm1[(n[i]+1):(n[i+1])]) , yy=(y1_two[i,]+y1[n[i]+1]), OLED=paste("OLED #", i))
# data_fit <- data.frame(xx=exp(tm_p+tm1[n[i]+1]), yy=(fit[i,]+y1[n[i]+1]), OLED=paste("OLED #", i))
# 
# for(i in 2:6){
#   data_origin <- rbind(data_origin, data.frame(xx=exp(tm1[(n[i]+1):(n[i+1])]), yy=(y1_two[i,]+y1[n[i]+1]), OLED=paste("OLED #", i)))
#   data_fit <- rbind(data_fit, data.frame(xx=exp(tm_p+tm1[n[i]+1]), yy=(fit[i,]+y1[n[i]+1]), OLED=paste("OLED #", i)))
# }
# 
# p_d <- ggplot(data_origin,aes(x=xx, y=yy, color=OLED, shape=OLED,linetype=OLED))+
#   geom_point(size=0.5,show.legend = FALSE)+ylab("Luminorsity")+coord_trans(x="log") +
#   theme(title= element_text(size=8),axis.text=element_text(size=8), 
#         axis.title=element_text(size = 8),legend.key.size=unit(0.3,'cm'),
#         legend.text = element_text(size = 7),legend.title=element_blank())+
#   geom_line(data=data_fit, aes(x=xx, y=yy, color=OLED),size=0.2,show.legend = FALSE)+
#   xlab("Time")+ggtitle("CPWP")


postscript("D:/My Document/2023/paper/ED_change_point/figure/case/fit_compare.eps", width = 8, height=4)
grid.newpage()  ##?½?ҳ??
pushViewport(viewport(layout = grid.layout(2,2))) ####??ҳ???ֳ?2*2????
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
print(p_ed, vp = vplayout(1,1))   
print(p_in, vp = vplayout(1,2))
print(p_two, vp = vplayout(2,1))  
print(p_bi, vp = vplayout(2,2))  
dev.off()

model1 <- ind_fit1 <- data.frame()
for(g in 1:6){
  model1 <- rbind(model1, data_origin[which(data_origin$OLED==paste("OLED # ", g, sep="")),])
  
  fit1 <- data_fit[which(data_fit$OLED==paste("OLED # ", g, sep="")),]
  fit1 <- cbind(fit1, model=rep("CPWP", nrow(fit1)))
  
  fit2 <- data_plot_a_ed[which(data_plot_a_ed$OLED==paste("OLED # ", g, sep="")),]
  fit2 <- cbind(fit2, model=rep("CPTEDP", nrow(fit2)))
  
  fit3 <- data_plot_a_two[which(data_plot_a_two$OLED==paste("OLED # ", g, sep="")),]
  fit3 <- cbind(fit3, model=rep("TPLCP", nrow(fit3)))
  
  # fit4 <- data_log_bi[which(data_log_bi$OLED==paste("OLED # ", g, sep="")),]
  # fit4 <- cbind(fit4, model=rep("BE", nrow(fit4)))
  
  ind_fit1 <- rbind(ind_fit1, rbind( fit2, fit1, fit3))
}


# for(g in 1:6){
#   model1 <- data_origin[which(data_origin$OLED==paste("OLED # ", g, sep="")),]
#   
#   fit1 <- data_fit[which(data_fit$OLED==paste("OLED # ", g, sep="")),]
#   fit1 <- cbind(fit1, model=rep("CPWP", nrow(fit1)))
#   
#   fit2 <- data_plot_a_ed[which(data_plot_a_ed$OLED==paste("OLED # ", g, sep="")),]
#   fit2 <- cbind(fit2, model=rep("CPTEDP", nrow(fit2)))
#   
#   fit3 <- data_plot_a_two[which(data_plot_a_two$OLED==paste("OLED # ", g, sep="")),]
#   fit3 <- cbind(fit3, model=rep("TPLCP", nrow(fit3)))
#   
#   # fit4 <- data_log_bi[which(data_log_bi$OLED==paste("OLED # ", g, sep="")),]
#   # fit4 <- cbind(fit4, model=rep("BE", nrow(fit4)))
#   
#   ind_fit1 <- rbind( fit2, fit1, fit3)
#   
#   p_d1[[g]] <- ggplot(model1, aes(x=xx, y=yy))+
#     geom_point(size=1)+ylab("Luminorsity") + ylim(20,80)+
#     theme(legend.position = 'right', axis.text=element_text(size=7), 
#           plot.title = element_text(size = 10),axis.title=element_text(size = 9), 
#           legend.key.size=unit(0.4,'cm'),legend.text = element_text(size = 8),
#           legend.title=element_blank())+
#     geom_line(data=ind_fit1, aes(x=xx, y=yy, color=model, linetype = model),
#               size=0.4)+xlim(0,8)+xlab("Time")+
#     ggtitle(paste("OLED # ", g, sep=""))
#   
# }

postscript("D:/My Document/2023/paper/ED_change_point/figure/case/fit_unit.eps", width = 8, height=4)
# ggarrange(
#   p_d1[[1]], p_d1[[2]], p_d1[[3]],
#   p_d1[[4]], p_d1[[5]], p_d1[[6]],
#   ncol = 2, nrow = 3, 
#   common.legend = T
# )
ggplot(model1, aes(x=xx, y=yy))+ facet_wrap(model1$OLED)+
  geom_point(size=1)+ylab("Luminorsity") + ylim(20,80)+
  theme(legend.position = 'right', axis.text=element_text(size=7), 
        plot.title = element_text(size = 10),axis.title=element_text(size = 9), 
        legend.key.size=unit(0.4,'cm'),legend.text = element_text(size = 8),
        legend.title=element_blank())+
  geom_line(data=ind_fit1, aes(x=xx, y=yy, color=model, linetype = model),
            size=0.4)+xlim(0,8)+xlab("Time")

dev.off()
