rm(list=ls())
library(R2OpenBUGS)
##BS distribtion
library(VGAM)
## IG distribution
library(statmod)
library(ggplot2)
library(xtable)
beta_h_true <- c(6.719841, 7.082128, 6.626296, 7.713431, 7.147360, 6.633075, 7.217985, 7.330189, 7.257497, 6.863426)
beta_l_true <- c(2.478067, 2.123279, 1.803546, 1.299650, 2.355734, 1.985791, 1.994880, 2.298467, 2.259693, 2.187808)
tau_true <- c(15.50334, 15.42839, 15.04084, 13.91039, 15.33949, 14.96926, 14.91467, 14.19444, 14.73811, 15.22892)


## Wiener process
lambda_true <- 4
para_true1 <- para_true <- c(beta_h_true, beta_l_true, tau_true, lambda_true)
  
est_wie <- est_gamma <- est_ig <- est_ed <- list()
for(g in 1:100){
  est_wie[[g]] <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/wiener/","bugs_wie",g, "/", "log.txt", sep=""))[[1]]
  est_gamma[[g]] <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/wiener/","bugs_gamma",g, "/", "log.txt", sep=""))[[1]]
  est_ig[[g]] <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/wiener/","bugs_ig",g, "/", "log.txt", sep=""))[[1]]
  est_ed[[g]] <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/wiener/","bugs_ed",g, "/", "log.txt", sep=""))[[1]]
}


est_arr <- array(NA, dim=c(31, 12, 100))
for(g in 1:100){
  est_arr[1:10, 1:3, g] <- est_wie[[g]][(0:9)*3+1, c(1, 3, 5)] 
  est_arr[11:20, 1:3, g] <- est_wie[[g]][(0:9)*3+2, c(1, 3, 5)]
  est_arr[21:30, 1:3, g] <- est_wie[[g]][(0:9)*3+3, c(1, 3, 5)] 
  est_arr[31, 1:3, g] <- est_wie[[g]][31, c(1, 3, 5)] 
  
  est_arr[1:10, 4:6, g] <- est_gamma[[g]][(0:9)*3+1, c(1, 3, 5)] 
  est_arr[11:20, 4:6, g] <- est_gamma[[g]][(0:9)*3+2, c(1, 3, 5)]
  est_arr[21:30, 4:6, g] <- est_gamma[[g]][(0:9)*3+3, c(1, 3, 5)] 
  est_arr[31, 4:6, g] <- est_gamma[[g]][31, c(1, 3, 5)] 

  est_arr[1:10, 7:9, g] <- est_ig[[g]][(0:9)*3+1, c(1, 3, 5)] 
  est_arr[11:20, 7:9, g] <- est_ig[[g]][(0:9)*3+2, c(1, 3, 5)]
  est_arr[21:30, 7:9, g] <- est_ig[[g]][(0:9)*3+3, c(1, 3, 5)] 
  est_arr[31, 7:9, g] <- est_ig[[g]][31, c(1, 3, 5)] 
  
  est_arr[1:10, 10:12, g] <- est_ed[[g]][(0:9)*3+1, c(1, 3, 5)] 
  est_arr[11:20, 10:12, g] <- est_ed[[g]][(0:9)*3+2, c(1, 3, 5)]
  est_arr[21:30, 10:12, g] <- est_ed[[g]][(0:9)*3+3, c(1, 3, 5)] 
  est_arr[31, 10:12, g] <- est_ed[[g]][31, c(1, 3, 5)] 
}

est_summary <- matrix(NA, nrow=31, ncol=12)

est_summary[,1] <- apply(apply(est_arr[, 1, ], 2, function(x) abs(x-para_true) ), 1, mean)
est_summary[,4] <- apply(apply(est_arr[, 4, ], 2, function(x) abs(x-para_true) ), 1, mean)
est_summary[,7] <- apply(apply(est_arr[, 7, ], 2, function(x) abs(x-para_true) ), 1, mean)
est_summary[,10] <- apply(apply(est_arr[, 10, ], 2, function(x) abs(x-para_true) ), 1, mean)

est_summary[,2] <- apply(apply(est_arr[, 1, ], 2, function(x) (x-para_true)^2 ), 1, mean)
est_summary[,5] <- apply(apply(est_arr[, 4, ], 2, function(x) (x-para_true)^2 ), 1, mean)
est_summary[,8] <- apply(apply(est_arr[, 7, ], 2, function(x) (x-para_true)^2 ), 1, mean)
est_summary[,11] <- apply(apply(est_arr[, 10, ], 2, function(x) (x-para_true)^2 ), 1, mean)

est_summary[,3] <- apply(apply(est_arr[, 2, ], 2, function(x) x<para_true )&apply(est_arr[, 3, ], 2, function(x) x>para_true ), 1, mean)
est_summary[,6] <- apply(apply(est_arr[, 5, ], 2, function(x) x<para_true )&apply(est_arr[, 6, ], 2, function(x) x>para_true ), 1, mean)
est_summary[,9] <- apply(apply(est_arr[, 8, ], 2, function(x) x<para_true )&apply(est_arr[, 9, ], 2, function(x) x>para_true ), 1, mean)
est_summary[,12] <- apply(apply(est_arr[, 11, ], 2, function(x) x<para_true)&apply(est_arr[, 12, ], 2, function(x) x>para_true), 1, mean)

est_summary1 <- est_summary



## Gamma process
lambda_true <- 1
para_true2 <- para_true <- c(beta_h_true, beta_l_true, tau_true, lambda_true)

est_wie <- est_gamma <- est_ig <- est_ed <- list()
for(g in 1:100){
  est_wie[[g]] <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/gamma/","bugs_wie",g, "/", "log.txt", sep=""))[[1]]
  est_gamma[[g]] <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/gamma/","bugs_gamma",g, "/", "log.txt", sep=""))[[1]]
  est_ig[[g]] <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/gamma/","bugs_ig",g, "/", "log.txt", sep=""))[[1]]
  est_ed[[g]] <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/gamma/","bugs_ed",g, "/", "log.txt", sep=""))[[1]]
}


est_arr <- array(NA, dim=c(31, 12, 100))
for(g in 1:100){
  est_arr[1:10, 1:3, g] <- est_wie[[g]][(0:9)*3+1, c(1, 3, 5)] 
  est_arr[11:20, 1:3, g] <- est_wie[[g]][(0:9)*3+2, c(1, 3, 5)]
  est_arr[21:30, 1:3, g] <- est_wie[[g]][(0:9)*3+3, c(1, 3, 5)] 
  est_arr[31, 1:3, g] <- est_wie[[g]][31, c(1, 3, 5)] 
  
  est_arr[1:10, 4:6, g] <- est_gamma[[g]][(0:9)*3+1, c(1, 3, 5)] 
  est_arr[11:20, 4:6, g] <- est_gamma[[g]][(0:9)*3+2, c(1, 3, 5)]
  est_arr[21:30, 4:6, g] <- est_gamma[[g]][(0:9)*3+3, c(1, 3, 5)] 
  est_arr[31, 4:6, g] <- est_gamma[[g]][31, c(1, 3, 5)] 
  
  est_arr[1:10, 7:9, g] <- est_ig[[g]][(0:9)*3+1, c(1, 3, 5)] 
  est_arr[11:20, 7:9, g] <- est_ig[[g]][(0:9)*3+2, c(1, 3, 5)]
  est_arr[21:30, 7:9, g] <- est_ig[[g]][(0:9)*3+3, c(1, 3, 5)] 
  est_arr[31, 7:9, g] <- est_ig[[g]][31, c(1, 3, 5)] 
  
  est_arr[1:10, 10:12, g] <- est_ed[[g]][(0:9)*3+1, c(1, 3, 5)] 
  est_arr[11:20, 10:12, g] <- est_ed[[g]][(0:9)*3+2, c(1, 3, 5)]
  est_arr[21:30, 10:12, g] <- est_ed[[g]][(0:9)*3+3, c(1, 3, 5)] 
  est_arr[31, 10:12, g] <- est_ed[[g]][31, c(1, 3, 5)] 
}

est_summary <- matrix(NA, nrow=31, ncol=12)

est_summary[,1] <- apply(apply(est_arr[, 1, ], 2, function(x) abs(x-para_true) ), 1, mean)
est_summary[,4] <- apply(apply(est_arr[, 4, ], 2, function(x) abs(x-para_true) ), 1, mean)
est_summary[,7] <- apply(apply(est_arr[, 7, ], 2, function(x) abs(x-para_true) ), 1, mean)
est_summary[,10] <- apply(apply(est_arr[, 10, ], 2, function(x) abs(x-para_true) ), 1, mean)

est_summary[,2] <- apply(apply(est_arr[, 1, ], 2, function(x) (x-para_true)^2 ), 1, mean)
est_summary[,5] <- apply(apply(est_arr[, 4, ], 2, function(x) (x-para_true)^2 ), 1, mean)
est_summary[,8] <- apply(apply(est_arr[, 7, ], 2, function(x) (x-para_true)^2 ), 1, mean)
est_summary[,11] <- apply(apply(est_arr[, 10, ], 2, function(x) (x-para_true)^2 ), 1, mean)

est_summary[,3] <- apply(apply(est_arr[, 2, ], 2, function(x) x<para_true )&apply(est_arr[, 3, ], 2, function(x) x>para_true ), 1, mean)
est_summary[,6] <- apply(apply(est_arr[, 5, ], 2, function(x) x<para_true )&apply(est_arr[, 6, ], 2, function(x) x>para_true ), 1, mean)
est_summary[,9] <- apply(apply(est_arr[, 8, ], 2, function(x) x<para_true )&apply(est_arr[, 9, ], 2, function(x) x>para_true ), 1, mean)
est_summary[,12] <- apply(apply(est_arr[, 11, ], 2, function(x) x<para_true)&apply(est_arr[, 12, ], 2, function(x) x>para_true), 1, mean)


est_summary2 <- est_summary


## IG process
lambda_true <- 10
para_true3 <- para_true <- c(beta_h_true, beta_l_true, tau_true, lambda_true)

est_wie <- est_gamma <- est_ig <- est_ed <- list()
for(g in 1:100){
  est_wie[[g]] <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/ig/","bugs_wie",g, "/", "log.txt", sep=""))[[1]]
  est_gamma[[g]] <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/ig/","bugs_gamma",g, "/", "log.txt", sep=""))[[1]]
  est_ig[[g]] <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/ig/","bugs_ig",g, "/", "log.txt", sep=""))[[1]]
  est_ed[[g]] <- bugs.log(paste("D:/My Document/Rcode/degradation-ED/ed_change_point/simulation/ig/","bugs_ed",g, "/", "log.txt", sep=""))[[1]]
}


est_arr <- array(NA, dim=c(31, 12, 100))
for(g in 1:100){
  est_arr[1:10, 1:3, g] <- est_wie[[g]][(0:9)*3+1, c(1, 3, 5)] 
  est_arr[11:20, 1:3, g] <- est_wie[[g]][(0:9)*3+2, c(1, 3, 5)]
  est_arr[21:30, 1:3, g] <- est_wie[[g]][(0:9)*3+3, c(1, 3, 5)] 
  est_arr[31, 1:3, g] <- est_wie[[g]][31, c(1, 3, 5)] 
  
  est_arr[1:10, 4:6, g] <- est_gamma[[g]][(0:9)*3+1, c(1, 3, 5)] 
  est_arr[11:20, 4:6, g] <- est_gamma[[g]][(0:9)*3+2, c(1, 3, 5)]
  est_arr[21:30, 4:6, g] <- est_gamma[[g]][(0:9)*3+3, c(1, 3, 5)] 
  est_arr[31, 4:6, g] <- est_gamma[[g]][31, c(1, 3, 5)] 
  
  est_arr[1:10, 7:9, g] <- est_ig[[g]][(0:9)*3+1, c(1, 3, 5)] 
  est_arr[11:20, 7:9, g] <- est_ig[[g]][(0:9)*3+2, c(1, 3, 5)]
  est_arr[21:30, 7:9, g] <- est_ig[[g]][(0:9)*3+3, c(1, 3, 5)] 
  est_arr[31, 7:9, g] <- est_ig[[g]][31, c(1, 3, 5)] 
  
  est_arr[1:10, 10:12, g] <- est_ed[[g]][(0:9)*3+1, c(1, 3, 5)] 
  est_arr[11:20, 10:12, g] <- est_ed[[g]][(0:9)*3+2, c(1, 3, 5)]
  est_arr[21:30, 10:12, g] <- est_ed[[g]][(0:9)*3+3, c(1, 3, 5)] 
  est_arr[31, 10:12, g] <- est_ed[[g]][31, c(1, 3, 5)] 
}

est_summary <- matrix(NA, nrow=31, ncol=12)

est_summary[,1] <- apply(apply(est_arr[, 1, ], 2, function(x) abs(x-para_true) ), 1, mean)
est_summary[,4] <- apply(apply(est_arr[, 4, ], 2, function(x) abs(x-para_true) ), 1, mean)
est_summary[,7] <- apply(apply(est_arr[, 7, ], 2, function(x) abs(x-para_true) ), 1, mean)
est_summary[,10] <- apply(apply(est_arr[, 10, ], 2, function(x) abs(x-para_true) ), 1, mean)

est_summary[,2] <- apply(apply(est_arr[, 1, ], 2, function(x) (x-para_true)^2 ), 1, mean)
est_summary[,5] <- apply(apply(est_arr[, 4, ], 2, function(x) (x-para_true)^2 ), 1, mean)
est_summary[,8] <- apply(apply(est_arr[, 7, ], 2, function(x) (x-para_true)^2 ), 1, mean)
est_summary[,11] <- apply(apply(est_arr[, 10, ], 2, function(x) (x-para_true)^2 ), 1, mean)

est_summary[,3] <- apply(apply(est_arr[, 2, ], 2, function(x) x<para_true )&apply(est_arr[, 3, ], 2, function(x) x>para_true ), 1, mean)
est_summary[,6] <- apply(apply(est_arr[, 5, ], 2, function(x) x<para_true )&apply(est_arr[, 6, ], 2, function(x) x>para_true ), 1, mean)
est_summary[,9] <- apply(apply(est_arr[, 8, ], 2, function(x) x<para_true )&apply(est_arr[, 9, ], 2, function(x) x>para_true ), 1, mean)
est_summary[,12] <- apply(apply(est_arr[, 11, ], 2, function(x) x<para_true)&apply(est_arr[, 12, ], 2, function(x) x>para_true), 1, mean)


est_summary3 <- est_summary

xtable(cbind(est_summary1), digits = 3)
xtable(cbind(est_summary2), digits = 3)
xtable(cbind(est_summary3), digits = 3)
