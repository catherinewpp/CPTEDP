rm(list=ls())
library(ggplot2)
library(scales)
library(grid)
library(xtable)
library(rlang)
library(plotly)
source("D:/My Document/Rcode/xtab.R")
source("D:/My Document/bugscode/Change-point-Codes/R2OpenBUGS/led/read_led.R")

est_gp_in <- read.table("D:/My Document/bugscode/Change-point-Codes/Codes/led_data/gp_con_inde.txt")

result_e <- result_sd <- matrix(NA, 6, 3) 
for(i in 1:6){
  result_e[i, ] <- est_gp_in[((i-1)*3+1:3), 1]
  
  result_sd[i, ] <- est_gp_in[((i-1)*3+1:3), 2]
}

result_e[,1] <- -result_e[,1]
result_e[,2] <- -result_e[,2]

xtable(result_e, digits = 3)
xtable(result_sd, digits = 3)
xtab(result_e, result_sd, digits = 3)


est_ed <- bugs.log("D:/My Document/Rcode/degradation-ED/ed_change_point/ed_dep/log.txt")$stats

result_e <- result_sd <- matrix(NA, 6, 3) 
for(i in 1:6){
  result_e[i, ] <- est_ed[((i-1)*3+9+1:3), 1]
  
  result_sd[i, ] <- est_ed[((i-1)*3+9+1:3), 2]
}

result_e[,1] <- result_e[,1]
result_e[,2] <- result_e[,2]

xtable(result_e, digits = 3)
xtable(result_sd, digits = 3)
xtab(result_e, result_sd, digits = 3)


est_two <- read.table("D:/My Document/bugscode/Change-point-Codes/Codes/led_data/beta_two.txt")

result_e <- result_sd <- matrix(NA, 6, 3) 
for(i in 1:6){
  result_e[i, ] <- est_two[((i-1)*3+9+1:3), 1]
  
  result_sd[i, ] <- est_two[((i-1)*3+9+1:3), 2]
}

xtable(result_e, digits = 3)
xtable(result_sd, digits = 3)
xtab(result_e, result_sd, digits = 3)

est_bi <- read.table("D:/My Document/bugscode/Change-point-Codes/Codes/led_data/gamma_bi.txt")
result_e <- result_sd <- matrix(NA, 6, 3) 
for(i in 1:6){
  result_e[i, ] <- est_bi[((i-1)*3+1:3), 1]
  
  result_sd[i, ] <- est_bi[((i-1)*3+1:3), 2]
}

result_e[,1] <- -result_e[,1]
result_e[,2] <- -result_e[,2]

xtable(result_e, digits = 3)
xtable(result_sd, digits = 3)
xtab(result_e, result_sd, digits = 3)

