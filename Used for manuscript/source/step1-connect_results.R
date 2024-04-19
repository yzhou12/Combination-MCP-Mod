rm(list=ls())
library(multcomp)
library(dplyr)
library(ggplot2)
setwd("../data")
############### combine multiple jobs results ##################
all.model <- c("flat", "linear", "linlog", "quadratic", "exponential","emax", "sigmoid")
for(k in c(9,16)){
  for(n in c(5,10,15,20)){
    for(model in all.model){
      all <- NULL
      coef.mcp <- coef.aov <- NULL
      for(i in 1:10){
        load(paste0("k",k,"/additive/output_n",n,"/out_", model, "_", i, ".RData"))
        all <- rbind(all, out.B)
        ### get parameter estimate
        coef.mcp <- rbind(coef.mcp, out.coef.mcp)
        coef.aov <- rbind(coef.aov, out.coef.aov)
      }
      colnames(all) <- c("detect.DR.aov", "detect.DR.mcp", "detect.dose.aov", "detect.dose.mcp", 
                         "select.model.mcp", "select.model.aov", 
                         paste0("mu", 1:k, ".mcp"),
                         paste0("mu", 1:k, ".aov"),
                         paste0("T_", can.model.set),
                         "used.model.mcp", "used.model.aov")
      out <- all
      save(doses, can.model.set, out, coef.mcp, coef.aov,
           file = paste0("k",k,"/additive/output_n",n,"/out_", model, ".RData"))
    }
    print(n)
  }
}

