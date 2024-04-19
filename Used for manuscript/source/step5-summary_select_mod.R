rm(list=ls())
getwd()
library(dplyr)
options(digits = 4)
setwd('../data')
################################################################
### get the distribution of selected model for each scenario ###
load(paste0("k",9,"/additive/output_n",5,"/out_emax.RData"))
add.sum.select.mod <- NULL
for(model in can.model.set){
  sum.all <- NULL
  for(n in c(5,10,15,20)){
    k=9
    load(paste0("k",k,"/additive/output_n",n,"/out_", model, ".RData"))
    out <- as.data.frame(out)
    sum(out$detect.dose.mcp)
    dat.canmod <- data.frame(select.model.mcp = 1:6, sel.model.mcp = can.model.set)
    out1 <- left_join(out, dat.canmod) 
    total <- sum(is.na(out1$select.model.mcp)==0)
    sum9 <- c(sum(out1$sel.model.mcp%in%can.model.set[1]),
              sum(out1$sel.model.mcp%in%can.model.set[2]),
              sum(out1$sel.model.mcp%in%can.model.set[3]),
              sum(out1$sel.model.mcp%in%can.model.set[4]),
              sum(out1$sel.model.mcp%in%can.model.set[5]),
              sum(out1$sel.model.mcp%in%can.model.set[6]))/total*100

    k=16
    load(paste0("k",k,"/additive/output_n",n,"/out_", model, ".RData"))
    out <- as.data.frame(out)
    dat.canmod <- data.frame(select.model.mcp = 1:6, sel.model.mcp = can.model.set)
    out1 <- left_join(out, dat.canmod) 
    total <- sum(is.na(out1$select.model.mcp)==0)
    sum16 <- c(sum(out1$sel.model.mcp%in%can.model.set[1]),
              sum(out1$sel.model.mcp%in%can.model.set[2]),
              sum(out1$sel.model.mcp%in%can.model.set[3]),
              sum(out1$sel.model.mcp%in%can.model.set[4]),
              sum(out1$sel.model.mcp%in%can.model.set[5]),
              sum(out1$sel.model.mcp%in%can.model.set[6]))/total*100    
    sum.all <- rbind(sum.all, rbind(c(n, 9, sum9), c(n, 16, sum16)))
  }
  sum.all <- as.data.frame(sum.all)
  sum.all$True_model <- model
  add.sum.select.mod <- rbind(add.sum.select.mod, sum.all)
}
names(add.sum.select.mod) <- c("n", "k", can.model.set, "TrueModel")
write.csv(add.sum.select.mod, file = paste0("./result-additive/dist_selected_mod_update.csv"), quote = F, row.names = F)
##########################END###################################

rm(list=ls())
################################################################
### get the distribution of used model for each scenario ###
load(paste0("k",9,"/additive/output_n",5,"/out_emax.RData"))
add.sum.select.mod <- NULL
for(model in can.model.set){
  sum.all <- NULL
  for(n in c(5,10,15,20)){
    k=9
    load(paste0("k",k,"/additive/output_n",n,"/out_", model, ".RData"))
    out <- as.data.frame(out)
    sum(out$detect.dose.mcp)
    dat.canmod <- data.frame(used.model.mcp = 1:6, use.model.mcp = can.model.set)
    out1 <- left_join(out, dat.canmod) 
    total <- sum(is.na(out1$used.model.mcp)==0)
    sum9 <- c(sum(out1$use.model.mcp%in%can.model.set[1]),
              sum(out1$use.model.mcp%in%can.model.set[2]),
              sum(out1$use.model.mcp%in%can.model.set[3]),
              sum(out1$use.model.mcp%in%can.model.set[4]),
              sum(out1$use.model.mcp%in%can.model.set[5]),
              sum(out1$use.model.mcp%in%can.model.set[6]))/total*100
    
    k=16
    load(paste0("k",k,"/additive/output_n",n,"/out_", model, ".RData"))
    out <- as.data.frame(out)
    dat.canmod <- data.frame(used.model.mcp = 1:6, use.model.mcp = can.model.set)
    out1 <- left_join(out, dat.canmod) 
    total <- sum(is.na(out1$used.model.mcp)==0)
    sum16 <- c(sum(out1$use.model.mcp%in%can.model.set[1]),
               sum(out1$use.model.mcp%in%can.model.set[2]),
               sum(out1$use.model.mcp%in%can.model.set[3]),
               sum(out1$use.model.mcp%in%can.model.set[4]),
               sum(out1$use.model.mcp%in%can.model.set[5]),
               sum(out1$use.model.mcp%in%can.model.set[6]))/total*100    
    sum.all <- rbind(sum.all, rbind(c(n, 9, sum9), c(n, 16, sum16)))
  }
  sum.all <- as.data.frame(sum.all)
  sum.all$True_model <- model
  add.sum.select.mod <- rbind(add.sum.select.mod, sum.all)
}
names(add.sum.select.mod) <- c("n", "k", can.model.set, "TrueModel")
write.csv(add.sum.select.mod, file = paste0("./result-additive/dist_used_mod_update.csv"), quote = F, row.names = F)
##########################END###################################
