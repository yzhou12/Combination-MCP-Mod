getwd()
rm(list=ls())
library(multcomp)
library(dplyr)
library(ggplot2)
############### combine multiple jobs results ##################
all.model <- c("flat", "linear", "linlog", "quadratic", "exponential","emax", "sigmoid")
dir.create("../data/result-additive")
############### Flat ##################
f <- function(x){mean(x, na.rm = T)}
load(paste0("./k9/additive/output_n5/out_flat.RData"))
pr.n1 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100, pr.dose = apply(out, 2, f)[3:4] * 100, 
                    method = c("MultiComp", "2D MCPMod"), n = 5)
load(paste0("./k9/additive/output_n10/out_flat.RData"))
pr.n2 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100, pr.dose = apply(out, 2, f)[3:4] * 100, 
                    method = c("MultiComp", "2D MCPMod"), n = 10)
load(paste0("./k9/additive/output_n15/out_flat.RData"))
pr.n3 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100, pr.dose = apply(out, 2, f)[3:4] * 100, 
                    method = c("MultiComp", "2D MCPMod"), n = 15)
load(paste0("./k9/additive/output_n20/out_flat.RData"))
pr.n4 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100, pr.dose = apply(out, 2, f)[3:4] * 100, 
                    method = c("MultiComp", "2D MCPMod"), n = 20)
data.flat.1 <- rbind(pr.n1, pr.n2, pr.n3, pr.n4)
data.flat.1$dose <- "9 dose levels"
load(paste0("./k16/additive/output_n5/out_flat.RData"))
pr.n1 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100, pr.dose = apply(out, 2, f)[3:4] * 100, 
                    method = c("MultiComp", "2D MCPMod"), n = 5)
load(paste0("./k16/additive/output_n10/out_flat.RData"))
pr.n2 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100, pr.dose = apply(out, 2, f)[3:4] * 100, 
                    method = c("MultiComp", "2D MCPMod"), n = 10)
load(paste0("./k16/additive/output_n15/out_flat.RData"))
pr.n3 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100, pr.dose = apply(out, 2, f)[3:4] * 100, 
                    method = c("MultiComp", "2D MCPMod"), n = 15)
load(paste0("./k16/additive/output_n20/out_flat.RData"))
pr.n4 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100, pr.dose = apply(out, 2, f)[3:4] * 100, 
                    method = c("MultiComp", "2D MCPMod"), n = 20)
data.flat.2 <- rbind(pr.n1, pr.n2, pr.n3, pr.n4)
data.flat.2$dose <- "16 dose levels"
data.flat <- rbind(data.flat.1, data.flat.2)
#################################
inter.effect <- "additive"
#inter.effect <- "synergistic"
#inter.effect <- "antagonistic"
if(inter.effect == "additive"){
  all.model <- c("linear", "linlog", "quadratic", "exponential","emax", "sigmoid")}
if(inter.effect == "synergistic"){
  all.model <- c("linear", "linlog", "quadratic", "exponential","emax", "sigmoid")}
if(inter.effect == "antagonistic"){
  all.model <- c("linear", "linlog", "quadratic", "exponential","emax", "sigmoid")}

model.list <- all.model
alpha <- 0.05
# clinical relevant effect Delta = 1.2 (here we use -y as the response to fit the MCP-Mod):
Delta <- 1.2
# standard deviation of random error in DR model:
dose.d1 <- c(0, 0.2, 0.6, 1)
dose.d2 <- c(0, 0.2, 0.6, 1)
d1 <- rep(dose.d1, each = length(dose.d1))
d2 <- rep(dose.d2, length(dose.d2))
doses <- data.frame(d1 = d1, d2 = d2)
# dose levels: k
doses$level <- 1:nrow(doses)
k <- max(doses$level)

data <- NULL
tab.mod <- NULL
for(i in 1:6){
  model <- model.list[i]
  source("../source/true_model.R")
  # true response mu
  mu <- apply(doses[,1:2], 1, gen_arm)
  f2 <- function(x, prob){quantile(x, probs = prob, na.rm = T)}
  ################# Pr(DR), Pr(TD) #################
  load(paste0("./k9/", inter.effect, "/output_n5/out_", model, ".RData"))
  prmod1 <- c(sum(out[,5] == which(all.model[1:6] == model), na.rm = T)/nrow(out), 
              sum(out[,6] == which(all.model[1:6] == model), na.rm = T)/nrow(out), n = 5, dose = 9)
  pr.n1 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100, pr.dose = apply(out, 2, f)[3:4] * 100, 
                      method = c("MultiComp", "2D MCPMod"), n = 5, model = model, dose = "9 dose levels")
  
  load(paste0("./k9/", inter.effect, "/output_n10/out_", model, ".RData"))
  prmod2 <- c(sum(out[,5] == which(all.model[1:6] == model), na.rm = T)/nrow(out), 
              sum(out[,6] == which(all.model[1:6] == model), na.rm = T)/nrow(out), n = 10, dose = 9)
  pr.n2 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100, pr.dose = apply(out, 2, f)[3:4] * 100, 
                      method = c("MultiComp", "2D MCPMod"), n = 10, model = model, dose = "9 dose levels")
  
  load(paste0("./k9/", inter.effect, "/output_n15/out_", model, ".RData"))
  prmod3 <- c(sum(out[,5] == which(all.model[1:6] == model), na.rm = T)/nrow(out), 
              sum(out[,6] == which(all.model[1:6] == model), na.rm = T)/nrow(out), n = 15, dose = 9)
  pr.n3 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100, pr.dose = apply(out, 2, f)[3:4] * 100, 
                      method = c("MultiComp", "2D MCPMod"), n = 15, model = model, dose = "9 dose levels")
  
  load(paste0("./k9/", inter.effect, "/output_n20/out_", model, ".RData"))
  prmod4 <- c(sum(out[,5] == which(all.model[1:6] == model), na.rm = T)/nrow(out), 
              sum(out[,6] == which(all.model[1:6] == model), na.rm = T)/nrow(out), n = 20, dose = 9)
  pr.n4 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100, pr.dose = apply(out, 2, f)[3:4] * 100, 
                      method = c("MultiComp", "2D MCPMod"), n = 20, model = model, dose = "9 dose levels")
  
  load(paste0("./k16/", inter.effect, "/output_n5/out_", model, ".RData"))
  prmod5 <- c(sum(out[,5] == which(all.model[1:6] == model), na.rm = T)/nrow(out), 
              sum(out[,6] == which(all.model[1:6] == model), na.rm = T)/nrow(out), n = 5, dose = 16)
  pr.n5 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100, pr.dose = apply(out, 2, f)[3:4] * 100, 
                      method = c("MultiComp", "2D MCPMod"), n = 5, model = model, dose = "16 dose levels")
  
  load(paste0("./k16/", inter.effect, "/output_n10/out_", model, ".RData"))
  prmod6 <- c(sum(out[,5] == which(all.model[1:6] == model), na.rm = T)/nrow(out), 
              sum(out[,6] == which(all.model[1:6] == model), na.rm = T)/nrow(out), n = 10, dose = 16)
  pr.n6 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100, pr.dose = apply(out, 2, f)[3:4] * 100, 
                      method = c("MultiComp", "2D MCPMod"), n = 10, model = model, dose = "16 dose levels")
  load(paste0("./k16/", inter.effect, "/output_n15/out_", model, ".RData"))
  prmod7 <- c(sum(out[,5] == which(all.model[1:6] == model), na.rm = T)/nrow(out), 
              sum(out[,6] == which(all.model[1:6] == model), na.rm = T)/nrow(out), n = 15, dose = 16)
  pr.n7 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100, pr.dose = apply(out, 2, f)[3:4] * 100, 
                      method = c("MultiComp", "2D MCPMod"), n = 15, model = model, dose = "16 dose levels")
  
  load(paste0("./k16/", inter.effect, "/output_n20/out_", model, ".RData"))
  prmod8 <- c(sum(out[,5] == which(all.model[1:6] == model), na.rm = T)/nrow(out), 
              sum(out[,6] == which(all.model[1:6] == model), na.rm = T)/nrow(out), n = 20, dose = 16)
  pr.n8 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100, pr.dose = apply(out, 2, f)[3:4] * 100, 
                      method = c("MultiComp", "2D MCPMod"), n = 20, model = model, dose = "16 dose levels")
  
  data <- rbind(data, pr.n1, pr.n2, pr.n3, pr.n4, pr.n5, pr.n6, pr.n7, pr.n8)
  mod <- as.data.frame(rbind(prmod1, prmod2, prmod3, prmod4, prmod5, prmod6, prmod7, prmod8))
  mod$model <- model
  tab.mod <- rbind(tab.mod, mod)
}

save(data.flat, data, tab.mod, file="result-additive/all_result.RData")

################# DO NOT RUN #################

### Table of PoC power:
data.dr <- NULL
data.do <- NULL
inter.effect <- "additive"
all.model2 <- c("flat", all.model)
for(model in all.model2){
  ################# Pr(DR), Pr(TD) #################
  load(paste0("./k9/", inter.effect, "/output_n5/out_", model, ".RData"))
  pr.n1 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100,
                      method = c("MultiComp", "2D MCPMod"), n = 5, model = model, dose = "9 dose levels")
  do.n1 <- data.frame(pr.dose = apply(out, 2, f)[3:4] * 100, 
                      method = c("MultiComp", "2D MCPMod"), n = 5, model = model, dose = "9 dose levels")
  load(paste0("./k9/", inter.effect, "/output_n10/out_", model, ".RData"))
  pr.n2 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100,
                      method = c("MultiComp", "2D MCPMod"), n = 10, model = model, dose = "9 dose levels")
  do.n2 <- data.frame(pr.dose = apply(out, 2, f)[3:4] * 100, 
                      method = c("MultiComp", "2D MCPMod"), n = 10, model = model, dose = "9 dose levels")
  
  load(paste0("./k9/", inter.effect, "/output_n15/out_", model, ".RData"))
  pr.n3 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100,
                      method = c("MultiComp", "2D MCPMod"), n = 15, model = model, dose = "9 dose levels")
  do.n3 <- data.frame(pr.dose = apply(out, 2, f)[3:4] * 100, 
                      method = c("MultiComp", "2D MCPMod"), n = 15, model = model, dose = "9 dose levels")
  
  load(paste0("./k9/", inter.effect, "/output_n20/out_", model, ".RData"))
  pr.n4 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100,
                      method = c("MultiComp", "2D MCPMod"), n = 20, model = model, dose = "9 dose levels")
  do.n4 <- data.frame(pr.dose = apply(out, 2, f)[3:4] * 100, 
                      method = c("MultiComp", "2D MCPMod"), n = 20, model = model, dose = "9 dose levels")
  
  load(paste0("./k16/", inter.effect, "/output_n5/out_", model, ".RData"))
  pr.n5 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100, 
                      method = c("MultiComp", "2D MCPMod"), n = 5, model = model, dose = "16 dose levels")
  do.n5 <- data.frame(pr.dose = apply(out, 2, f)[3:4] * 100, 
                      method = c("MultiComp", "2D MCPMod"), n = 5, model = model, dose = "16 dose levels")
  
  load(paste0("./k16/", inter.effect, "/output_n10/out_", model, ".RData"))
  pr.n6 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100,
                      method = c("MultiComp", "2D MCPMod"), n = 10, model = model, dose = "16 dose levels")
  do.n6 <- data.frame(pr.dose = apply(out, 2, f)[3:4] * 100, 
                      method = c("MultiComp", "2D MCPMod"), n = 10, model = model, dose = "16 dose levels")
  
  load(paste0("./k16/", inter.effect, "/output_n15/out_", model, ".RData"))
  pr.n7 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100,
                      method = c("MultiComp", "2D MCPMod"), n = 15, model = model, dose = "16 dose levels")
  do.n7 <- data.frame(pr.dose = apply(out, 2, f)[3:4] * 100, 
                      method = c("MultiComp", "2D MCPMod"), n = 15, model = model, dose = "16 dose levels")
  
  load(paste0("./k16/", inter.effect, "/output_n20/out_", model, ".RData"))
  pr.n8 <- data.frame(pr.dr = apply(out, 2, f)[1:2] * 100,
                      method = c("MultiComp", "2D MCPMod"), n = 20, model = model, dose = "16 dose levels")
  do.n8 <- data.frame(pr.dose = apply(out, 2, f)[3:4] * 100, 
                      method = c("MultiComp", "2D MCPMod"), n = 20, model = model, dose = "16 dose levels")
  
  data.dr <- rbind(data.dr, pr.n1, pr.n2, pr.n3, pr.n4, pr.n5, pr.n6, pr.n7, pr.n8)
  data.do <- rbind(data.do, do.n1, do.n2, do.n3, do.n4, do.n5, do.n6, do.n7, do.n8)
}
data.dr

data.d9.aov <- data.dr[data.dr$dose == "9 dose levels" & data.dr$method == "MultiComp",]
data.d9.mcp <- data.dr[data.dr$dose == "9 dose levels" & data.dr$method == "2D MCPMod",]
data.d16.aov <- data.dr[data.dr$dose == "16 dose levels" & data.dr$method == "MultiComp",]
data.d16.mcp <- data.dr[data.dr$dose == "16 dose levels" & data.dr$method == "2D MCPMod",]
a <- data.frame(model = data.d9.aov[data.d9.mcp$n == 5, 4],
                data.d9.mcp[data.d9.mcp$n == 5, 1], 
                data.d9.mcp[data.d9.mcp$n == 10, 1], 
                data.d9.mcp[data.d9.mcp$n == 15, 1], 
                data.d9.mcp[data.d9.mcp$n == 20, 1], 
                data.d9.aov[data.d9.aov$n == 5, 1], 
                data.d9.aov[data.d9.aov$n == 10, 1], 
                data.d9.aov[data.d9.aov$n == 15, 1],
                data.d9.aov[data.d9.aov$n == 20, 1])
colnames(a) <- c("model k9",  "n5mcp", "n10mcp", "n15mcp", "n20mcp", "n5aov", "n10aov", "n15aov", "n20aov")
b <- data.frame(model = data.d16.aov[data.d16.mcp$n == 5, 4],
                data.d16.mcp[data.d16.mcp$n == 5, 1], 
                data.d16.mcp[data.d16.mcp$n == 10, 1], 
                data.d16.mcp[data.d16.mcp$n == 15, 1], 
                data.d16.mcp[data.d16.mcp$n == 20, 1], 
                data.d16.aov[data.d16.aov$n == 5, 1], 
                data.d16.aov[data.d16.aov$n == 10, 1], 
                data.d16.aov[data.d16.aov$n == 15, 1],
                data.d16.aov[data.d16.aov$n == 20, 1])
colnames(b) <- c("model k16",  "n5mcp", "n10mcp", "n15mcp", "n20mcp", "n5aov", "n10aov", "n15aov", "n20aov")

write.csv(rbind(t(a), t(b)), file = paste0("./result-", inter.effect,"/prob-DR.csv"), quote = F)


### plot of bias of estimated effect \hat\mu_ij
all.model <- c("flat", "linear", "linlog", "quadratic", "exponential","emax", "sigmoid")
inter.effect <- 'additive'
for(k in c(9, 16)){
  if(k==9){
    dose.d1 <- c(0, 0.6, 1)
    dose.d2 <- c(0, 0.6, 1)
  }
  if(k==16){
    dose.d1 <- c(0, 0.2, 0.6, 1)
    dose.d2 <- c(0, 0.2, 0.6, 1)
  }
  d1 <- rep(dose.d1, each = length(dose.d1))
  d2 <- rep(dose.d2, length(dose.d2))
  doses <- data.frame(d1 = d1, d2 = d2)
  for(model in all.model[1:7]){
    ########### pre-specifiy##########
    f <- function(x){mean(x, na.rm = T)}
    source("../source/true_model.R")
    # true response mu
    mu <- apply(doses[,1:2], 1, gen_arm)
    f2 <- function(x, prob){quantile(x, probs = prob, na.rm = T)}
    ################# quantiles of prediction error #################
    #par(mfrow = c(1,4))
    ylimit <- c(-1, 1)
    jpeg(paste0("../data/result-", inter.effect,"/bias_mu_k", k, "_", model, ".jpeg"), width = 10, height = 4, units = "in", res = 500)
    par(mfrow = c(1, 4))
    for(n in c(5, 10, 15, 20)){
      load(paste0("../data/k", k, "/", inter.effect, "/output_n", n, "/out_", model, ".RData"))
      plot(apply(out[,7:(6+k)], 2, f2, prob = 0.5) - mu, type = "none",  ylim = ylimit,
           main = paste0("Comb MCP-Mod (n=", n, ")"), #ylab = "Prediction error",
           ylab = expression(hat(mu)-mu), xlab = "Dose level") # mcp
      lines(apply(out[,7:(6+k)], 2, f2, prob = 0.5) - mu)
      lines(apply(out[,7:(6+k)], 2, f2, prob = 0.25) - mu, lty = 2)
      lines(apply(out[,7:(6+k)], 2, f2, prob = 0.75) - mu, lty = 2)
      lines(apply(out[,7:(6+k)], 2, f2, prob = 0.05) - mu, lty = 3)
      lines(apply(out[,7:(6+k)], 2, f2, prob = 0.95) - mu, lty = 3)
      abline(h = 0, col = 2, lty = 4)
    }
    #plot(apply(out[,(7+k):(7+k+k-1)], 2, f2, prob = 0.5) - mu, ylim = ylimit, type = "none", 
    #     main = "ANOVA (n=10)", ylab = "Prediction error", xlab = "Dose level") # mcp
    #lines(apply(out[,(7+k):(7+k+k-1)], 2, f2, prob = 0.5) - mu)
    #lines(apply(out[,(7+k):(7+k+k-1)], 2, f2, prob = 0.25) - mu, lty = 2)
    #lines(apply(out[,(7+k):(7+k+k-1)], 2, f2, prob = 0.75) - mu, lty = 2)
    #lines(apply(out[,(7+k):(7+k+k-1)], 2, f2, prob = 0.05) - mu, lty = 3)
    #lines(apply(out[,(7+k):(7+k+k-1)], 2, f2, prob = 0.95) - mu, lty = 3)
    #abline(h = 0, col = 2, lty = 4)
    legend("topright", legend = c("50%", "25%, 75%", "5%, 95%"), lty = 1:3)
    dev.off()
    ##################################
  }
  
}

################# DO NOT RUN #################

