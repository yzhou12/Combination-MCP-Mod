getwd()
#"C:/Users/ZHOUY195/OneDrive - Pfizer/Desktop/Copy from C/Research/Combination/Manuscript_Simulation/simulation-revision-Jun2023/simu_14Jun2023_largebias/data"
rm(list=ls()) 
library(ggplot2)
library(dplyr)
library(reshape2)

load("result-additive/all_result.RData")
data.flat$method <- ifelse(data.flat$method=="2D MCPMod", "Comb MCP-Mod", "Multiple Comparisons")
data$method <- ifelse(data$method=="2D MCPMod", "Comb MCP-Mod", "Multiple Comparisons")
data$dose<-factor(data$dose, levels=c("9 dose levels", "16 dose levels"))
data$model <- factor(data$model, levels=c("linear", "linlog", "quadratic", "exponential","emax", "sigmoid"))
#flat
jpeg("./result-additive/flat-dr.jpeg", width = 5, height = 2, units = "in", res = 500)
ggplot(data = data.flat, aes(x = method, y = pr.dr, fill = dose)) + 
  geom_bar(stat = "identity",width = 0.7, color="black", position=position_dodge()) + 
  ylim(0, 6) + scale_fill_manual(values=c('#999999','#E69F00')) + 
  geom_hline(yintercept = 5, col = "red") + 
  ylab("Pr(DR|flat) (%)") + xlab("") +
  coord_flip() +
  facet_wrap( ~ n, ncol = 4, labeller = label_bquote("n = " * .(n))) + 
  theme_bw()+
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8,face="bold"),
        text = element_text(size=8))
dev.off()
jpeg("./result-additive/flat-dose.jpeg", width = 5, height = 2, units = "in", res = 500)
ggplot(data = data.flat, aes(x = method, y = pr.dose, fill = dose)) + 
  geom_bar(stat = "identity",width = 0.7, color="black", position=position_dodge()) + 
  ylim(0, 6) + scale_fill_manual(values=c('#999999','#E69F00')) + 
  ylab("Pr(TD|flat) (%)") + xlab("") +
  geom_hline(yintercept = 5, col = "red") + 
  coord_flip() +
  facet_wrap( ~ n, ncol = 4, labeller = label_bquote("n = " * .(n))) + 
  theme_bw()+
  theme(legend.position = "bottom")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=8,face="bold"),
        text = element_text(size=8))
dev.off()
#active
inter.effect <- "additive"
jpeg(paste0("./result-", inter.effect,"/all-dr-", inter.effect,"-step.jpeg"), width = 10, height = 3, units = "in", res = 500)
ggplot(data = data, aes(x = method, y = pr.dr)) + 
  geom_point(aes(shape=dose)) + ylab("Pr(DR) (%)") + xlab("") +
  #ylim(0,100) + 
  scale_shape_manual(values=c(4, 0)) +
  coord_flip() +
  facet_grid(cols = vars(model), rows = vars(n), labeller = label_bquote("n = " * .(n)))  + 
  #scale_y_continuous(breaks = seq(90, 100, by = 5)) +
  theme_bw()+
  theme(legend.position = "bottom")
dev.off()
#table of pr(DR)
dat <- data %>% select(-pr.dose)
tab.dr <- dcast(dat, dose+n+method~model, value.var = 'pr.dr')

jpeg(paste0("./result-", inter.effect,"/all-dose-", inter.effect,"-step.jpeg"), width = 10, height = 3, units = "in", res = 500)
ggplot(data = data, aes(x = method, y = pr.dose)) + 
  geom_point(aes(shape=dose)) +  ylab("Pr(TD) (%)") + xlab("") +
  #ylim(0,100) +
  scale_shape_manual(values=c(4, 0)) +
  coord_flip() +
  facet_grid(cols = vars(model), rows = vars(n), labeller = label_bquote("n = " * .(n))) + 
  #scale_y_continuous(breaks = seq(70, 100, by = 5)) +
  theme_bw()+
  theme(legend.position = "bottom")
dev.off()
#table of pr(dose)
dat <- data %>% select(-pr.dr)
tab.dose <- dcast(dat, dose+n+method~model, value.var = 'pr.dose')
#save tables
write.csv(tab.dr, file = "./result-additive/Pr_DR.csv", row.names = F, col.names = F)
write.csv(tab.dose, file = "./result-additive/Pr_TD.csv", row.names = F, col.names = F)

#table pr(identify correct model)
inter.effect <- "additive"
tab.mod$model <- factor(tab.mod$model, levels=c("linear", "linlog", "quadratic", "exponential","emax", "sigmoid"))
tab.mod[,1:2] <- tab.mod[,1:2] * 100
colnames(tab.mod)[1:2] <- c("MCPMod", "MultiComp")
tab.mod <- tab.mod %>% select(-'MultiComp') %>% filter(model!='step')
tab.mod$MCPMod <- tab.mod$MCPMod
res.all <- dcast(tab.mod, n+dose ~ model, value.var = 'MCPMod')
tab9 <- res.all %>% filter(dose==9)
tab16 <- res.all %>% filter(dose==16)

pr.mod.res <- rbind(tab9, tab16)

write.csv(pr.mod.res, file = paste0("./result-", inter.effect,"/detect-mod_update.csv"), quote = F, row.names = F)


### plot of bias of estimated effect \hat\mu_ij
all.model <- c("linear", "linlog", "quadratic", "exponential","emax", "sigmoid")
inter.effect <- 'additive'
for(model in all.model[1:6]){
  ########### pre-specifiy##########
  f <- function(x){mean(x, na.rm = T)}
  source("../source/true_model.R")
  f2 <- function(x, prob){quantile(x, probs = prob, na.rm = T)}
  ################# quantiles of prediction error #################
  ylimit <- c(-1, 1)
  jpeg(paste0("../data/result-", inter.effect,"/bias_mu_", model, ".jpeg"), 
       width = 8, height = 6, units = "in", res = 500)
  par(mfrow = c(2, 4))
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
    # true response mu
    mu <- apply(doses[,1:2], 1, gen_arm)
    
    for(n in c(5, 10, 15, 20)){
      load(paste0("../data/k", k, "/", inter.effect, "/output_n", n, "/out_", model, ".RData"))
      plot(apply(out[,7:(6+k)], 2, f2, prob = 0.5) - mu, type = "none",  ylim = ylimit,
           main = paste0("n=", n, ", k=", k), #ylab = "Prediction error",
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
    ##################################
  }
  legend("topright", legend = c("50%", "25%, 75%", "5%, 95%"), lty = 1:3)
  dev.off()
  
}
### table of bias (as mean of muhat-mu)
all.model <- c("flat", "linear", "linlog", "quadratic", "exponential","emax", "sigmoid")
inter.effect <- 'additive'
res.bias.k9 <- res.bias.k16 <- NULL
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
  for(n in c(5, 10, 15, 20)){
    ########### pre-specifiy##########
    f <- function(x){mean(x, na.rm = T)}
    source("../source/true_model.R")
    # true response mu
    mu <- apply(doses[,1:2], 1, gen_arm)
    f2 <- function(x, prob){quantile(x, probs = prob, na.rm = T)}
    ################# quantiles of prediction error #################
    #par(mfrow = c(1,4))
    ylimit <- c(-1, 1)
    for(model in all.model[1:7]){
      load(paste0("../data/k", k, "/", inter.effect, "/output_n", n, "/out_", model, ".RData"))
      bias.mean.mcp <- as.data.frame(t(apply(out[,7:(6+k)], 2, f) - mu))
      #bias.mean.aov <- mean(apply(out[,(7+k):(7+k+k-1)], 2, f) - mu)
      res1 <- data.frame(dose=k, n=n, model = model)
      res <- cbind(res1, bias.mean.mcp)
      names(res) <- c("dose", "n", "model", paste0("(", paste(doses$d1, doses$d2, sep = ", "), ")"))
      if(k==9) res.bias.k9 <- rbind(res.bias.k9, res)
      if(k==16) res.bias.k16 <- rbind(res.bias.k16, res)
    }
    ##################################
  }
}
write.csv(res.bias.k9, file = paste0("./result-", inter.effect,"/bias_mu_mean_k9.csv"), quote = F, row.names = F)
write.csv(res.bias.k16, file = paste0("./result-", inter.effect,"/bias_mu_mean_k16.csv"), quote = F, row.names = F)


############# parameter estimate
all.model <- c("flat", "linear", "linlog", "quadratic", "exponential","emax", "sigmoid")
inter.effect <- 'additive'

model <- 'linear'
coef.true <- c(0, 0.75, 0.75)
ylimit <- c(-1,1.2)
source("fun_step3-plot.R")

model <- 'linlog'
coef.true <- c(0, 1.082, 1.082)
ylimit <- c(-1,1.6)
source("fun_step3-plot.R")

model <- 'quadratic'
coef.true <- c(0, 2.508, 2.508, - 2.09, - 2.09)
ylimit <- c(-4,4)
source("fun_step3-plot.R")

model <- 'exponential'
coef.true <- c(0, 0.067, 0.067, 0.4, 0.4)
ylimit <- c(-0.6,0.6)
source("fun_step3-plot.R")

model <- 'emax'
coef.true <- c(0, 0.9, 0.9, 0.2, 0.2)
ylimit <- c(-1,4)
source("fun_step3-plot.R")

model <- 'sigmoid'
coef.true <- c(0, 0.87, 0.87, 0.4^2, 0.4^2)
ylimit <- c(-1,4)
source("fun_step3-plot.R")




