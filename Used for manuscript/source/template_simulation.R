jobId <- as.numeric(commandArgs(TRUE))
library(multcomp)
library(dplyr)
################## Simulation setting ################## 
# set simulation times as 10,000
B <- 1000
# Possible dose for drug A and drug B
dose.d1 <- c(0, 0.6, 1)
dose.d2 <- c(0, 0.6, 1)
d1 <- rep(dose.d1, each = length(dose.d1))
d2 <- rep(dose.d2, length(dose.d2))
doses <- data.frame(d1 = d1, d2 = d2)
# set underlying model
dr.model <- c( "flat", "linear", "linlog", "quadratic", "exponential","emax","sigmoid",
               "step")
model <- dr.model[7]
n <- 5
### response generating models (true models) are specified here
if(model == "step"){
  gen_arm <- function(x){
    ifelse(x[1] + x[2] <= 0.8, 0, 1.5 ) + rnorm(n = 1, mean = 0, sd = sigma)}
}
if(model == "flat"){
  gen_arm <- function(x){rnorm(n = 1, mean = 0, sd = sigma)}}
### Additive effect:
if(model == "linear"){
  gen_arm <- function(x){
    0.75 * x[1] + 0.75 * x[2] + rnorm(n = 1, mean = 0, sd = sigma)}
}
if(model == "linlog"){
  gen_arm <- function(x){
    1.082 * log(x[1] + 1) + 1.082 * log(x[2] + 1) + rnorm(n = 1, mean = 0, sd = sigma)}
}
if(model == "quadratic"){
  gen_arm <- function(x){
    2.09*2*0.6*x[1] - 2.09*x[1]^2 + 2.09*2*0.6*x[2] - 2.09*x[2]^2 + rnorm(n = 1, mean = 0, sd = sigma)}
}
if(model == "exponential"){
  gen_arm <- function(x){
    0.067*(exp(x[1]/0.4) - 1) + 0.067*(exp(x[2]/0.4) - 1) + rnorm(n = 1, mean = 0, sd = sigma)}
}
if(model == "emax"){
  gen_arm <- function(x){
    0.9*x[1] /(0.2 + x[1]) + 0.9*x[2]/(0.2 + x[2]) + rnorm(n = 1, mean = 0, sd = sigma)}
}
if(model == "sigmoid"){
  gen_arm <- function(x){0.87*x[1]^2 /(0.4^2 + x[1]^2) + 
      0.87*x[2]^2/(0.4^2 + x[2]^2)+ rnorm(n = 1, mean = 0, sd = sigma)}
}
# dose levels: k
doses$level <- 1:nrow(doses)
k <- max(doses$level)
dir.create(paste0("../data/k", k)) # in case there is no output directory
setwd(paste0("../data/k", k))
stat.comp.level <- doses$level[which(doses$d1 + doses$d2 == 0 | doses$d1 * doses$d2 != 0)] # levels used to conduct statistical comparison in ANOVA
################## Assumptions ################## 
# significane level (FWER)
alpha <- 0.05
# clinical relevant effect Delta = 1.2
Delta <- 1.2
# standard deviation of response in DR model:
sigma <- 1.4
################## Prior in MCP-Mod ################## 
# candidate model set in MCP-Mod
can.model.set <- c("linear", "linlog", "quadratic", "exponential", "emax", "sigmoid")
n.models <- length(can.model.set) #number of candidate model fitted
# Prior parameters in candidate model
fun.linear <- function(x){0.75*x[1] + 0.75*x[2]}
fun.linlog <- function(x){1.082 * log(x[1] + 1) + 1.082* log(x[2] + 1)}
fun.quadratic <- function(x){2.09*2*0.6*x[1] - 2.09*x[1]^2 + 2.09*2*0.6*x[2] - 2.09*x[2]^2 }
fun.exponential <- function(x){0.067*(exp(x[1]/0.4) - 1) + 0.067*(exp(x[2]/0.4) - 1)}
fun.emax <- function(x){0.9*x[1] /(0.2 + x[1]) + 0.9*x[2]/(0.2 + x[2])}
fun.sigmoid <- function(x){0.87*x[1]^2 /(0.4^2 + x[1]^2) + 0.87*x[2]^2/(0.4^2 + x[2]^2)}
# calculate optimal contrast coefficients
mu1 <- apply(doses[,1:2], 1, fun.linear)
mu2 <- apply(doses[,1:2], 1, fun.linlog)
mu3 <- apply(doses[,1:2], 1, fun.quadratic)
mu4 <- apply(doses[,1:2], 1, fun.exponential)
mu5 <- apply(doses[,1:2], 1, fun.emax)
mu6 <- apply(doses[,1:2], 1, fun.sigmoid)
mu <- cbind(mu1, mu2, mu3, mu4, mu5, mu6)
colnames(mu) <- can.model.set
rownames(mu) <- rownames(doses)
get.copt <- function(x){
  a <- (x - mean(x))/sqrt(sum((x - mean(x))^2))
  if(sum(x*a) >= 0){opt.contrast <- a}
  if(sum(x*a) < 0){opt.contrast <- -a}
  return(opt.contrast)
}  
copt <- apply(mu, 2, get.copt)
corr <- t(copt) %*% copt
muMat <- data.frame(mu)
contMat <-  copt # contrast coefficients \boldsymbol{c}_opt
corMat <-  corr # contrast correlation matrix \boldsymbol{R}
# critical value for multivariate t distribution T_1,...,T_M
CV <- qmvt(1 - alpha, df = n*nrow(doses) - nrow(doses), tail = "lower", corr = corMat)$quantile
inter.effect <- "additive"
dir.create(inter.effect)
setwd(inter.effect)
dir.create(paste0("output_n", n)) # in case there is no output directory
setwd(paste0("output_n", n))
for(loop in jobId){
  # name <- list.files(".")
  # num.ex <- grep(paste0("out_", model,"*"), name)
  #  if(loop %in% num.ex == 0){
  out.B <- NULL
  out.coef.mcp <- NULL
  out.coef.aov <- NULL
  for(i.rep in 1:B){
    ################## Run simulation ###################
    ### vector of estimated mean treatment effect mu_m,ij
    mu.hat.mcp <- rep(NA, k)
    mu.hat.aov <- rep(NA, k)
    ### ID of selected model
    order.select.model.mcp <- NA
    order.select.model.aov <- NA
    order.used.model.mcp <- NA
    order.used.model.aov <- NA
    ### generate the simulated dataset
    data <- data.frame(d1 = rep(dose.d1, each = length(dose.d1)*n),
                       d2 = rep(rep(dose.d2, each = n), length(dose.d2)),
                       level = rep(1:nrow(doses), each = n))
    ### generate response y
    data$resp <- apply(data[, 1:2], 1, gen_arm)
    ### fit linear models using lm to avoid convergence problem:
    flag.fit <- rep(0, length(can.model.set))
    fit.linear <- lm(resp ~ d1 + d2, data = data)
    flag.fit[1] <- 1
    fit.linlog <- lm(resp ~ log(d1+1) + log(d2+1), data = data)
    flag.fit[2] <- 1
    ### fit non-linear models using nls, might be convergence error:
    fit.quadratic <- NA
    if(class(try(
      nls(resp~E0+delta1*d1+delta2*d1^2+delta3*d2+delta4*d2^2,
          start= list(E0=0, delta1=0, delta2=0, delta3=0, delta4=0),
          data=data), silent = T
    ))!='try-error'){
      fit.quadratic <- nls(resp~E0+delta1*d1+delta2*d1^2+delta3*d2+delta4*d2^2,
                           start= list(E0=0, delta1=0, delta2=0, delta3=0, delta4=0),
                           data=data)
      flag.fit[3] <- 1
    }
    fit.exponential <- NA
    if(class(try(
      nls(resp~E0+alpha1*(exp(d1/delta1)-1)+alpha2*(exp(d2/delta2)-1),
          start= list(E0=0.5, alpha1=0.5, delta1=8, alpha2=0.5, delta2=8),
          data=data), silent = T
    ))!='try-error'){
      fit.exponential <-nls(resp~E0+alpha1*(exp(d1/delta1)-1)+alpha2*(exp(d2/delta2)-1),
                            start= list(E0=0.5, alpha1=0.5, delta1=8, alpha2=0.5, delta2=8),
                            data=data)
      flag.fit[4]<-1
    }
    fit.emax <- NA
    if(class(try(
      nls(resp ~ E0+ emax.1*d1/(d1+ed50.1) + emax.2* d2/(d2+ed50.2), 
          start= list(E0=0, emax.1=1.5, emax.2=1.5, ed50.1=0.7, ed50.2=0.7),
          data= data), silent = T
    ))!='try-error'){
      fit.emax <- nls(resp ~ E0+ emax.1*d1/(d1+ed50.1) + emax.2* d2/(d2+ed50.2), 
                      start= list(E0=0, emax.1=1.5, emax.2=1.5, ed50.1=0.7, ed50.2=0.7),
                      data= data)
      flag.fit[5]<-1
    }
    fit.sigmoid <- NA
    if(class(try(
      nls(resp ~ E0+ alpha1*d1^2/(d1^2+delsq1) + alpha2*d2^2/(d2^2+delsq2), 
          start= list(E0=0, alpha1=0.5, delsq1=0.2, alpha2=0.5, delsq2=0.2),
          data= data), silent = T
    ))!='try-error'){
      fit.sigmoid <- nls(resp ~ E0+ alpha1*d1^2/(d1^2+delsq1) + alpha2*d2^2/(d2^2+delsq2), 
                         start= list(E0=0, alpha1=0.5, delsq1=0.2, alpha2=0.5, delsq2=0.2),
                         data= data)
      flag.fit[6]<-1
    }
    list.fit <- list(fit.linear, fit.linlog, fit.quadratic, fit.exponential, fit.emax, fit.sigmoid)
    ### calculate AIC for models can be fitted
    aic <- rep(NA, length(can.model.set))
    for(i in 1:length(can.model.set)){
      if(flag.fit[i]==1){
        eval(parse(text=paste0("aic[i]<- AIC(fit.", can.model.set[i], ")")))}
    }
    names(aic) <- can.model.set
    ### MCP: single contrast test
    dat <- data %>% group_by(level) %>% summarise(ybar = mean(resp), n = n())
    new.data <- data %>% group_by(level) %>% mutate(ybar = mean(resp))
    ### calculate T-stat
    s <- sqrt(sum((new.data$resp - new.data$ybar)^2)/(n*nrow(doses) - nrow(doses)))
    get_t.stat <- function(x){sum(dat$ybar * x) / (s*sqrt(sum(x ^ 2 /n)))}
    T.stats <- apply(contMat, 2, get_t.stat)
    reject.null <- T.stats > CV
    n.model.sig <- sum(reject.null) # number of significant models
    coef.mcp <- NULL
    ### if no model is significant
    if(n.model.sig ==  0){
      detect.DR.mcp  <- 0 
      detect.dose.mcp <- 0
    }
    ### if as least on model is significant
    if(n.model.sig > 0){
      ### Pr(DR) = percent of detect DR 
      detect.DR.mcp <- 1
      ### select best model: the one with the largest T statistics
      select.model.mcp <- names(which(T.stats==max(T.stats)))
      order.select.model.mcp <- which(can.model.set == select.model.mcp) #the order of selected model in candidate model set
      eval(parse(text = paste0("fit.selected <- fit.", select.model.mcp)))
      ### use the model in fitting (in case the selected model is not converged)
      #within the significant models, search the largest T stat
      order.Tstats <- order(T.stats, decreasing = T)
      fit.used <- NA
      for(i.mod in order.Tstats){
        if(class(list.fit[[i.mod]])!="logical" & T.stats[i.mod]>CV){
          fit.used <- list.fit[[i.mod]]
          order.used.model.mcp <- i.mod
          used.model.mcp <- can.model.set[i.mod]
          break
        }
      }
      if(is.na(order.used.model.mcp )==1){
        detect.dose.mcp <- 0
      }
      ### target dose (this part can be improved)
      if(is.na(order.used.model.mcp )==0){
        data.all.d <- data.frame(d1 = rep(seq(0, 1, by = 0.001), each = 1001), 
                                 d2 = rep(seq(0, 1, by = 0.001), 1001),
                                 level = NA, resp = NA)
        
        eval(parse(text = paste0("data.all.d$yhat <- predict(fit.", used.model.mcp, ", data.all.d)")))
        eval(parse(text = paste0("mu.hat.mcp <- predict(fit.", used.model.mcp, ", doses)")))
        ### Pr(dose) = percent of obtaining TD
        detect.dose.mcp <- max(data.all.d$yhat) > Delta + mu.hat.mcp[1] #if > 0, then exist the target dose
        ### evaluate the parameter estimation
        if(used.model.mcp==model){
          if(model=='linear'){
            out.coef <- fit.linear$coefficients 
            #(Intercept), d1, d2. true: 0, 0.75, 0.75
          }
          if(model=='linlog'){
            out.coef <- fit.linlog$coefficients 
            #(Intercept), d1.lg, d2.lg. true: 0, 1.082, 1.082
          }
          if(model=='quadratic'){
            out.coef <- summary(fit.quadratic)$coefficients[,1]
            #(Intercept), d1, d1.sq, d2, d2.sq. true: 0, 2.508, - 2.09, 2.508, - 2.09
          }
          if(model=='exponential'){
            out.coef <- summary(fit.exponential)$coefficients[,1]
            #(Intercept), d1.exp, d2.exp, Delta_d1, Delta_d2. true: 0, 0.067, 0.067, 0.4, 0.4 
          }
          if(model=='emax'){
            out.coef <- summary(fit.emax)$coefficients[,1]
            #(Intercept), d1.ed, d2.ed, ED50_d1, ED50_d2. true:0, 0.9, 0.9, 0.2, 0.2
          }
          if(model=='sigmoid'){
            out.coef <- summary(fit.sigmoid)$coefficients[,1]
            #(Intercept), d1.sig, d2.sig, ED50_d1, ED50_d2. true: 0, 0.87, 0.87, 0.4, 0.4  
          }
          coef.mcp <- out.coef 
        }
      }
    }
    
    ### ANOVA with Dunnett's test:
    new.data <- data[data$level %in% stat.comp.level,]
    new.data$group <- factor(rep(paste0("b", 0:(length(stat.comp.level) - 1)), each = n))
    data.aov <- aov(resp ~ group, data = new.data) # intercept?
    ### Dunnett's test
    data.mc <-  glht(data.aov, linfct = mcp(group = "Dunnett"), alternative = "greater")
    ### Pr(DR) = percent of detect DR 
    detect.DR.aov <- sum(as.numeric(summary(data.mc)$test$pvalues) < alpha) > 0
    coef.aov <- NULL
    if(detect.DR.aov == 1){
      ### Pr(dose) = percent of find TD
      detect.dose.aov <- sum((as.numeric(summary(data.mc)$test$coefficients) > Delta) * 
                               (as.numeric(summary(data.mc)$test$pvalues) < alpha)) > 0
      if(length(which(aic == min(aic, na.rm = T)))==1){
        select.model.aov <- can.model.set[which(aic == min(aic, na.rm = T))] #select the model with smallest AIC
      }
      if(length(which(aic == min(aic, na.rm = T)))>1){
        select.model.aov <- can.model.set[sample(which(aic == min(aic, na.rm = T)), 1)] #select the model with smallest AIC
      }
      order.select.model.aov <- which(can.model.set == select.model.aov) #the order of selected model in candidate model set
      used.model.aov <- select.model.aov 
      order.used.model.aov <- order.select.model.aov 
      eval(parse(text = paste0("mu.hat.aov <- predict(fit.", used.model.aov, ", doses)")))
      ### evaluate the parameter estimation
      if(used.model.aov==model){
        if(model=='linear'){
          out.coef <- fit.linear$coefficients 
          #(Intercept), d1, d2. true: 0, 0.75, 0.75
        }
        if(model=='linlog'){
          out.coef <- fit.linlog$coefficients 
          #(Intercept), d1.lg, d2.lg. true: 0, 1.082, 1.082
        }
        if(model=='quadratic'){
          out.coef <- summary(fit.quadratic)$coefficients[,1]
          #(Intercept), d1, d1.sq, d2, d2.sq. true: 0, 2.508, - 2.09, 2.508, - 2.09
        }
        if(model=='exponential'){
          out.coef <- summary(fit.exponential)$coefficients[,1]
          #(Intercept), d1.exp, d2.exp, Delta_d1, Delta_d2. true: 0, 0.067, 0.067, 0.4, 0.4 
        }
        if(model=='emax'){
          out.coef <- summary(fit.emax)$coefficients[,1]
          #(Intercept), d1.ed, d2.ed, ED50_d1, ED50_d2. true:0, 0.9, 0.9, 0.2, 0.2
        }
        if(model=='sigmoid'){
          out.coef <- summary(fit.sigmoid)$coefficients[,1]
          #(Intercept), d1.sig, d2.sig, ED50_d1, ED50_d2. true: 0, 0.87, 0.87, 0.4, 0.4  
        }
        coef.aov <- out.coef 
      }
    }
    if(detect.DR.aov == 0){
      detect.dose.aov <- 0}
    # results to be saved into 'out' object
    out <- c(detect.DR.aov, detect.DR.mcp, detect.dose.aov, detect.dose.mcp, 
             order.select.model.mcp, order.select.model.aov,
             mu.hat.mcp, mu.hat.aov, T.stats,
             order.used.model.mcp, order.used.model.aov)
    names(out) <- c("detect.DR.aov", "detect.DR.mcp", "detect.dose.aov", "detect.dose.mcp", 
                    "select.model.mcp", "select.model.aov", 
                    paste0("mu", 1:k, ".mcp"),
                    paste0("mu", 1:k, ".aov"),
                    paste0("T_", can.model.set),
                    "used.model.mcp", "used.model.aov")
    out.B <- rbind(out.B, out)
    out.coef.mcp <- rbind(out.coef.mcp, coef.mcp)
    out.coef.aov <- rbind(out.coef.aov, coef.aov)
    #save(doses, can.model.set, out, file = paste0("out_", model,"_", ((jobId-1)*B+i.rep), ".RData"))
    rm(data, 
       detect.DR.aov, detect.DR.mcp, detect.dose.aov, detect.dose.mcp, 
       order.select.model.mcp, order.select.model.aov,
       mu.hat.mcp, mu.hat.aov,
       T.stats, order.used.model.mcp, order.used.model.aov,
       coef.mcp, coef.aov)
  }
  names(out.B) <- c("detect.DR.aov", "detect.DR.mcp", "detect.dose.aov", "detect.dose.mcp", 
                    "select.model.mcp", "select.model.aov", 
                    paste0("mu", 1:k, ".mcp"),
                    paste0("mu", 1:k, ".aov"),
                    paste0("T_", can.model.set))
  save(doses, can.model.set, out.B, 
       out.coef.mcp, out.coef.aov, file = paste0("out_", model,"_", loop, ".RData"))
  #}
}
