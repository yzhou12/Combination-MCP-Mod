library(multcomp)
library(dplyr)

################## Assumptions ################## 
alpha <- 0.05 # FWER
Delta <- 1.2 # clinical relevant effect Delta = 1.2
sigma <- 1.4 # standard deviation of response in DR model
# dose level for drug A and drug B
dose.d1 <- c(0, 0.6, 1)
dose.d2 <- c(0, 0.6, 1)
d1 <- rep(dose.d1, each = length(dose.d1))
d2 <- rep(dose.d2, length(dose.d2))
doses <- data.frame(d1 = d1, d2 = d2)
doses$level <- 1:nrow(doses)
k <- max(doses$level)
stat.comp.level <- doses$level[which(doses$d1 + doses$d2 == 0 | doses$d1 * doses$d2 != 0)] # levels used to conduct statistical comparison in ANOVA
# specify the sample size per arm
n <- 5
################## END ################## 

################## Prior in MCP-Mod ################## 
# specify the candidate model set in MCP-Mod
can.model.set <- c("linear", "linlog", "quadratic", "exponential", "emax", "sigmoid")
n.models <- length(can.model.set) #number of candidate model fitted
# Prior parameters/guessetmiation in candidate model
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
# critical value for multivariate t distribution T_1,...,T_M:
CV <- qmvt(1 - alpha, df = n*nrow(doses) - nrow(doses), tail = "lower", corr = corMat)$quantile
################## END ################## 

################## Generate synthetic trial data ################## 
# choose the underlying dose-response model (true model)
dr.model <- c( "linear", "linlog", "quadratic", "exponential","emax","sigmoid")
model <- dr.model[5]
# specify the response generating models (true models) with additive effect:
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
# generate the simulated dataset:
data <- data.frame(d1 = rep(dose.d1, each = length(dose.d1)*n),
                   d2 = rep(rep(dose.d2, each = n), length(dose.d2)),
                   level = rep(1:nrow(doses), each = n))
# generate response y:
data$resp <- apply(data[, 1:2], 1, gen_arm)
################## END ################## 



################## MCP step in Comb MCP-Mod (POC) ################## 
### pre-enter: estimated mean treatment effect mu_m,ij
mu.hat.mcp <- rep(NA, k)
### pre-enter: order of selected and used model
select.model.mcp <- NA
used.model.mcp <- NA
    
    ### fit linear models using lm :
    # generate an index set indicating whether a model can be fitted
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
      detect.DR.mcp  <- FALSE 
      detect.dose.mcp <- FALSE 
    }
    ### if as least on model is significant
    if(n.model.sig > 0){
      detect.DR.mcp <- TRUE
      ### select best model: the one with the largest T statistics
      select.model.mcp <- names(which(T.stats==max(T.stats)))
      eval(parse(text = paste0("fit.selected <- fit.", select.model.mcp)))
    }
      ################## END ################## 
     
      ################## Mod step in Comb MCP-Mod (find target dose) ################## 
      
      ### note!!! the selected model (fit.selected) can be NA due to convergence issue
      ### the model used for target dose (fit.used) can be different from fit.selected
      ### within the significant and models which can be fitted, search the model with the largest T stat
    if(n.model.sig > 0){  
    order.Tstats <- order(T.stats, decreasing = T)
      fit.used <- NA
      for(i.mod in order.Tstats){
        if(class(list.fit[[i.mod]])!="logical" & T.stats[i.mod]>CV){
          fit.used <- list.fit[[i.mod]]
          used.model.mcp <- can.model.set[i.mod]
          break
        }
      }
      if(is.na(used.model.mcp )==1){
        detect.dose.mcp <- FALSE
      }
      ### find the target dose set
      if(is.na(used.model.mcp )==0){
        data.all.d <- data.frame(d1 = rep(seq(0, 1, by = 0.001), each = 1001), 
                                 d2 = rep(seq(0, 1, by = 0.001), 1001),
                                 level = NA, resp = NA)
        
        eval(parse(text = paste0("data.all.d$yhat <- predict(fit.", used.model.mcp, ", data.all.d)")))
        eval(parse(text = paste0("mu.hat.mcp <- predict(fit.", used.model.mcp, ", doses)")))
        detect.dose.mcp <- max(data.all.d$yhat) > Delta + mu.hat.mcp[1] #if > 0, then exist the target dose
        }
      }
    ################## END ################## 
    
    ################## Multiple Comparison ################## 
    ### pre-enter: estimated mean treatment effect mu_m,ij
    mu.hat.aov <- rep(NA, k)
    ### pre-enter: order of selected and used model
    used.model.aov <- NA
    select.model.aov <- NA
    ### ANOVA with Dunnett's test:
    new.data <- data[data$level %in% stat.comp.level,]
    new.data$group <- factor(rep(paste0("b", 0:(length(stat.comp.level) - 1)), each = n))
    data.aov <- aov(resp ~ group, data = new.data) # intercept?
    ### Dunnett's test
    data.mc <-  glht(data.aov, linfct = mcp(group = "Dunnett"), alternative = "greater")
    detect.DR.aov <- sum(as.numeric(summary(data.mc)$test$pvalues) < alpha) > 0
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

      used.model.aov <- select.model.aov 
      eval(parse(text = paste0("mu.hat.aov <- predict(fit.", used.model.aov, ", doses)")))
      ### evaluate the parameter estimation
    }
    if(detect.DR.aov == 0){
      detect.dose.aov <- FALSE}
    ################## END ################## 
    
    ################## Display results ################## 
    CombMCPMod.result <- c(detect.DR.mcp, detect.dose.mcp,
             select.model.mcp, used.model.mcp, 
             detect.DR.aov, detect.dose.aov, 
             select.model.aov,used.model.aov)
    names(CombMCPMod.result) <- c("detect.POC.Comb", "detect.TargetDose.Comb", 
                    "select.model.Comb", "used.model.Comb", 
                    "detect.POC.MC", "detect.TargetDose.MC",
                    "select.model.MC", "used.model.MC")
    # Comb MCP-Mod and Multiple Comparison results
    CombMCPMod.result
    estimate.result <- c(mu.hat.mcp, mu.hat.aov)
    names(estimate.result) <- c(paste0("mu", 1:k, ".Comb"),
                                       paste0("mu", 1:k, ".Comb"))
    # estimated response at each dose level using the fitted 'used model':
    estimate.result
    