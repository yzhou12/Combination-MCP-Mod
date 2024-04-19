fun.linear <- function(x){0.75*x[1] + 0.75*x[2]}
fun.linlog <- function(x){1.082 * log(x[1] + 1) + 1.082* log(x[2] + 1)}
fun.quadratic <- function(x){2.09*2*0.6*x[1] - 2.09*x[1]^2 + 2.09*2*0.6*x[2] - 2.09*x[2]^2 }
fun.exponential <- function(x){0.067*(exp(x[1]/0.4) - 1) + 0.067*(exp(x[2]/0.4) - 1)}
fun.emax <- function(x){0.9*x[1] /(0.2 + x[1]) + 0.9*x[2]/(0.2 + x[2])}
fun.sigmoid <- function(x){0.87*x[1]^2 /(0.4^2 + x[1]^2) + 0.87*x[2]^2/(0.4^2 + x[2]^2)}

fit.linear <- lm(resp ~ d1 + d2 + d1*d2, data = data1)
fit.linlog <- lm(resp ~ d1.lg + d2.lg +  d1.lg * d2.lg, data = data1)
fit.quadratic <- lm(resp ~ d1 + d1.sq + d2 + d2.sq + d1 * d2, data = data1)
fit.exponential <- lm(resp ~ d1.exp + d2.exp + d1.exp * d2.exp, data = data1)
fit.emax <- lm(resp ~ d1.ed + d2.ed + d1.ed * d2.ed, data= data1)
fit.sigmoid <- lm(resp ~ d1.sig + d2.sig + d1.sig * d2.sig, data= data1)


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
model <- dr.model[5]
n <- 20
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
#dir.create(paste0("../data/k", k)) # in case there is no output directory
#setwd(paste0("../data/k", k))
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
################## Run simulation ###################
### vector of estimated mean treatment effect mu_m,ij
mu.hat.mcp <- rep(NA, k)
mu.hat.aov <- rep(NA, k)
### ID of selected model
order.select.model.mcp <- NA
order.select.model.aov <- NA
### generate the simulated dataset
data <- data.frame(d1 = rep(dose.d1, each = length(dose.d1)*n),
                   d2 = rep(rep(dose.d2, each = n), length(dose.d2)),
                   level = rep(1:nrow(doses), each = n))
### generate response y
data$resp <- apply(data[, 1:2], 1, gen_arm)
### nonlinear parameter estimated by optimization: exponential, emax model
### exponential model, scale parameter delta:
f1 <- function(x){
  data.f1 <- data
  data.f1$d1.exp <- (exp(data.f1$d1/x[1]) - 1)
  data.f1$d2.exp <- (exp(data.f1$d2/x[2]) - 1)
  try.fit.exp <- lm(resp ~ d1.exp + d2.exp + d1.exp*d2.exp, data = data.f1)
  a <- summary(try.fit.exp)
  return(a$sigma)
}
dmin <- nlminb(c(1, 1), f1, lower = 1*0.05, upper = 1*8)
delta <- dmin$par
### emax: ED50
f2 <- function(x){
  data.f2 <- data
  data.f2$d1.ed <- data.f2$d1 / (x[1] + data.f2$d1)
  data.f2$d2.ed <- data.f2$d2 / (x[2] + data.f2$d2)
  try.fit.emax <-  lm(resp ~ d1.ed + d2.ed + d1.ed * d2.ed, data= data.f2)
  a <- summary(try.fit.emax)
  return(a$sigma)
}
edmin <- nlminb(c(1, 1), f2, lower = 1*0.0005, upper = 1*4)
ed50 <- edmin$par
### sigmoid: ED50
f3 <- function(x){
  data.f3 <- data
  data.f3$d1.ed <- data.f3$d1^2 / (x[1]^2 + data.f3$d1^2)
  data.f3$d2.ed <- data.f3$d2^2 / (x[2]^2 + data.f3$d2^2)
  try.fit.emax <-  lm(resp ~ d1.ed + d2.ed + d1.ed * d2.ed, data= data.f3)
  a <- summary(try.fit.emax)
  return(a$sigma)
}
edmin.sigmoid <- nlminb(c(1, 1), f3, lower = 1*0.0005, upper = 1*4)
ed50.sigmoid <- edmin.sigmoid$par
#write.table(x = paste0("loop ", i.rep, " , model-" , model, " : delta = ", paste(delta, collapse = ", "), 
#                       ", ED50 = ", paste(ed50, collapse = ", "), ", ED50.sig = ", paste(ed50.sigmoid, collapse = ", ")),
#            file = paste0("error_", model, "_", i.rep, ".txt"), row.names = F, col.names = F, quote = F)
### data transformation
trans.data <- function(x){
  new.x <- cbind(x, log(x[,1] + 1), log(x[,2] + 1), 
                 x[,1]^2, x[,2]^2,
                 (exp(x[,1]/delta[1]) - 1), (exp(x[,2]/delta[2]) - 1),
                 (x[,1] /  (ed50[1] + x[,1])), (x[,2] /  (ed50[2] + x[,2])),
                 (x[,1]^2 /  (ed50.sigmoid[1]^2 + x[,1]^2)), (x[,2]^2 /  (ed50.sigmoid[2]^2 + x[,2]^2))
  )
  colnames(new.x) <- c("d1", "d2", "level", "resp","d1.lg", "d2.lg", 
                       "d1.sq", "d2.sq", 
                       "d1.exp", "d2.exp", 
                       "d1.ed", "d2.ed", 
                       "d1.sig", "d2.sig")
  return(new.x)
}
data1 <- trans.data(data)
#write.table(x = rbind(rep(i.rep, ncol(data1)), data1),
#            file = paste0("errordata_", model, "_", i.rep, ".txt"), row.names = F, col.names = F, quote = F)
### fit all models using lm to avoid convergence problem:
fit.linear <- lm(resp ~ d1 + d2 + d1*d2, data = data1)
fit.linlog <- lm(resp ~ d1.lg + d2.lg +  d1.lg * d2.lg, data = data1)
fit.quadratic <- lm(resp ~ d1 + d1.sq + d2 + d2.sq + d1 * d2, data = data1)
fit.exponential <- lm(resp ~ d1.exp + d2.exp + d1.exp * d2.exp, data = data1)
fit.emax <- lm(resp ~ d1.ed + d2.ed + d1.ed * d2.ed, data= data1)
fit.sigmoid <- lm(resp ~ d1.sig + d2.sig + d1.sig * d2.sig, data= data1)
### calculate AIC for all models
eval(parse(text=paste0("aic <- AIC(", paste("fit.", can.model.set, sep = "", collapse = ", "), ")")))
### correct df for exponential, emax and sigmoid model
aic$df[4] <- aic$df[4] + 2
aic$AIC[4] <- aic$AIC[4] + 2*2
aic$df[5] <- aic$df[5] + 2
aic$AIC[5] <- aic$AIC[5] + 2*2
aic$AIC[6] <- aic$AIC[6] + 2*2

### MCP: single contrast test
dat <- data %>% group_by(level) %>% summarise(ybar = mean(resp), n = n())
new.data <- data %>% group_by(level) %>% mutate(ybar = mean(resp))
### calculate T-stat
s <- sqrt(sum((new.data$resp - new.data$ybar)^2)/(n*nrow(doses) - nrow(doses)))
get_t.stat <- function(x){sum(dat$ybar * x) / (s*sqrt(sum(x ^ 2 /n)))}
T.stats <- apply(contMat, 2, get_t.stat)
reject.null <- T.stats > CV
n.model.sig <- sum(reject.null) # number of significant models
### if no model is significant
if(n.model.sig ==  0){
  detect.DR.mcp  <- 0 
  detect.dose.mcp <- 0
}
### if as least on model is significant
if(n.model.sig > 0){
  ### Pr(DR) = percent of detect DR 
  detect.DR.mcp <- 1
  if(n.model.sig > 1){
    ### select best model if more than 2 models are significant using AIC
    aic.mcp <- aic
    if(n.model.sig < 6){aic.mcp[which(reject.null == 0),]$AIC <- Inf}
    select.model.mcp <- can.model.set[which(aic.mcp$AIC == min(aic.mcp$AIC))] #select the model with smallest AIC
  }
  if(n.model.sig == 1){
    select.model.mcp <- can.model.set[which(reject.null == 1)]
  }
  order.select.model.mcp <- which(can.model.set == select.model.mcp) #the order of selected model in candidate model set
  ### target dose (this part can be improved)
  data.all.d <- data.frame(d1 = rep(seq(0, 1, by = 0.001), each = 1001), 
                           d2 = rep(seq(0, 1, by = 0.001), 1001),
                           level = NA, resp = NA)
  data.all.d1 <- trans.data(data.all.d)
  data.doses1 <- trans.data(cbind(doses, resp = NA))
  eval(parse(text = paste0("data.all.d1$yhat <- predict(fit.", select.model.mcp, ", data.all.d1)")))
  eval(parse(text = paste0("mu.hat.mcp <- predict(fit.", select.model.mcp, ", data.doses1)")))
  ### Pr(dose) = percent of obtaining TD
  detect.dose.mcp <- max(data.all.d1$yhat) > Delta + mu.hat.mcp[1] #if > 0, then exist the target dose
  ### evaluate the parameter estimation
  if(model=='linear'){
    out.coef <- fit.linear$coefficients 
    #(Intercept), d1, d2, d1:d2. true: 0, 0.75, 0.75, 0
  }
  if(model=='linlog'){
    out.coef <- fit.linlog$coefficients 
    #(Intercept), d1.lg, d2.lg, d1.lg:d2.lg. true: 0, 1.082, 1.082, 0
  }
  if(model=='quadratic'){
    out.coef <- fit.quadratic$coefficients
    #(Intercept), d1, d1.sq, d2, d2.sq, d1:d2. true: 0, 2.508, - 2.09, 2.508, - 2.09, 0
  }
  if(model=='exponential'){
    out.coef <- c(fit.exponential$coefficients, delta)
    #(Intercept), d1.exp, d2.exp, d1.exp:d2.exp, Delta_d1, Delta_d2. true: 0, 0.067, 0.067, 0, 0.4, 0.4 
  }
  if(model=='emax'){
    out.coef <- c(fit.emax$coefficients, ed50)
    #(Intercept), d1.ed, d2.ed, d1.ed:d2.ed, ED50_d1, ED50_d2. true:0, 0.9, 0.9, 0, 0.2, 0.2
  }
  if(model=='sigmoid'){
    out.coef <- c(fit.sigmoid$coefficients, ed50.sigmoid)
    #(Intercept), d1.sig, d2.sig, d1.sig:d2.sig, ED50_d1, ED50_d2. true: 0, 0.87, 0.87, 0, 0.4, 0.4  
  }
  out.coef.mcp <- out.coef
}

### ANOVA with Dunnett's test:
new.data <- data[data$level %in% stat.comp.level,]
new.data$group <- factor(rep(paste0("b", 0:(length(stat.comp.level) - 1)), each = n))
data.aov <- aov(resp ~ group, data = new.data) # intercept?
### Dunnett's test
data.mc <-  glht(data.aov, linfct = mcp(group = "Dunnett"), alternative = "greater")
### Pr(DR) = percent of detect DR 
detect.DR.aov <- sum(as.numeric(summary(data.mc)$test$pvalues) < alpha) > 0
if(detect.DR.aov == 1){
  ### Pr(dose) = percent of find TD
  detect.dose.aov <- sum((as.numeric(summary(data.mc)$test$coefficients) > Delta) * 
                           (as.numeric(summary(data.mc)$test$pvalues) < alpha)) > 0
  select.model.aov <- can.model.set[which(aic$AIC == min(aic$AIC))] #select the model with smallest AIC
  order.select.model.aov <- which(can.model.set == select.model.aov) #the order of selected model in candidate model set
  data.doses2 <- trans.data(cbind(doses, resp = NA))
  eval(parse(text = paste0("mu.hat.aov <- predict(fit.", select.model.aov, ", data.doses2)")))
}
if(detect.DR.aov == 0){
  detect.dose.aov <- 0}
# results to be saved into 'out' object
out <- c(detect.DR.aov, detect.DR.mcp, detect.dose.aov, detect.dose.mcp, 
         order.select.model.mcp, order.select.model.aov,
         mu.hat.mcp, mu.hat.aov, T.stats)
names(out) <- c("detect.DR.aov", "detect.DR.mcp", "detect.dose.aov", "detect.dose.mcp", 
                "select.model.mcp", "select.model.aov", 
                paste0("mu", 1:k, ".mcp"),
                paste0("mu", 1:k, ".aov"),
                paste0("T_", can.model.set))
out.B <- rbind(out.B, out)