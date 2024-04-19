library(multcomp)
library(dplyr)
################## k=16 ################## 
# Possible dose for drug A and drug B
dose.d1 <- c(0, 0.2, 0.6, 1)
dose.d2 <- c(0, 0.2, 0.6, 1)
d1 <- rep(dose.d1, each = length(dose.d1))
d2 <- rep(dose.d2, length(dose.d2))
doses <- data.frame(d1 = d1, d2 = d2)
# dose levels: k
doses$level <- 1:nrow(doses)
k <- max(doses$level)
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
cor.k16 <- as.data.frame(corMat)
################## k=9 ################## 
# Possible dose for drug A and drug B
dose.d1 <- c(0, 0.6, 1)
dose.d2 <- c(0, 0.6, 1)
d1 <- rep(dose.d1, each = length(dose.d1))
d2 <- rep(dose.d2, length(dose.d2))
doses <- data.frame(d1 = d1, d2 = d2)
# dose levels: k
doses$level <- 1:nrow(doses)
k <- max(doses$level)
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
cor.k9 <- as.data.frame(corMat)
cor.k9$dose <- 9
cor.k16$dose <- 16
res <- rbind(cor.k9, cor.k16)
write.csv(res, file = "../data/result-additive/contrast_correlation.csv", row.names = F, col.names = F)
