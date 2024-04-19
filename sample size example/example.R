MCPpower <- function(n.range, dose.comb, alpha, sigma, can.model.set){
  library(mvtnorm)
  # candidate model set
  eval(parse(text = paste0("models <- list(", paste("fun.", can.model.set, sep = "", collapse = ", ") ,")")))
  n.models <- length(models) #number of candidate model fitted
  # obtain total dose level (k) 
  doses <- dose.comb
  doses$level <- 1:nrow(doses)
  k <- max(doses$level)
  # calculate optimal contrast coefficients
  for(i in 1:n.models){
    eval(parse(text = paste0("mu", i, " <- apply(doses[,1:2], 1, models[[", i, "]])")))
  }
  eval(parse(text = paste0("mu <- cbind(", paste("mu", 1:n.models, sep = "", collapse = ", "), ")")))
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
  # (PoC) power calculation
  get_power <- function(model, n){
    CV <- qmvt(1 - alpha, df = n*k - k, tail = "lower", corr = corMat)$quantile
    #non-central parameter
    delta.m <- NULL
    for(l in 1:n.models){
      eval(parse(text = paste0("temp <- sqrt(n)* sum(contMat[,l] * muMat$", model, ")/sigma")))
      delta.m <- c(delta.m, temp)
    }
    #power
    as.numeric(1 - pmvt(lower = rep(-Inf, n.models), upper = rep(CV, n.models), df = n*k - k, 
                        corr = corMat, type = "shifted", delta = delta.m))
  }
  data <- NULL
  for(size in n.range){
    power <- NULL
    for(m in can.model.set){
      power <- c(power, get_power(m, n = size))
    }
    data <- rbind(data, c(size, min(power), mean(power), max(power)))
  }
  colnames(data) <- c("n", "min", "mean", "max")
  data <- data.frame(data)
  return(data)
}

library(mvtnorm)
alpha <- 0.05
sigma <- 1.4
# Possible dose for drug A and drug B
dose.d1 <- c(0, 0.6, 1)
dose.d2 <- c(0, 0.6, 1)
d1 <- rep(dose.d1, each = length(dose.d1))
d2 <- rep(dose.d2, length(dose.d2))
doses <- data.frame(d1 = d1, d2 = d2)
# dose levels: k
doses$level <- 1:nrow(doses)
k <- max(doses$level)
################ 2D MCP-Mod ################
# Input guesstimate of model parameters
fun.linear <- function(x){0.75*x[1] + 0.75*x[2]}
fun.linlog <- function(x){1.082 * log(x[1] + 1) + 1.082* log(x[2] + 1)}
fun.quadratic <- function(x){2.09*2*0.6*x[1] - 2.09*x[1]^2 + 2.09*2*0.6*x[2] - 2.09*x[2]^2 }
fun.exponential <- function(x){0.067*(exp(x[1]/0.4) - 1) + 0.067*(exp(x[2]/0.4) - 1)}
fun.emax <- function(x){0.9*x[1] /(0.2 + x[1]) + 0.9*x[2]/(0.2 + x[2])}
# calculate power for a given sample size range 'n.range':
data <- MCPpower(n.range = 2:30, dose.comb = doses, alpha = alpha, sigma = sigma, 
         can.model.set = c("linear", "linlog","quadratic", "exponential", "emax"))
data # minimum sample size using min, mean, max: n = 9, 9, 7
jpeg("ss_mcp.jpg", res = 200, width = 5, height = 4, unit = "in")
plot(data$n, data$min*100, type = "none", xlab = "Sample size (n/arm)", ylab = "Power (%)")
lines(data$n, data$min*100, col = 2, lwd = 1.5, lty = 2)
lines(data$n, data$mean*100, col = 4, lwd = 1.5)
lines(data$n, data$max*100, col = 3, lwd = 1.5, lty = 4)
abline(h = 80, lty = 3)
legend("bottomright", legend = c("min", "mean", "max"), col = c(2,4,3), lwd = 1.5,
       lty = c(2,1,4))
dev.off()

################ Bonferroni Correction ################
alpha = 0.05/4 # there are 4 pairwise comparison: 4 active arms v.s. placebo
power = 0.8
sigma = 1.4
Delta = c(0.9, 1.2, 1.5) # To detect the maximum effect
n <- ceiling(2 * ((qnorm(1-alpha) + qnorm(power))*sigma/Delta)^2)
n #minimum sample size to detect TV 0.9, 1.2, 1.5: n= 46, 26, 17

n <- 2:50
Delta <- 0.9
power <- 1 - pnorm(qnorm(1-alpha) - Delta/(sigma*sqrt(1/n*2)))
dat <- data.frame(n = n)
dat$power1 <- power
Delta <- 1.2
power <- 1 - pnorm(qnorm(1-alpha) - Delta/(sigma*sqrt(1/n*2)))
dat$power2 <- power
Delta <- 1.5
power <- 1 - pnorm(qnorm(1-alpha) - Delta/(sigma*sqrt(1/n*2)))
dat$power3 <- power

jpeg("ss_MCP-MC.jpg", res = 400, width = 5, height = 4, unit = "in")
plot(data$n, data$min, type = "none", xlab = "Sample size (n/arm)", ylab = "Power (%)",
     ylim = c(0,100))
#lines(data$n, data$min, col = 2, lwd = 1.5, lty = 2)
lines(data$n, data$mean * 100, col = 6, lwd = 1.5)
#lines(data$n, data$max, col = 3, lwd = 1.5, lty = 4)
abline(h = 80, lty = 3)
#plot(dat$n, dat$power1, type = "n", xlab = "n/arm", ylab = "power", ylim = c(0,1))
lines(dat$n, dat$power1 * 100, col = 2, lwd = 1.5, lty = 2)
lines(dat$n, dat$power2 * 100, col = 3, lwd = 1.5, lty = 3)
lines(dat$n, dat$power3 * 100, col = 4, lwd = 1.5, lty = 4)
legend("bottomright", 
       legend = c("MCP-Mod (mean)",expression(Delta=="0.9"), expression(Delta=="1.2"), expression(Delta=="1.5")),
       col = c(6,2,3,4), lwd = 1.5,
       lty = c(1, 2,3,4), cex = 0.7)
dev.off()

