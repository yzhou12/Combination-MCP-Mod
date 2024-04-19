### dose levels for drug A and B in 3D plots:
d1 <- seq(0, 1, by = 0.05)
d2 <- d1
### treatment groups dose levels: 9 dose level as an example
dose.d1 <- c(0, 0.6, 1)
dose.d2 <- c(0, 0.6, 1)
doses <- data.frame(d1 = rep(dose.d1, each = length(dose.d1)), 
                    d2 = rep(dose.d2, length(dose.d2)))

jpeg("true-additive.jpg", res = 500, width = 8, height = 6, unit = "in")
par(mfrow=c(2,4))
par(mai=c(0.2, 0, 0.2, 0))
### Additive effect:
### flat
gen.flat <- function(x, y){
  0*x*y
}
persp(x = d1, y = d2, z = outer(d1, d2, gen.flat), phi = 45, theta = 45,
      xlab = "d1", ylab = "d2", zlim = c(0, 1.5),
      zlab = "resp", ticktype = "detailed", nticks = 2, main = "Flat")
### Linear
gen.linear <- function(x, y){
  0.75*x + 0.75*y
}
persp(x = d1, y = d2, z = outer(d1, d2, gen.linear), phi = 45, theta = 45,
      xlab = "d1", ylab = "d2", zlim = c(0, 1.5),
      zlab = "resp", ticktype = "detailed", nticks = 2, main = "Linear")
### Linear log
gen.linlog <- function(x, y){
  1.082 * log(x + 1) + 1.082* log(y + 1)
}
persp(x = d1, y = d2, z = outer(d1, d2, gen.linlog), phi = 45, theta = 45,
      xlab = "d1", ylab = "d2", zlim = c(0, 1.5),
      zlab = "resp", ticktype = "detailed", nticks = 2, main = "Linear in log") 
### Quadratic
gen.quadratic <- function(x, y){
  2.09*2*0.6*x - 2.09*x^2 + 2.09*2*0.6*y - 2.09*y^2 
}
persp(x = d1, y = d2, z = outer(d1, d2, gen.quadratic), phi = 45, theta = 45,
      xlab = "d1", ylab = "d2", zlim = c(0, 1.5),
      zlab = "resp", ticktype = "detailed", nticks = 2, main = "Quadratic") 
### Exponential
gen.exponential <- function(x, y){
  0.067*(exp(x/0.4) - 1) + 0.067*(exp(y/0.4) - 1)
}
persp(x = d1, y = d2, z = outer(d1, d2, gen.exponential), phi = 45, theta = 45,
      xlab = "d1", ylab = "d2", zlim = c(0, 1.5),
      zlab = "resp", ticktype = "detailed", nticks = 2, main = "Exponential") 
### Emax
gen.emax <- function(x, y){
  0.9*x /(0.2 + x) + 0.9*y/(0.2 + y)
}
persp(x = d1, y = d2, z = outer(d1, d2, gen.emax), phi = 45, theta = 45,
      xlab = "d1", ylab = "d2", zlim = c(0, 1.5),
      zlab = "resp", ticktype = "detailed", nticks = 2, main = "Emax") 
### Sigmoid
gen.sigmoid <- function(x, y){
  0.87*x^2 /(0.4^2 + x^2) + 0.87*y^2/(0.4^2 + y^2)
}
persp(x = d1, y = d2, z = outer(d1, d2, gen.sigmoid), phi = 45, theta = 45,
      xlab = "d1", ylab = "d2", zlim = c(0, 1.6),
      zlab = "resp", ticktype = "detailed", nticks = 2, main = "Sigmoid") 
#max(outer(d1, d2, gen.sigmoid))
### 6 step functions
### dose levels for drug A and B in 3D plots:
d1 <- seq(0, 1, by = 0.05)
d2 <- d1
### treatment groups dose levels: 9 dose level as an example
dose.d1 <- c(0, 0.6, 1)
dose.d2 <- c(0, 0.6, 1)
doses <- data.frame(d1 = rep(dose.d1, each = length(dose.d1)), 
                    d2 = rep(dose.d2, length(dose.d2)))

#gen.step <- function(x, y){
#  ifelse(x + y <= a, 0, 1.5)
#}
#a <- 0.8
#persp(x = d1, y = d2, z = outer(d1, d2, gen.step), phi = 45, theta = 45,
#      xlab = "d1", ylab = "d2", zlim = c(0, 1.5),
#      zlab = "resp", ticktype = "detailed", nticks = 2, main = "step")
dev.off()

