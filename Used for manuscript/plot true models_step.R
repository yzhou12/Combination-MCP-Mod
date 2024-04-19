getwd()
jpeg("true-step.jpg", res = 500, width = 6, height = 3, unit = "in")
par(mfrow=c(1,3))
par(mai=c(0.2, 0, 0.2, 0))
#max(outer(d1, d2, gen.sigmoid))
### step function
### dose levels for drug A and B in 3D plots:
d1 <- seq(0, 1, by = 0.02)
d2 <- d1
### treatment groups dose levels: 9 dose level as an example
dose.d1 <- c(0, 0.6, 1)
dose.d2 <- c(0, 0.6, 1)
doses <- data.frame(d1 = rep(dose.d1, each = length(dose.d1)), 
                    d2 = rep(dose.d2, length(dose.d2)))
gen.step <- function(x, y){
  ifelse(x + y <= a, 0, 1.5)
}
a <- 0.4
persp(x = d1, y = d2, z = outer(d1, d2, gen.step), phi = 45, theta = 45,
      xlab = "d1", ylab = "d2", zlim = c(0, 1.5),
      zlab = "resp", ticktype = "detailed", nticks = 2, main = "step 1")

gen.step <- function(x, y){
  ifelse(x + y <= a, 0, 1.5)
}
a <- 0.8
persp(x = d1, y = d2, z = outer(d1, d2, gen.step), phi = 45, theta = 45,
      xlab = "d1", ylab = "d2", zlim = c(0, 1.5),
      zlab = "resp", ticktype = "detailed", nticks = 2, main = "step 2")

gen.step <- function(x, y){
  ifelse(x + y <= a, 0, 1.5)
}
a <- 1.2
persp(x = d1, y = d2, z = outer(d1, d2, gen.step), phi = 45, theta = 45,
      xlab = "d1", ylab = "d2", zlim = c(0, 1.5),
      zlab = "resp", ticktype = "detailed", nticks = 2, main = "step 3")

dev.off()

