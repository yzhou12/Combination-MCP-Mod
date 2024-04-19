jpeg(paste0("../data/result-", inter.effect,"/bias_parameter_", model, ".jpeg"), 
     width = 8, height = 6, units = "in", res = 500)
par(mfrow = c(2, 4))
for(k in c(9,16)){
  for(n in c(5, 10, 15, 20)){
    load(paste0("../data/k", k, "/", inter.effect, "/output_n", n, "/out_", model, ".RData"))
    if(model == "quadratic"){
      coef.mcp <- as.data.frame(coef.mcp)
      coef.mcp <- coef.mcp %>% select(E0, delta1, delta3, delta2, delta4)
    }
    if(model == "exponential"){
      coef.mcp <- as.data.frame(coef.mcp)
      coef.mcp <- coef.mcp %>% select(E0, alpha1, alpha2, delta1, delta2)
    }
    if(model == "sigmoid"){
      coef.mcp <- as.data.frame(coef.mcp)
      coef.mcp <- coef.mcp %>% select(E0, alpha1, alpha2, delsq1, delsq2)
    }
    plot(colMeans(coef.mcp) - coef.true, type = "none",  ylim = ylimit,
         main = paste0("n=", n, ", k=", k), #ylab = "Prediction error",
         ylab = expression(hat(theta)-theta), xlab = "Parameter",
         xaxt="n") # mcp
    if(model %in% c("linear", "linlog")){
      axis(1, at = 1:3)}
    if(model %in% c("quadratic", "exponential","emax", "sigmoid")){
      axis(1, at = 1:5)}
    lines(apply(coef.mcp, 2, quantile, prob=0.5) - coef.true)
    lines(apply(coef.mcp, 2, quantile, prob=0.25) - coef.true, lty = 2)
    lines(apply(coef.mcp, 2, quantile, prob=0.75) - coef.true, lty = 2)
    lines(apply(coef.mcp, 2, quantile, prob=0.05) - coef.true, lty = 3)
    lines(apply(coef.mcp, 2, quantile, prob=0.95) - coef.true, lty = 3)
    abline(h = 0, col = 2, lty = 4)
  }
}
legend("topright", legend = c("50%", "25%, 75%", "5%, 95%"), lty = 1:3)
dev.off()