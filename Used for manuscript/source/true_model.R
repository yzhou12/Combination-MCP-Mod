### response generating models (true models) without error term
if(model == "step"){
  gen_arm <- function(x){
    ifelse(x[1] + x[2] <= 0.8, 0, 1.5 )}
}

if(model == "flat"){
  gen_arm <- function(x){0}
}
### Additive effect:
if(inter.effect == "additive"){
  if(model == "linear"){
    gen_arm <- function(x){
      0.75 * x[1] + 0.75 * x[2]}
  }
  if(model == "linlog"){
    gen_arm <- function(x){
      1.082 * log(x[1] + 1) + 1.082 * log(x[2] + 1)}
  }
  if(model == "quadratic"){
    gen_arm <- function(x){
      2.09*2*0.6*x[1] - 2.09*x[1]^2 + 2.09*2*0.6*x[2] - 2.09*x[2]^2}
  }
  if(model == "exponential"){
    gen_arm <- function(x){
      0.067*(exp(x[1]/0.4) - 1) + 0.067*(exp(x[2]/0.4) - 1)}
  }
  if(model == "emax"){
    gen_arm <- function(x){
      0.9*x[1] /(0.2 + x[1]) + 0.9*x[2]/(0.2 + x[2])}
  }
  if(model == "sigmoid"){
    gen_arm <- function(x){
      0.87*x[1]^2 /(0.4^2 + x[1]^2) + 
        0.87*x[2]^2/(0.4^2 + x[2]^2)}
  }
}

### Synergistic effect:
if(inter.effect == "synergistic"){
  if(model == "linear"){
    gen_arm <- function(x){
      0.5 * x[1] + 0.7 * x[2] + 0.3 * x[1] * x[2]}
  }
  if(model == "linlog"){
    gen_arm <- function(x){
      0.82 * log(x[1] + 1) + 1.11 * log(x[2] + 1) + 0.33 * log(x[1] + 1) * log(x[2] + 1)}
  }
  if(model == "quadratic"){
    gen_arm <- function(x){
      1.77*2*0.6*x[1] - 1.77*x[1]^2 + 2*2*0.6*x[2]  - 2*x[2]^2 + 0.4*x[1]*x[2]}
  }
  if(model == "exponential"){
    gen_arm <- function(x){
      0.03*(exp(x[1]/0.4) - 1) + 0.037*(exp(x[2]/0.4) - 1) + 0.006 * (exp(x[1]/0.4) - 1) * (exp(x[2]/0.4) - 1)}
  }
  if(model == "emax"){
    gen_arm <- function(x){
      0.75*x[1] /(0.2 + x[1]) + 0.87*x[2]/(0.2 + x[2]) + 0.21 *(x[1]/(0.2 + x[1]))*(x[2]/(0.2 + x[2]))}
  }
}
### Antagonistic effect
if(inter.effect == "antagonistic"){
  if(model == "linear"){
    gen_arm <- function(x){
      0.75 * x[1] + 0.95 * x[2] - 0.2 * x[1] * x[2]}
  }
  if(model == "linlog"){
    gen_arm <- function(x){
      0.98 * log(x[1] + 1) + 1.39 * log(x[2] + 1) - 0.3 * log(x[1] + 1) * log(x[2] + 1)}
  }
  if(model == "quadratic"){
    gen_arm <- function(x){
      2.12*2*0.6*x[1] - 2.12*x[1]^2 + 2.45*2*0.6*x[2] - 2.45*x[2]^2 - 0.4*x[1]*x[2]}
  }
  if(model == "exponential"){
    gen_arm <- function(x){
      0.099*(exp(x[1]/0.4) - 1) + 0.103*(exp(x[2]/0.4) - 1) - 0.006 * (exp(x[1]/0.4) - 1) * (exp(x[2]/0.4) - 1)}
  }
  if(model == "emax"){
    gen_arm <- function(x){
      0.87*x[1] /(0.2 + x[1]) + 1.11*x[2]/(0.2 + x[2]) - 0.21 *(x[1]/(0.2 + x[1]))*(x[2]/(0.2 + x[2]))}
  }
}