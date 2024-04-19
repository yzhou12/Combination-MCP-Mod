dir.create("../script")
dir.create("../data")
setwd("../script")
dr.model <- c( "flat", "linear", "linlog", "quadratic", "exponential","emax","sigmoid")
j=1
for(i in 1:7){
  data1 <- scan(paste0("../source/template_simulation.R"), what = character(), sep = "\n")
  data1[16] <- paste0("model <- dr.model[", i, "]")
  for(n in c(5,10,15,20)){
    data1[17] <- paste0("n <- ", n)
    write.table(data1, paste0("all_simulation_", dr.model[i], "_n", n, "_k9.R"), row.names = F, col.names = F, quote = F)
    data2 <- scan("../source/test.sh", what = character(), sep = "\n")
    data2[6] <- paste0("Rscript all_simulation_", dr.model[i], "_n", n, "_k9.R $LSB_JOBINDEX")
    write.table(data2, paste0("s_", j, ".sh"), row.names = F, col.names = F, quote = F)
    j=j+1
  }
}
for(i in 1:7){
  data1 <- scan(paste0("../source/template_simulation.R"), what = character(), sep = "\n")
  data1[8] <- "dose.d1 <- c(0, 0.2, 0.6, 1)"
  data1[9] <- "dose.d2 <- c(0, 0.2, 0.6, 1)"
  data1[16] <- paste0("model <- dr.model[", i, "]")
  for(n in c(5,10,15,20)){
    data1[17] <- paste0("n <- ", n)
    write.table(data1, paste0("all_simulation_", dr.model[i], "_n", n, "_k16.R"), row.names = F, col.names = F, quote = F)
    data2 <- scan("../source/test.sh", what = character(), sep = "\n")
    data2[6] <- paste0("Rscript all_simulation_", dr.model[i], "_n", n, "_k16.R $LSB_JOBINDEX")
    write.table(data2, paste0("s_", j, ".sh"), row.names = F, col.names = F, quote = F)
    j=j+1
  }
}

res <- "module load ib\nmodule load R/4.1.2"
for(i in 1:56){
  res <- c(res, paste0("bsub -J s", i, "[1-10] < s_", i,".sh"))
}
write.table(res, "command.txt", row.names = F, col.names = F, quote = F)
