tau <- 1
setting <- 1 # different simulation setting
source("setting_para.R")
method <- c("om", "crc", "cru")
table1 <- data.para
table2 <- data.para
table3 <- data.para
table4 <- data.para
j <- 1
for(s in slist){
  load(paste0("data/tau_s", setting, ".RData"))
  colnames(result.om.all) <- c("tauvu", "taucu", "tauva", "tauca", "alpha0.reg", "alpha1.reg", "mah", "r2.c")
  colnames(result.crc.all) <- c("tauvu", "taucu", "tauva", "tauca", "alpha0.reg", "alpha1.reg", "mah", "r2.c")
  colnames(result.cru.all) <- c("tauvu", "taucu", "tauva", "tauca", "alpha0.reg", "alpha1.reg", "mah", "r2.c")
  # CE estimator:
  for(k in 1:4){
    for(m in 1:3){
      eval(parse(text = paste0("table", k, "$bia.", method[m], "[j] <- mean(result.", method[m], ".all[,", k, "] - tau)")))
      eval(parse(text = paste0("table", k, "$SD.", method[m], "[j] <- sd(result.", method[m], ".all[,", k, "])")))
      eval(parse(text = paste0("table", k, "$var.", method[m], "[j] <- var(result.", method[m], ".all[,", k, "])")))
    }
    eval(parse(text = paste0("table", k, "$priv[j] <- (var(result.crc.all[,", k, "]) - var(result.om.all[,", k, "]))/var(result.crc.all[,", k, "]) * 100")))
  }
  result.om.all <- as.data.frame(result.om.all)
  j <- j + 1
}
# table1, table2, table3, table4 summarize performance of 4 types of estimators
