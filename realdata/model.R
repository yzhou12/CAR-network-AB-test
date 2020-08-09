library(igraph)
library(Matrix)
library(dplyr)
B <- 1:1000 # number of simulations
setting <- 1 # different simulation setting
###################### PARAMETER SETTING ######################
sd.err <- 2 # sd of random error
mu1 <- 1 # treatment effect of A
mu0 <- 0 # treatment effect of B
set.beta <- 1 # cluster effect 
p <- 0.85 # biased coin probability
set.feature <- c(1, 2, 5, 7) # id of response model covariates in cluster feature matrix
p.variable <- length(set.feature)
source("setting_para.R")
alpha1 <- -alpha0 #network effect from B
work.feature <- set.feature[1:p.work]
beta <- rep(set.beta, p.variable)
ncont <- p.work
###################### END ######################

###################### FUNCTIONS ######################
fun_sample_nj <- function(k1, k2){
  k <- k1:k2
  if(k1 == k2){
    nj <- k
  }
  if((k1 != k2) && length(k)%%2 == 0){ #even
    nj <- sample(k, 1, prob = (1/(abs(k - (k1+k2)/2)))/sum(1/(abs(k - (k1+k2)/2))))
  }
  if((k1 != k2) && length(k)%%2 == 1){ #odd
    nj <- sample(k, 1, prob = (1/(abs(k - (k1+k2)/2) + 0.5))/sum(1/(abs(k - (k1+k2)/2) + 0.5)))
  }
  return(nj)     
}
delta <- function(vec.x, com.group){
  # delta(ei,ej) = 1 if ci=cj
  return(com.group$membership[vec.x[1]] == com.group$membership[vec.x[2]])
}
DiffCal.Mah <- function(x.data, n.assign){
  newdata <- x.data[1:n.assign, ]
  subdata_A <- newdata[which(newdata$groupA == 1), 1:p.variable]
  subdata_B <- newdata[which(newdata$groupA == 0), 1:p.variable]
  xbar_A <- as.matrix(colMeans(subdata_A))
  xbar_B <- as.matrix(colMeans(subdata_B))
  temp <- as.matrix(newdata[, 1:p.variable])
  if(det(cov(temp)) == 0){
    Dc <- 0
  }
  if(det(cov(temp)) != 0){
    Dc <- n.assign* (nrow(subdata_B)/n.assign) * (nrow(subdata_A)/n.assign) *t(xbar_A - xbar_B) %*% solve(cov(temp), tol = 1e-99) %*% (xbar_A - xbar_B)
  }
  M <- Dc
  return(M)
}
#dichotomize continuous covariate by median
DICH.coun <- function(x){
  return(ifelse(x <= median(x), 0, 1))
}
#get dummy variables for each margin
margin_dummy <- function(dataline){
  l.mar <- as.factor(0:1)
  x.l.mar <- l.mar[dataline + 1]
  dummy.mar <- dummy(x.l.mar, drop = F)
} 
DiffCal.Margin <- function(x.data, n.assign, dummy.mar){
  x <- x.data[1:n.assign, ]
  # overall
  Doverall <- sum(x$groupA)-sum(x$groupB)
  # marginal
  x.dummy.mar <- dummy.mar[1:n.assign,]
  margin_A <- t(as.matrix(x$groupA)) %*% as.matrix(x.dummy.mar)
  margin_B <- t(as.matrix(x$groupB)) %*% as.matrix(x.dummy.mar)
  diff_margin <- as.vector(margin_A - margin_B)
  
  Dnew <- diff_margin
  #str.name <- paste0("names(Dnew) <- c(", paste0(paste0("\"Dz", 1:nbina, "0\", ", "\"Dz", 1:nbina ,"1\" "), collapse = ", "), "\")")
  #eval(parse(text = str.name))
  return(Dnew)
}
get_result <- function(method){
  mah <- DiffCal.Mah(data.c, nrow(data.c), p.work)
  design.x <- as.matrix(data.n[, 1:p.variable])  #INPUT: \bX groupA groupB cluster
  fun.nz1 <- function(vec){
    return(vec%*% as.matrix(data.n$groupA))
  }
  fun.nz0 <- function(vec){
    return(vec%*% as.matrix(data.n$groupB))
  }
  temp.nz1 <- apply(g.all.mat, 1, fun.nz1)
  temp.nz0 <- apply(g.all.mat, 1, fun.nz0)
  temp.nz0 <- temp.nz0*data.n$groupA
  temp.nz1 <- temp.nz1*data.n$groupB
  temp.y <-cbind(data.n$groupA, data.n$groupB) %*% c(mu1, mu0) + alpha0 * temp.nz1 + alpha1 * temp.nz0 +
    design.x %*% beta + rnorm(n = nrow(data.n), mean = 0, sd = sd.err)
  clus <- data.n
  clus$y <- temp.y
  clus$nz1 <- temp.nz1
  clus$nz0 <- temp.nz0
  clus.c <- clus %>% group_by(cluster) %>% summarise(weiy = mean(y))
  clus.c <- left_join(data.c, clus.c, by = "cluster")
  tau1 <- mean(clus[clus$groupA == 1,]$y) - mean(clus[clus$groupB == 1,]$y) #node-level un-adjusted
  tau2 <- mean(clus.c[clus.c$groupA == 1,]$weiy) - mean(clus.c[clus.c$groupB == 1,]$weiy) #cluster-level un-adjusted
  type3 <- clus[clus$nz0 == 0 & clus$nz1 == 0,] #no bad edge
  type3.c <- type3 %>% group_by(cluster) %>% summarise(weiy = mean(y, na.rm = T))
  type3.c <- left_join(data.c, type3.c, by = "cluster")
  tau3 <- mean(type3[type3$groupA == 1,]$y, na.rm = T) - mean(type3[type3$groupB == 1,]$y, na.rm = T) #node-level adjusted
  tau4 <- mean(type3.c[type3.c$groupA == 1,]$weiy, na.rm = T) - mean(type3.c[type3.c$groupB == 1,]$weiy, na.rm = T) #cluster-level un-adjusted
  if(method != "CRU"){
    psey <- as.matrix(clus.c[,1:p.work]) %*% rep(set.beta, p.work) + cbind(clus.c$groupA, clus.c$groupB) %*% c(mu1, mu0) + rnorm(n = nrow(data.c), mean = 0, sd = sd.err)
    clus.c$psey <- psey
    eval(parse(text = paste0("lmc.a <- lm(data = clus.c[clus.c$groupA == 1, ], psey ~ ", 
                             paste("c", 1:p.work, collapse = " + ",sep = "" ), ")")))
    eval(parse(text = paste0("lmc.b <- lm(data = clus.c[clus.c$groupB == 1, ], psey ~ ", 
                             paste("c", 1:p.work, collapse = " + ",sep = "" ), ")")))
    r2.c <- mean(summary(lmc.b)$r.squared, summary(lmc.a)$r.squared)
  }
  if(method == "CRU"){
    r2.c <- NA
  }
  clus$al1 <- clus$nz0 * clus$groupA
  clus$al0 <- clus$nz1 * clus$groupB
  eval(parse(text = paste0("est.fit <- lm(y ~ groupA + groupB + al0 + al1 + ", 
                           paste("c", 1:p.variable, collapse = " + ",sep = "" )," - 1, clus)")))
  alpha0.reg <- coef(est.fit)[3]
  alpha1.reg <- coef(est.fit)[4]
  result <- c(tau1, tau2, tau3, tau4, alpha0.reg, alpha1.reg, r2, mah, r2.c)
  return(result)
}
###################### END ######################

dir.create("data")
result.om.all <- NULL
result.crc.all <- NULL
result.cru.all <- NULL
for(l in B){
  ###################### GEN DATA ######################
  load("realitymining.RData")
  N <- nrow(table.info.true)
  data.id <- sample(1:nrow(data), nrow(data), replace = F)
  data <- data[data.id,]
  colnames(data) <- paste("c", 1:p.variable, sep = "")
  data$groupA <- 0
  data$groupB <- 0
  ###################### END ######################
  ##################### MS #####################
  data.x <- data[, c(1:p.work, (p.variable + 1), ncol(data))]
  data <- data[order(as.numeric(rownames(data))),]
  data <- CAR(data.x) 
  data.c$groupA <- data$groupA
  data.c$groupB <- data$groupB
  data.n <- left_join(data.group.true, data.c, by = "cluster")
  data.n <- data.n[,c(3:(3 + p.variable + 1),2)] 
  result.om <- get_result("OM")
  ###################### END ######################
  ##################### CRC #####################
  data$groupA <- 0
  data$groupB <- 0
  id <- sample(1:nrow(data), nrow(data)/2 , replace = F)
  data[id,]$groupA <- 1
  data$groupB <- 1- data$groupA
  data.c$groupA <- data$groupA
  data.c$groupB <- data$groupB
  data.n <- left_join(data.group.true, data.c, by = "cluster")
  data.n <- data.n[,c(3:(3 + p.variable + 1),2)] # \bXi, groupA, groupB, cluster 
  result.crc <- get_result("CRC")
  ###################### END ######################
  
  ##################### CRU #####################
  data.n$groupA <- 0 
  data.n$groupB <- 0
  Nnode <- nrow(data.n)
  if(Nnode %% 2 == 0){
    groupA.node <- sample(1:Nnode, Nnode/2, replace = F)
  }
  if(Nnode %% 2 == 1){  #odd number of clusters
    random <- sample(c(0,1), 1)
    if(random == 1){
      groupA.node <- sample(1:Nnode, (Nnode-1)/2, replace = F)}
    if(random == 0){
      groupA.node <- sample(1:Nnode, (Nnode+1)/2, replace = F)}
  }
  data.n$groupA[groupA.node] <- 1
  data.n$groupB <- 1 - data.n$groupA
  data.c$groupA <- 0
  data.c$groupB <- 0
  result.cru <- get_result("CRU")
  ###################### END ######################
  result.om.all <- rbind(result.om.all, result.om)
  result.crc.all <- rbind(result.crc.all, result.crc)
  result.cru.all <- rbind(result.cru.all, result.cru)
  save(result.om.all, result.crc.all, result.cru.all, file = paste0("data/tau_s", setting, ".RData"))
}
