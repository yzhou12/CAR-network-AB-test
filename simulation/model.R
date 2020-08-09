library(igraph)
library(Matrix)
library(dplyr)
B <- 1:1000 # number of simulations
setting <- 1 # different simulation setting
###################### PARAMETER SETTING ######################
N <- 500 # number of clusters
k1 <- 10 # lower bound of cluster size
k2 <- 30 # upper bound of cluster size
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
fun.sample2 <- function(x){
  return(sample(x, 2, replace = F))
}
delta <- function(vec.x, com.group){
  # delta(ei,ej) = 1 if ci=cj
  return(com.group$membership[vec.x[1]] == com.group$membership[vec.x[2]])
}
DiffCal.Mah <- function(x.data, n.assign, n.working.cov){
  newdata <- x.data[1:n.assign, ]
  subdata_A <- newdata[which(newdata$groupA == 1), 1:n.working.cov]
  subdata_B <- newdata[which(newdata$groupA == 0), 1:n.working.cov]
  xbar_A <- as.matrix(colMeans(subdata_A))
  xbar_B <- as.matrix(colMeans(subdata_B))
  temp <- as.matrix(newdata[, 1:n.working.cov])
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
CAR <- function(input.data){
  m <- ceiling(p.work/2) + 1 
  sample.m <- sample(c(0,1), m, replace = T) # 1 means T1=1, T2=0
  id <- c(rbind(sample.m, 1-sample.m))
  input.data[1:(2*m),]$groupA <- id
  input.data[1:(2*m),]$groupB <- 1-id
  for (j in m:(N%/%2-1)){
    i <- 2*j+1
    xA <- input.data
    xA[i,]$groupA <- 1
    xA[i + 1, ]$groupB <- 1
    xB <- input.data
    xB[i,]$groupB <- 1
    xB[i + 1,]$groupA <- 1
    
    D_A <- DiffCal.Mah(xA, (i+1), p.work)
    D_B <- DiffCal.Mah(xB, (i+1), p.work)
    
    if(D_A < D_B){Tnew = rbinom(1, 1, p)}
    if(D_A > D_B){Tnew = rbinom(1, 1, 1-p)}
    if(D_A == D_B){Tnew = rbinom(1, 1, 0.5)}
    input.data[i,]$groupA <- Tnew
    input.data[i,]$groupB <- 1-Tnew
    input.data[i + 1,]$groupA <- 1-Tnew
    input.data[i + 1,]$groupB <- Tnew
  }
  if(N%%2 != 0){
    input.data[N,]$groupA <- rbinom(1,1,0.5)
    input.data[N,]$groupB <- 1 - data[N,]$groupA
  }
  return(input.data)
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
  member.true <- NULL
  for(i.com in 1:N){
    size.com <- fun_sample_nj(k1, k2)
    g <- watts.strogatz.game(dim=1, size = size.com, nei= sample(c(2,3,4,5),1), p= runif(1))
    if(i.com == 1){g.all <- g}
    if(i.com != 1){g.all <- g.all + g}
    member.true <- c(member.true, rep(i.com, size.com))
  }
  nN <- length(V(g.all)) #number of nodes
  pool <- matrix(rep(1:nN, prob.recon*nN), nrow = nN) #each column is 1:nN
  id2 <- apply(pool, 2, fun.sample2)
  g.all <- add_edges(g.all, as.vector(id2))
  g.all.mat <- as_adjacency_matrix(g.all) 
  g.all.E <- get.edgelist(g.all)
  group.true <- make_clusters(g.all, membership = member.true)
  data.group.true <- data.frame(label = as.numeric(V(g.all)), cluster= group.true$membership)
  length.comm.true <- max(data.group.true$cluster)
  info.all.true <- NULL
  for(i in 1:length.comm.true){
    Vset.i <- data.group.true[data.group.true$cluster==i,]$label
    set.id <- which((g.all.E[,1] %in% Vset.i) |(g.all.E[,2] %in% Vset.i))
    if(length(set.id) != 1){Eset.i <- g.all.E[set.id,]}
    if(length(set.id) == 1){Eset.i <- t(g.all.E[set.id,])}
    edge.bew <- apply(Eset.i, 1, delta, com.group = group.true)
    badlink <- sum(!edge.bew)
    sub.i <- subgraph(g.all, Vset.i)
    info <- c(nN = length(V(sub.i)), nE = length(E(sub.i)), avgD = mean(degree(sub.i)), maxD = max(degree(sub.i)), 
              denE = edge_density(sub.i), tran = transitivity(sub.i), badE = badlink, nCC = count_components(sub.i))
    info.all.true <- rbind(info.all.true, info)
  }
  rownames(info.all.true) <- 1:length.comm.true
  table.info.true <- info.all.true 
  data <- as.data.frame(scale(table.info.true[, set.feature]))
  colnames(data) <- paste("c", 1:p.variable, sep = "")
  data$groupA <- 0
  data$groupB <- 0
  ##################### MS #####################
  data.x <- data[, c(1:p.work, (p.variable + 1), ncol(data))]
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
