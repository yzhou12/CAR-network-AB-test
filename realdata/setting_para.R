a.list <- rep(c(0, 0.25, 0.5, 0.75, 1), 4)
p.list <- rep(rep(c(2, 4), each = 5), 2)
r.list <- rep(c(0.1, 0.5), each = 10)

data.para <- data.frame(alpha = a.list, p = p.list, r = r.list, set = 1:20)

prob.recon <- data.para$r[setting] #reconnect probability (r)
alpha0 <- data.para$alpha[setting] #network effect from A
p.work <- data.para$p[setting] # number of working covariates
