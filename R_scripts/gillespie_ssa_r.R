library(GillespieSSA)
library(purrr)

#define parameters
pop.init <- c(1,1) #number of initiating genomes from genotypes 1 and 2
mu <- 0 #mixing rate 
tf <- 8 #final time

#per-segment activity for each genotype
load("runs/run_091024/gen_list_abc.RData")
load("runs/run_091024/gen_list_ab.RData")
load("runs/run_091024/gen_list_c.RData")
rep.parms.abc <- list()
for(i in seq_along(gen.list.abc)){
  rep.parms.abc[[i]] <- c(x.a1 = gen.list.abc[[i]][[1]][1,1],
                          x.a2 = gen.list.abc[[i]][[2]][1,1],
                          x.b1 = gen.list.abc[[i]][[1]][2,1],
                          x.b2 = gen.list.abc[[i]][[2]][2,1],
                          x.c1 = gen.list.abc[[i]][[1]][3,1],
                          x.c2 = gen.list.abc[[i]][[2]][3,1],
                          r.a1 = gen.list.abc[[i]][[1]][1,2],
                          r.a2 = gen.list.abc[[i]][[2]][1,2],
                          r.b1 = gen.list.abc[[i]][[1]][2,2],
                          r.b2 = gen.list.abc[[i]][[2]][2,2],
                          r.c1 = gen.list.abc[[i]][[1]][3,2],
                          r.c2 = gen.list.abc[[i]][[2]][3,2])
}

rep.parms.ab <- list()
for(i in seq_along(gen.list.ab)){
  rep.parms.ab[[i]] <- c(x.a1 = gen.list.ab[[i]][[1]][1,1],
                         x.a2 = gen.list.ab[[i]][[2]][1,1],
                         x.b1 = gen.list.ab[[i]][[1]][2,1],
                         x.b2 = gen.list.ab[[i]][[2]][2,1],
                         x.c1 = gen.list.ab[[i]][[1]][3,1],
                         x.c2 = gen.list.ab[[i]][[2]][3,1],
                         r.a1 = gen.list.ab[[i]][[1]][1,2],
                         r.a2 = gen.list.ab[[i]][[2]][1,2],
                         r.b1 = gen.list.ab[[i]][[1]][2,2],
                         r.b2 = gen.list.ab[[i]][[2]][2,2],
                         r.c1 = gen.list.ab[[i]][[1]][3,2],
                         r.c2 = gen.list.ab[[i]][[2]][3,2])
}

rep.parms.c <- list()
for(i in seq_along(gen.list.ab)){
  rep.parms.c[[i]] <- c(x.a1 = gen.list.c[[i]][[1]][1,1],
                        x.a2 = gen.list.c[[i]][[2]][1,1],
                        x.b1 = gen.list.c[[i]][[1]][2,1],
                        x.b2 = gen.list.c[[i]][[2]][2,1],
                        x.c1 = gen.list.c[[i]][[1]][3,1],
                        x.c2 = gen.list.c[[i]][[2]][3,1],
                        r.a1 = gen.list.c[[i]][[1]][1,2],
                        r.a2 = gen.list.c[[i]][[2]][1,2],
                        r.b1 = gen.list.c[[i]][[1]][2,2],
                        r.b2 = gen.list.c[[i]][[2]][2,2],
                        r.c1 = gen.list.c[[i]][[1]][3,2],
                        r.c2 = gen.list.c[[i]][[2]][3,2])
}



#initial segment counts
segment.counts <- rep(c(pop.init[1], pop.init[2]), 3)
names(segment.counts) <- c("A1","A2","B1","B2","C1","C2")

#replication state change matrix
rep.effects <- as.data.frame(matrix(nrow = 6, ncol = 6))
colnames(rep.effects) <- c("A1","A2","B1","B2","C1","C2") 
rownames(rep.effects) <- c("plus.a1","plus.a2","plus.b1","plus.b2","plus.c1","plus.c2")
rep.effects[1,] <- c(1,0,0,0,0,0)
rep.effects[2,] <- c(0,1,0,0,0,0)
rep.effects[3,] <- c(0,0,1,0,0,0)
rep.effects[4,] <- c(0,0,0,1,0,0)
rep.effects[5,] <- c(0,0,0,0,1,0)
rep.effects[6,] <- c(0,0,0,0,0,1)

rep.effects <- t(rep.effects)

#replication propensities
rep.props <- c(
  "(((1 - (mu * (B2/(B1+B2)))) * (1-abs(r.a1 - r.b1)) * x.a1) +
    ((mu * (B2/(B1+B2))) * (1-abs(r.a1 - r.b2)) * x.a1)) * A1", #A1
  "(((1 - (mu * (B1/(B1+B2)))) * (1-abs(r.a2 - r.b2)) * x.a2) +
    ((mu * (B1/(B1+B2))) * (1-abs(r.a2 - r.b1)) * x.a2)) * A2", #A2
  "(((1 - (mu * (A2/(A1+A2)))) * (1-abs(r.b1 - r.a1)) * x.b1) +
    ((mu * (A2/(A1+A2))) * (1-abs(r.b1 - r.a2)) * x.b1)) * B1", #B1
  "(((1 - (mu * (A1/(A1+A2)))) * (1-abs(r.b2 - r.a2)) * x.b2) +
    ((mu * (A1/(A1+A2))) * (1-abs(r.b2 - r.a1)) * x.b2)) * B2", #B2
  "x.c1 * C1", #C1
  "x.c2 * C2"  #C2
)

#replication simulations
set.seed(1)
replication.simulation <- function(parms){
  out <- ssa(
    x0 = segment.counts,
    a = rep.props,
    nu = rep.effects,
    parms = parms,
    tf = tf,
    method = ssa.otl(),
    verbose = FALSE
  ) 
  return(out)
}

rep.mu0.abc <- map(rep.parms.abc, replication.simulation)
rep.mu0.abc <- map(rep.mu0.abc, function(x) x$data)
rep.mu0.abc <- map(rep.mu0.abc, as.data.frame)
save(rep.mu0.abc, file = "runs/run_091024/rep_mu0_abc.RData")

rep.mu0.ab <- map(rep.parms.ab, replication.simulation)
rep.mu0.ab <- map(rep.mu0.ab, function(x) x$data)
rep.mu0.ab <- map(rep.mu0.ab, as.data.frame)
save(rep.mu0.ab, file = "runs/run_091024/rep_mu0_ab.RData")

rep.mu0.c <- map(rep.parms.c, replication.simulation)
rep.mu0.c <- map(rep.mu0.c, function(x) x$data)
rep.mu0.c <- map(rep.mu0.c, as.data.frame)
save(rep.mu0.c, file = "runs/run_091024/rep_mu0_c.RData")


#genome packaging based off the last cycle of replication

#post-replication segment counts 
rep.abc <- rep.mu0.abc
rep.ab <- rep.mu0.ab
rep.c <- rep.mu0.c

#initial segment and genotype counts
genotype.counts <- rep(0,8)
names(genotype.counts) <- c("A1B1C1","A1B1C2","A1B2C1","A1B2C2",
                            "A2B1C1","A2B1C2","A2B2C1","A2B2C2")
total.counts.abc <- list()
for(i in seq_along(rep.abc)){
  total.counts.abc[[i]] <- c(A1 = rep.abc[[i]][dim(rep.abc[[i]])[1],]$A1,
                             A2 = rep.abc[[i]][dim(rep.abc[[i]])[1],]$A2,
                             B1 = rep.abc[[i]][dim(rep.abc[[i]])[1],]$B1,
                             B2 = rep.abc[[i]][dim(rep.abc[[i]])[1],]$B2,
                             C1 = rep.abc[[i]][dim(rep.abc[[i]])[1],]$C1,
                             C2 = rep.abc[[i]][dim(rep.abc[[i]])[1],]$C2,
                             genotype.counts)
}

total.counts.ab <- list()
for(i in seq_along(rep.ab)){
  total.counts.ab[[i]] <- c(A1 = rep.ab[[i]][dim(rep.ab[[i]])[1],]$A1,
                            A2 = rep.ab[[i]][dim(rep.ab[[i]])[1],]$A2,
                            B1 = rep.ab[[i]][dim(rep.ab[[i]])[1],]$B1,
                            B2 = rep.ab[[i]][dim(rep.ab[[i]])[1],]$B2,
                            C1 = rep.ab[[i]][dim(rep.ab[[i]])[1],]$C1,
                            C2 = rep.ab[[i]][dim(rep.ab[[i]])[1],]$C2,
                            genotype.counts)
}

total.counts.c <- list()
for(i in seq_along(rep.c)){
  total.counts.c[[i]] <- c(A1 = rep.c[[i]][dim(rep.c[[i]])[1],]$A1,
                           A2 = rep.c[[i]][dim(rep.c[[i]])[1],]$A2,
                           B1 = rep.c[[i]][dim(rep.c[[i]])[1],]$B1,
                           B2 = rep.c[[i]][dim(rep.c[[i]])[1],]$B2,
                           C1 = rep.c[[i]][dim(rep.c[[i]])[1],]$C1,
                           C2 = rep.c[[i]][dim(rep.c[[i]])[1],]$C2,
                           genotype.counts)
}

#packaging state change matrix
pack.effects <- as.data.frame(matrix(nrow = 8, ncol = 14))
colnames(pack.effects) <- c("A1","A2","B1","B2","C1","C2",
                            "A1B1C1","A1B1C2","A1B2C1","A1B2C2",
                            "A2B1C1","A2B1C2","A2B2C1","A2B2C2") 
rownames(pack.effects) <- c("plus.a1b1c1","plus.a1b1c2","plus.a1b2c1","plus.a1b2c2",
                            "plus.a2b1c1","plus.a2b1c2","plus.a2b2c1","plus.a2b2c2")
pack.effects[1,] <- c(-1,0,-1,0,-1,0,1,0,0,0,0,0,0,0)
pack.effects[2,] <- c(-1,0,-1,0,0,-1,0,1,0,0,0,0,0,0)
pack.effects[3,] <- c(-1,0,0,-1,-1,0,0,0,1,0,0,0,0,0)
pack.effects[4,] <- c(-1,0,0,-1,0,-1,0,0,0,1,0,0,0,0)
pack.effects[5,] <- c(0,-1,-1,0,-1,0,0,0,0,0,1,0,0,0)
pack.effects[6,] <- c(0,-1,-1,0,0,-1,0,0,0,0,0,1,0,0)
pack.effects[7,] <- c(0,-1,0,-1,-1,0,0,0,0,0,0,0,1,0)
pack.effects[8,] <- c(0,-1,0,-1,0,-1,0,0,0,0,0,0,0,1)

pack.effects <- t(pack.effects)


#packaging propensities
pack.props <- c(
  "(A1/(A1+A2)) * (B1/(B1+B2)) * (C1/(C1+C2)) * min(A1,B1,C1)", #A1B1C1
  "(A1/(A1+A2)) * (B1/(B1+B2)) * (C2/(C1+C2)) * min(A1,B1,C2) * mu", #A1B1C2
  "(A1/(A1+A2)) * (B2/(B1+B2)) * (C1/(C1+C2)) * min(A1,B2,C1) * mu", #A1B2C1
  "(A1/(A1+A2)) * (B2/(B1+B2)) * (C2/(C1+C2)) * min(A1,B2,C2) * mu", #A1B2C2
  "(A2/(A1+A2)) * (B1/(B1+B2)) * (C1/(C1+C2)) * min(A2,B1,C1) * mu", #A2B1C1
  "(A2/(A1+A2)) * (B1/(B1+B2)) * (C2/(C1+C2)) * min(A2,B1,C2) * mu", #A2B1C2
  "(A2/(A1+A2)) * (B2/(B1+B2)) * (C1/(C1+C2)) * min(A2,B2,C1) * mu", #A2B2C1
  "(A2/(A1+A2)) * (B2/(B1+B2)) * (C2/(C1+C2)) * min(A2,B2,C2)"  #A2B2C2
)


#packaging simulations
tf <- 10
set.seed(1)
packaging.simulation <- function(totals){
  out <- ssa(
    x0 = totals,
    a = pack.props,
    nu = pack.effects,
    parms = NULL,
    tf = tf,
    method = ssa.otl(),
    verbose = FALSE
  ) 
  return(out)
}

pack.mu0.abc <- map(total.counts.abc, packaging.simulation)
pack.mu0.abc <- map(pack.mu0.abc, function(x) x$data)
pack.mu0.abc <- map(pack.mu0.abc, as.data.frame)
save(pack.mu0.abc, file = "runs/run_091024/pack_mu0_abc.RData")

pack.mu0.ab <- map(total.counts.ab, packaging.simulation)
pack.mu0.ab <- map(pack.mu0.ab, function(x) x$data)
pack.mu0.ab <- map(pack.mu0.ab, as.data.frame)
save(pack.mu0.ab, file = "runs/run_091024/pack_mu0_ab.RData")

pack.mu0.c <- map(total.counts.c, packaging.simulation)
pack.mu0.c <- map(pack.mu0.c, function(x) x$data)
pack.mu0.c <- map(pack.mu0.c, as.data.frame)
save(pack.mu0.c, file = "runs/run_091024/pack_mu0_c.RData")


