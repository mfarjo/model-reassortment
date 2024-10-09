library(GillespieSSA)
library(purrr)

#load in starting genotype data
load("runs/run_100124/rdata/gen_lists_diffuse.RData")

#define general parameters
pop.init <- c(5,5) #number of initiating genomes from genotypes 1 and 2
mu.list <- c(0,0.1,0.5,0.9,1) #possible mixing rates
tf <- 6 #final time

#define replication parameters
build.rep.parms <- function(gen.list, mu){
  rep.parms <- list()
  for(i in seq_along(gen.list)){
    rep.parms[[i]] <- c(x.a1 = gen.list[[i]][[1]][1,1],
                        x.a2 = gen.list[[i]][[2]][1,1],
                        x.b1 = gen.list[[i]][[1]][2,1],
                        x.b2 = gen.list[[i]][[2]][2,1],
                        x.c1 = gen.list[[i]][[1]][3,1],
                        x.c2 = gen.list[[i]][[2]][3,1],
                        r.a1 = gen.list[[i]][[1]][1,2],
                        r.a2 = gen.list[[i]][[2]][1,2],
                        r.b1 = gen.list[[i]][[1]][2,2],
                        r.b2 = gen.list[[i]][[2]][2,2],
                        r.c1 = gen.list[[i]][[1]][3,2],
                        r.c2 = gen.list[[i]][[2]][3,2],
                        mu = mu)
  }
  return(rep.parms)
}
rep.parms <- list("mu0" = NA, "mu01" = NA, "mu05" = NA, "mu09" = NA, 
                  "mu1" = NA)
for(i in seq_along(mu.list)){
  rep.parms[[i]] <- map2(gen.lists.diff, mu.list[[i]], build.rep.parms)
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
replication.simulation <- function(parms){
  out <- ssa(
    x0 = segment.counts,
    a = rep.props,
    nu = rep.effects,
    parms = parms,
    tf = tf,
    method = ssa.otl(epsilon = 0.01),
    verbose = FALSE
  ) 
  return(out)
}

rep.per.scheme <- function(parms.scheme){
  rep.data <- list("abc" = NA, "ab" = NA, "c" = NA)
  for(i in seq_along(parms.scheme)){
    rep.data[[i]] <- map(parms.scheme[[i]], replication.simulation)
    rep.data[[i]] <- map(rep.data[[i]], function(x) x$data)
    rep.data[[i]] <- map(rep.data[[i]], as.data.frame)
  }
  return(rep.data)
}

set.seed(1)
rep.stats.init5.diff <- map(rep.parms, rep.per.scheme)
save(rep.stats.init5.diff, file = "runs/run_100124/rdata/rep_stats_init5_diff.RData")

#genome packaging based off the last cycle of replication
load("runs/run_100124/rdata/rep_stats_init5_diff.RData")

#define general parameters
mu.list <- c(0,0.1,0.5,0.9,1) #possible mixing rates

#define packaging parameters
pack.parms <- list("mu0" = NA, "mu01" = NA, "mu05" = NA, "mu09" = NA, "mu1" = NA)
for(i in seq_along(rep.stats.init5.diff)){
  pack.parms[[i]] <- c("mu" = mu.list[[i]])
}

#define initial segment and genotype counts
genotype.counts <- rep(0,8)
names(genotype.counts) <- c("A1B1C1","A1B1C2","A1B2C1","A1B2C2",
                            "A2B1C1","A2B1C2","A2B2C1","A2B2C2")

build.counts <- function(repx){
  repx.counts <- list()
  for(i in seq_along(repx)){
    repx.counts[[i]] <- c(A1 = repx[[i]][dim(repx[[i]])[1],]$A1,
                          A2 = repx[[i]][dim(repx[[i]])[1],]$A2,
                          B1 = repx[[i]][dim(repx[[i]])[1],]$B1,
                          B2 = repx[[i]][dim(repx[[i]])[1],]$B2,
                          C1 = repx[[i]][dim(repx[[i]])[1],]$C1,
                          C2 = repx[[i]][dim(repx[[i]])[1],]$C2,
                          genotype.counts)
  }
  return(repx.counts)
}

total.counts <- list("mu0" = NA, "mu01" = NA, "mu05" = NA, "mu09" = NA, 
                     "mu1" = NA)
for(i in seq_along(rep.stats.init5.diff)){
  total.counts[[i]] <- map(rep.stats.init5.diff[[i]], build.counts)
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
  "(A1/(A1+A2)) * (B1/(B1+B2)) * (C1/(C1+C2)) * sum(A1,B1,C1)", #A1B1C1
  "(A1/(A1+A2)) * (B1/(B1+B2)) * (C2/(C1+C2)) * sum(A1,B1,C2) * mu", #A1B1C2
  "(A1/(A1+A2)) * (B2/(B1+B2)) * (C1/(C1+C2)) * sum(A1,B2,C1) * mu", #A1B2C1
  "(A1/(A1+A2)) * (B2/(B1+B2)) * (C2/(C1+C2)) * sum(A1,B2,C2) * mu", #A1B2C2
  "(A2/(A1+A2)) * (B1/(B1+B2)) * (C1/(C1+C2)) * sum(A2,B1,C1) * mu", #A2B1C1
  "(A2/(A1+A2)) * (B1/(B1+B2)) * (C2/(C1+C2)) * sum(A2,B1,C2) * mu", #A2B1C2
  "(A2/(A1+A2)) * (B2/(B1+B2)) * (C1/(C1+C2)) * sum(A2,B2,C1) * mu", #A2B2C1
  "(A2/(A1+A2)) * (B2/(B1+B2)) * (C2/(C1+C2)) * sum(A2,B2,C2)"  #A2B2C2
)

#packaging simulations
pack.mu0 <- function(counts){
  out <- list()
  for(i in seq_along(counts)){
    out[[i]] <- ssa(
      x0 = counts[[i]],
      a = pack.props,
      nu = pack.effects,
      parms = pack.parms$mu0,
      tf = tf,
      method = ssa.d(),
      verbose = FALSE
    )
  }
  out <- map(out, function(x) as.data.frame(x$data))
  return(out)
}

pack.mu01 <- function(counts){
  out <- list()
  for(i in seq_along(counts)){
    out[[i]] <- ssa(
      x0 = counts[[i]],
      a = pack.props,
      nu = pack.effects,
      parms = pack.parms$mu01,
      tf = tf,
      method = ssa.d(),
      verbose = FALSE
    )
  }
  out <- map(out, function(x) as.data.frame(x$data))
  return(out)
}

pack.mu05 <- function(counts){
  out <- list()
  for(i in seq_along(counts)){
    out[[i]] <- ssa(
      x0 = counts[[i]],
      a = pack.props,
      nu = pack.effects,
      parms = pack.parms$mu05,
      tf = tf,
      method = ssa.d(),
      verbose = FALSE
    )
  }
  out <- map(out, function(x) as.data.frame(x$data))
  return(out)
}

pack.mu09 <- function(counts){
  out <- list()
  for(i in seq_along(counts)){
    out[[i]] <- ssa(
      x0 = counts[[i]],
      a = pack.props,
      nu = pack.effects,
      parms = pack.parms$mu09,
      tf = tf,
      method = ssa.d(),
      verbose = FALSE
    )
  }
  out <- map(out, function(x) as.data.frame(x$data))
  return(out)
}

pack.mu1 <- function(counts){
  out <- list()
  for(i in seq_along(counts)){
    out[[i]] <- ssa(
      x0 = counts[[i]],
      a = pack.props,
      nu = pack.effects,
      parms = pack.parms$mu1,
      tf = tf,
      method = ssa.d(),
      verbose = FALSE
    )
  }
  out <- map(out, function(x) as.data.frame(x$data))
  return(out)
}


#running
tf <- 5
set.seed(1)
pack.stats.mu0.init5 <- map(total.counts$mu0, pack.mu0)
set.seed(1)
pack.stats.mu01.init5 <- map(total.counts$mu01, pack.mu01)
set.seed(1)
pack.stats.mu05.init5 <- map(total.counts$mu05, pack.mu05)
set.seed(1)
pack.stats.mu09.init5 <- map(total.counts$mu09, pack.mu09)
set.seed(1)
pack.stats.mu1.init5 <- map(total.counts$mu1, pack.mu1)

pack.stats.init5.diff <- list("mu0" = pack.stats.mu0.init5, "mu01" = pack.stats.mu01.init5,
                              "mu05" = pack.stats.mu05.init5, "mu09" = pack.stats.mu09.init5,
                              "mu1" = pack.stats.mu1.init5)

save(pack.stats.init5.diff, file = "runs/run_100124/rdata/pack_stats_init5_diff.RData")


