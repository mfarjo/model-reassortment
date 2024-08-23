library(GillespieSSA)

#define parameters
pop.init <- c(10,10) #number of initiating genomes from genotypes 1 and 2
mu <- 0.5 #mixing rate 
tf <- 9 #final time

#per-segment activity for each genotype
load("runs/run_080224/gen_list.RData")
parm.list <- list()
for(i in seq_along(gen.list)){
  parm.list[[i]] <- c(x.a1 = gen.list[[i]][1,1],
                      x.a2 = gen.list[[i]][1,2],
                      x.b1 = gen.list[[i]][2,1],
                      x.b2 = gen.list[[i]][2,2],
                      x.c1 = gen.list[[i]][3,1],
                      x.c2 = gen.list[[i]][3,2])
}

#initial counts
segment.counts <- rep(c(pop.init[1], pop.init[2]), 3)
names(segment.counts) <- c("A1","A2","B1","B2","C1","C2")
genotype.counts <- c(pop.init[1], rep(0,6), pop.init[2])
names(genotype.counts) <- c("A1B1C1","A1B1C2","A1B2C1","A1B2C2","A2B1C1","A2B1C2",
                            "A2B2C1","A2B2C2")
total.counts <- c(segment.counts, genotype.counts)

#state change matrix
effects <- as.data.frame(matrix(nrow = 14, ncol = 14))
colnames(effects) <- c("A1","A2","B1","B2","C1","C2",
                       "A1B1C1","A1B1C2","A1B2C1","A1B2C2","A2B1C1","A2B1C2","A2B2C1","A2B2C2") 
rownames(effects) <- c("plus.a1","plus.a2","plus.b1","plus.b2","plus.c1","plus.c2",
                       "plus.a1b1c1","plus.a1b1c2","plus.a1b2c1","plus.a1b2c2","plus.a2b1c1",
                       "plus.a2b1c2","plus.a2b2c1","plus.a2b2c2")
effects[1,] <- c(1,0,0,0,0,0,0,0,0,0,0,0,0,0)
effects[2,] <- c(0,1,0,0,0,0,0,0,0,0,0,0,0,0)
effects[3,] <- c(0,0,1,0,0,0,0,0,0,0,0,0,0,0)
effects[4,] <- c(0,0,0,1,0,0,0,0,0,0,0,0,0,0)
effects[5,] <- c(0,0,0,0,1,0,0,0,0,0,0,0,0,0)
effects[6,] <- c(0,0,0,0,0,1,0,0,0,0,0,0,0,0)
effects[7,] <- c(-1,0,-1,0,-1,0,1,0,0,0,0,0,0,0)
effects[8,] <- c(-1,0,-1,0,0,-1,0,1,0,0,0,0,0,0)
effects[9,] <- c(-1,0,0,-1,-1,0,0,0,1,0,0,0,0,0)
effects[10,] <- c(-1,0,0,-1,0,-1,0,0,0,1,0,0,0,0)
effects[11,] <- c(0,-1,-1,0,-1,0,0,0,0,0,1,0,0,0)
effects[12,] <- c(0,-1,-1,0,0,-1,0,0,0,0,0,1,0,0)
effects[13,] <- c(0,-1,0,-1,-1,0,0,0,0,0,0,0,1,0)
effects[14,] <- c(0,-1,0,-1,0,-1,0,0,0,0,0,0,0,1)

effects <- t(effects)

#propensities
#need to decide whether segment replication propensity also depends on current population freq
#(multiply by eg (A1/(A1+A2)))
#or maybe propensities should increase with total number in system?
props <- c(
  "(((1 - (mu * (B2/(B1+B2)))) * (1-abs(x.a1 - x.b1)) * x.a1) +
    ((mu * (B2/(B1+B2))) * (1-abs(x.a1 - x.b2)) * x.a1)) * A1", #A1
  "(((1 - (mu * (B1/(B1+B2)))) * (1-abs(x.a2 - x.b2)) * x.a2) +
    ((mu * (B1/(B1+B2))) * (1-abs(x.a2 - x.b1)) * x.a2)) * A2", #A2
  "(((1 - (mu * (A2/(A1+A2)))) * (1-abs(x.b1 - x.a1)) * x.b1) +
    ((mu * (A2/(A1+A2))) * (1-abs(x.b1 - x.a2)) * x.b1)) * B1", #B1
  "(((1 - (mu * (A1/(A1+A2)))) * (1-abs(x.b2 - x.a2)) * x.b2) +
    ((mu * (A1/(A1+A2))) * (1-abs(x.b2 - x.a1)) * x.b2)) * B2", #B2
  "x.c1 * C1", #C1
  "x.c2 * C2", #C2
  "(A1/(A1+A2)) * (B1/(B1+B2)) * (C1/(C1+C2)) * min(A1,B1,C1)", #A1B1C1
  "(A1/(A1+A2)) * (B1/(B1+B2)) * (C2/(C1+C2)) * min(A1,B1,C2) * mu", #A1B1C2
  "(A1/(A1+A2)) * (B2/(B1+B2)) * (C1/(C1+C2)) * min(A1,B2,C1) * mu", #A1B2C1
  "(A1/(A1+A2)) * (B2/(B1+B2)) * (C2/(C1+C2)) * min(A1,B2,C2) * mu", #A1B2C2
  "(A2/(A1+A2)) * (B1/(B1+B2)) * (C1/(C1+C2)) * min(A2,B1,C1) * mu", #A2B1C1
  "(A2/(A1+A2)) * (B1/(B1+B2)) * (C2/(C1+C2)) * min(A2,B1,C2) * mu", #A2B1C2
  "(A2/(A1+A2)) * (B2/(B1+B2)) * (C1/(C1+C2)) * min(A2,B2,C1) * mu", #A2B2C1
  "(A2/(A1+A2)) * (B2/(B1+B2)) * (C2/(C1+C2)) * min(A2,B2,C2)"  #A2B2C2
)

#test
#out <- ssa(
  #x0 = total.counts,
  #a = props,
  #nu = effects,
  #parms = parm.list[[3]],
  #tf = tf,
  #method = ssa.otl(),
  #verbose = FALSE
#) 

#results <- as.data.frame(out$data)


#simulations
#try tau leaping
set.seed(1)
replication.simulation <- function(parms){
  out <- ssa(
    x0 = total.counts,
    a = props,
    nu = effects,
    parms = parms,
    tf = tf,
    method = ssa.otl(),
    verbose = FALSE
  ) 
  return(out)
}

results.mu05 <- map(parm.list, replication.simulation)
#this takes about 10 min which isnt good but is maybe better

results.mu05 <- map(results.mu05, function(x) x$data)
results.mu05 <- map(results.mu05, as.data.frame)
beep(2)
save(results.mu05, file = "runs/run_080224/mu05.RData")

#gen.results <- result.list[[5]][,c(1,8:15)] #this is just a test
#gen.results <- results[,c(1,8:15)]
#gen.gather <- gather(gen.results, key = "genotype", value = "count",-t)

#ggplot(data = gen.gather, aes(x = t, y = count, color = genotype))+
  #geom_line()
 





