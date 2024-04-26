library(tidyverse)
library(compiler)

#generate random normally distributed fitness values
#and generate random normally distributed small fitness offsets
fit.dist <- rnorm(100, mean = 1 , sd = 0.1)
off.dist <- rnorm(100, mean = 0.02, sd = 0.001)

#build starting genotype matrices for a virus with genes [A,B,C] 
#individual fitnesses of [A,B,C] for each genotype are drawn from fitness distribution above
#for each genotype [A,B,C], individual gene segments are "balanced" in their activity
#in that they are all within a small distribution of "offset" values 
build.gen <- function(count, fit, off){
  gen.mat <- matrix(data = NA, nrow = 3, ncol = count)
  rownames(gen.mat) <- c("A","B","C")
  colnames(gen.mat) <- c(1:count)
  gen.mat[1,] <- sample(fit, count)
  gen.mat[2,] <- gen.mat[1,] + (sample(c(1,-1), count, replace = T)) * sample(off,count)
  gen.mat[3,] <- gen.mat[1,] + (sample(c(1,-1), count, replace = T)) * sample(off,count)
  return(gen.mat)
}
gen.count <- 2

gen.list <- replicate(100, build.gen(gen.count, fit.dist, off.dist), simplify = F)
names(gen.list) <- paste0("pair",c(1:100))

#generate a range of possible starting population sizes (given transmission bottleneck)
pop.counts <- round(rnorm(100, mean = 2, sd = 0.5))

#generate starting viral populations (1 population per genotype)
#each population starts as a homogeneous group (of p individuals) with the same genotype
build.pop <- function(gen){
  geno <- list()
  p <- list()
  for(i in seq_along(gen[1,])){
    geno[[i]] <- as.data.frame(gen[,i])
    colnames(geno[[i]]) <- "x"
    p[[i]] <- sample(pop.counts, 1)
    geno[[i]] <- rep(list(geno[[i]]), p[[i]])
    names(geno[[i]]) <- paste0("i",c(1:p[[i]]))
  }
  names(geno) <- paste0("pop",c(1:gen.count))
  return(geno)
  
  
  #gen <- as.data.frame(gen)
  #colnames(gen) <- "x"
  #p <- sample(pop.counts, 1)
  #population <- rep(list(gen), p)
  #names(population) <- paste0("i",c(1:p))
  
}
pop.list <- map(gen.list, build.pop)

####Everything below gets iterated###
#for unit testing
pop <- pop.list[[1]]

#calculate fitness per genotype and replicate accordingly
w.calc <- function(gen){
  w <- mean(c((gen[1,])* (1-abs(gen[1,] - gen[2,])) , #fitness of A
              (gen[2,]) * (1-abs(gen[1,] - gen[2,])) , #fitness of B
              (gen[3,]))) #fitness of C
  return(w)
}

rep.pop <- function(pop,generations){
  w.list <- list()
  for(i in seq_along(pop)){
    w.list[[i]] <- matrix(data = NA, nrow = 3, ncol = 3)
    rownames(w.list[[i]]) <- c("A", "B", "C")
    colnames(w.list[[i]]) <- c("x", "w", "n")
    w.list[[i]][,1] <- pop[[i]]$x
    w.list[[i]][,2] <- w.calc(pop[[i]])
    w.list[[i]][,3] <- round(w.list[[i]][,2] * 2^generations)
  }
  names(w.list) <- names(pop)
  w.list <- lapply(w.list, as.data.frame)
  return(w.list)
}

rep.list <- map2(pop, 12, rep.pop)

#store population size and fitness
n.pull <- function(pop.rep){
  pop.n <- map(pop.rep, select, n)
  pop.n <- unlist(map(pop.n, slice, 1))
  pop.tot <- sum(pop.n)
  return(pop.tot)
}

w.pull <- function(pop.rep){
  pop.w <- map(pop.rep, select, w)
  pop.w <- unlist(map(pop.w, slice, 1))
  #pop.tot <- sum(pop.n)
  return(pop.w)
}

n.list <- map(rep.list, n.pull)

#create pools for self genotypes
pool.self <- function(pop.rep){
  a.self <- list()
  b.self <- list()
  c.self <- list()
  for(i in seq_along(pop.rep)){
    a.self[[i]] <- rep(pop.rep[[i]]$x[1], pop.rep[[i]]$n[1])
    b.self[[i]] <- rep(pop.rep[[i]]$x[2], pop.rep[[i]]$n[2])
    c.self[[i]] <- rep(pop.rep[[i]]$x[3], pop.rep[[i]]$n[3])
  }
  a.self <- unlist(a.self)
  b.self <- unlist(b.self)
  c.self <- unlist(c.self)
  self.stuff <- list("a" = a.self, "b" = b.self, "c" = c.self)
  return(self.stuff)
}

self.pools <- map(rep.list, pool.self)

#create gene pools for non-self genotypes
list.other <- function(pop.rep){
  other.flat <- list()
  for(i in seq_along(pop.rep)){
    other.flat[[i]] <- pop.rep[-i]
    other.flat[[i]] <- bind_rows(flatten(other.flat[[i]]))
  }
  return(other.flat)
}
other.list <- list.other(rep.list)

pool.other <- function(other){
  a.table <- other[grep("A", rownames(other)),]
  b.table <- other[grep("B", rownames(other)),]
  c.table <- other[grep("C", rownames(other)),]
  
  a.pool <- list()
  for(i in seq_along(a.table$x)){
    a.pool[[i]] <- rep(a.table[i,]$x, a.table[i,]$n)
  }
  a.pool <- unlist(a.pool)
  
  b.pool <- list()
  for(i in seq_along(b.table$x)){
    b.pool[[i]] <- rep(b.table[i,]$x, b.table[i,]$n)
  }
  b.pool <- unlist(b.pool)
  
  c.pool <- list()
  for(i in seq_along(c.table$x)){
    c.pool[[i]] <- rep(c.table[i,]$x, c.table[i,]$n)
  }
  c.pool <- unlist(c.pool)
  
  other.stuff <- list("a" = a.pool, "b" = b.pool, "c" = c.pool)
  return(other.stuff)
}
other.pools <- map(other.list, pool.other)

#reassort according to degree of mixing
#reassortment draws from entire self pool + [mu * size(other pools)] members of other pools

mix <- function(self, other, mixrate){ #combine self pool with other pool according to degree of mixing
  a.pool <- c(self$a, sample(other$a, (mixrate * length(other$a))))
  b.pool <- c(self$b, sample(other$b, (mixrate * length(other$b))))
  c.pool <- c(self$c, sample(other$c, (mixrate * length(other$c))))
  mix.stuff <- list("a" = a.pool, "b" = b.pool, "c" = c.pool)
  return(mix.stuff)
}
mixed.pools <- pmap(list(self.pools, other.pools, 0.3), mix)

reassort <- function(mixed){
  p <- sample(pop.counts,1)
  gen <- matrix(data = NA, nrow = 3, ncol = p)
  gen[1,] <- sample(mixed$a,p)
  gen[2,] <- sample(mixed$b,p)
  gen[3,] <- sample(mixed$c,p)
  rownames(gen) <- c("A","B","C")
  new.list <- list()
  for(i in seq_along(gen[1,])){
    new.list[[i]] <- as.data.frame(gen[,i])
    colnames(new.list[[i]]) <- "x"
  }
  names(new.list) <- paste0("i",1:p)
  return(new.list)
}
new.pops <- map(mixed.pools, reassort)

###NOW ITERATE###
iterate.fitness <- function(pops, gens, mu, cycles){
  rep.list <- list()
  n.list <- list()
  self.pools <- list()
  other.list <- list()
  other.pools <- list()
  mixed.pools <- list()
  new.pops <- list()
  for(i in 1:cycles){
    rep.list[[i]] <- map2(pops,gens,rep.pop)
    n.list[[i]] <- map(rep.list[[i]], n.pull)
    self.pools[[i]] <- map(rep.list[[i]],pool.self)
    other.list[[i]] <- list.other(rep.list[[i]])
    names(other.list[[i]]) <- names(self.pools[[i]])
    other.pools[[i]] <- map(other.list[[i]], pool.other)
    mixed.pools[[i]] <- pmap(list(self.pools[[i]], other.pools[[i]], mu), mix)
    new.pops[[i]] <- map(mixed.pools[[i]], reassort)
    pops <- new.pops[[i]]
  }
  return(rep.list)
}
iterate.fitness <- cmpfun(iterate.fitness)

gens <- 10
mu <- 1
cycles <- 10

test <- iterate.fitness(pop, gens,mu, cycles)

pop.stats <- function(fit.list){
  fit.w <- map(fit.list, w.pull)
  fit.n <- map(fit.list, n.pull)
  mean.w <- mean(unlist(fit.w)) #mean of metapopulation fitness
  max.w <- max(unlist(fit.w)) #max fitness from any given individual
  w.1 <- mean(fit.w$pop1) #avg fitness of pop 1
  w.2 <- mean(fit.w$pop2) #avg fitness of pop 2
  n.1 <- fit.n$pop1 #population growth for pop 1
  n.2 <- fit.n$pop2 #population growth for pop 2
  mat <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(mat) <- c("mean_w","max_w","mean_w_pop1","mean_w_pop2","n_pop1","n_pop2")
  mat[1,] <- c(mean.w,max.w,w.1,w.2,n.1,n.2)
  return(mat)
}

run <- function(pop.list, gens, mu, cycles){
  mu.fit <- list()
  mu.stats <- list()
  for(i in seq_along(pop.list)){
    mu.fit[[i]] <- iterate.fitness(pop.list[[i]], gens, mu, cycles)
    mu.stats[[i]] <- map(mu.fit[[i]], pop.stats)
    mu.stats[[i]] <- bind_rows(mu.stats[[i]])
    mu.stats[[i]] <- mu.stats[[i]] %>%
      mutate("cycle" = 1:cycles)%>%
      relocate(cycle, .before = mean_w) %>%
      mutate("n_pop1_track" = cumsum(n_pop1))%>%
      mutate("n_pop2_track" = cumsum(n_pop2)) %>%
      mutate("n_track" = n_pop1_track + n_pop2_track)
  }
  return(mu.stats)
}

#run simulations for different values of mu
gens <- 12
cycles <- 100


mu.0 <- run(pop.list, gens, mu = 0, cycles)
names(mu.0) <- paste0("pair",1:length(pop.list))
save(mu.0, file = "final_project/run_042524/mu0.RData")

mu.01 <- run(pop.list, gens, mu = 0.1, cycles)
names(mu.01) <- paste0("pair",1:length(pop.list))
save(mu.01, file = "final_project/run_042524/mu01.RData")

mu.03 <- run(pop.list, gens, mu = 0.3, cycles)
names(mu.03) <- paste0("pair",1:length(pop.list))
save(mu.03, file = "final_project/run_042524/mu03.RData")

mu.05 <- run(pop.list, gens, mu = 0.5, cycles)
names(mu.05) <- paste0("pair",1:length(pop.list))
save(mu.05, file = "final_project/run_042524/mu05.RData")

mu.07 <- run(pop.list, gens, mu = 0.7, cycles)
names(mu.07) <- paste0("pair",1:length(pop.list))
save(mu.07, file = "final_project/run_042524/mu07.RData")

mu.09 <- run(pop.list, gens, mu = 0.9, cycles)
names(mu.09) <- paste0("pair",1:length(pop.list))
save(mu.09, file = "final_project/run_042524/mu09.RData")

mu.1 <- run(pop.list, gens, mu = 1, cycles)
names(mu.1) <- paste0("pair",1:length(pop.list))
save(mu.1, file = "final_project/run_042524/mu1.RData")







