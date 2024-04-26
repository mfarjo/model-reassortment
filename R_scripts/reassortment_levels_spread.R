library(tidyverse)
library(compiler)
library(beepr)

#generate random normally distributed fitness values
#and generate random normally distributed small fitness offsets
set.seed(1)
fit.dist <- rnorm(100, mean = 1 , sd = 0.1)
off.dist <- rnorm(100, mean = 0.02, sd = 0.001)

#build starting genotype matrices for a virus with genes [A,B,C] 
#individual fitnesses of [A,B,C] for each genotype are drawn from fitness distribution above
#for each genotype [A,B,C], individual gene segments are "balanced" in their activity
#in that they are all within a small distribution of "offset" values 
gen.count <- 10 #number of starting genotypes
gen.mat <- matrix(data = NA, nrow = 3, ncol = gen.count)
rownames(gen.mat) <- c("A","B","C")
colnames(gen.mat) <- c(1:gen.count)
gen.mat[1,] <- sample(fit.dist, gen.count)
gen.mat[2,] <- gen.mat[1,] + (sample(c(1,-1), gen.count, replace = T)) * sample(off.dist,gen.count)
gen.mat[3,] <- gen.mat[1,] + (sample(c(1,-1), gen.count, replace = T)) * sample(off.dist,gen.count)
#save(gen.mat, file = "final_project/gen_mat.RData")


#generate starting viral populations (1 population per genotype)
#each population starts as a homogeneous group (of p individuals) with the same genotype
set.seed(1)
pops <- round(rnorm(100, mean = 11, sd = 1))
build.pop <- function(gen, psize){
  gen <- as.data.frame(gen)
  colnames(gen) <- "x"
  p <- sample(psize, 1)
  population <- rep(list(gen), p)
  names(population) <- paste0("i",c(1:p))
  return(population)
}
set.seed(1)
pop.list <- list()
for(i in seq_along(colnames(gen.mat))){
  pop.list[[i]] <- build.pop(gen.mat[,i], pops)
}
names(pop.list) <- paste0("pop",c(1:gen.count))

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

#create gene pools for non-self genotypes
list.other <- function(pop.rep){
  other.flat <- list()
  for(i in seq_along(pop.rep)){
    other.flat[[i]] <- pop.rep[-i]
    other.flat[[i]] <- bind_rows(flatten(other.flat[[i]]))
  }
  return(other.flat)
}

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

#mix populations together according to mu parameter
#entire self pool + [mu * size(other pools)] members of other pools

mix <- function(self, other, mixrate){ #combine self pool with other pool according to degree of mixing
  a.pool <- c(self$a, sample(other$a, (mixrate * length(other$a))))
  b.pool <- c(self$b, sample(other$b, (mixrate * length(other$b))))
  c.pool <- c(self$c, sample(other$c, (mixrate * length(other$c))))
  mix.stuff <- list("a" = a.pool, "b" = b.pool, "c" = c.pool)
  return(mix.stuff)
}

#get total subpopulation size for each subpop
size.up <- function(pop.rep){
  total.pop <- bind_rows(pop.rep)
  pop.size <- sum(total.pop[grep("A", rownames(total.pop)),]$n)
  return(pop.size)
}

#reassort according to degree of mixing
reassort <- function(mixed){
  p <- sample(pops,1)
  gen <- matrix(data = NA, nrow = 3, ncol = p)
  gen[1,] <- sample(mixed$a,p)
  gen[2,] <- sample(mixed$b,p)
  gen[3,] <- sample(mixed$c,p)
  new.list <- list()
  for(i in seq_along(gen[1,])){
    new.list[[i]] <- as.data.frame(gen[,i])
    colnames(new.list[[i]]) <- "x"
  }
  names(new.list) <- paste0("i",1:p)
  return(new.list)
}

#select next generation of subpopulations based on size
pop.select <- function(size){
  for(i in seq_along(size)){
    size[[i]] <- as.data.frame(size[[i]])
    colnames(size[[i]]) <- "n"
    size[[i]] <- size[[i]] %>%
      mutate("id" = names(size)[[i]]) %>%
      relocate(id, .before = n)
  }
  size <- bind_rows(size)
  size <- mutate(size, "p" = n/(sum(n)))
  rownames(size) <- c(1:gen.count)
  subset <- as.numeric(sample(rownames(size), size = gen.count, prob = size$p, replace = T))
  return(subset)
}

#ITERATION
#iterate over generations
iterate.fitness <- function(pops, gens, mu, cycles){
  rep.list <- list()
  means.list <- list()
  self.pools <- list()
  other.list <- list()
  other.pools <- list()
  mixed.pools <- list()
  size.list <- list()
  sub.list <- list()
  new.pops <- list()
  for(i in 1:cycles){
    rep.list[[i]] <- map2(pops,gens,rep.pop)
    self.pools[[i]] <- map(rep.list[[i]],pool.self)
    other.list[[i]] <- list.other(rep.list[[i]])
    names(other.list[[i]]) <- names(self.pools[[i]])
    other.pools[[i]] <- map(other.list[[i]], pool.other)
    mixed.pools[[i]] <- pmap(list(self.pools[[i]], other.pools[[i]], mu), mix)
    size.list[[i]] <- map(rep.list[[i]],size.up)
    sub.list[[i]] <- pop.select(size.list[[i]])
    new.pops[[i]] <- map(mixed.pools[[i]], reassort)
    new.pops[[i]] <- new.pops[[i]][sub.list[[i]]]
    names(new.pops[[i]]) <- paste0("pop",c(1:gen.count))
    pops <- new.pops[[i]]
  }
  return(rep.list)
}
iterate.fitness <- cmpfun(iterate.fitness)

pop.stats <- function(fit.list){
  fit.pops <- map(fit.list, bind_rows)
  fit.w <- map(fit.pops, select, w)
  fit.all <- bind_rows(fit.w)
  mean <- mean(fit.all$w)
  max <- max(fit.all$w)
  min <- min(fit.all$w)
  gap <- max-min
  mat <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 4))
  colnames(mat) <- c("mean","min","max","gap")
  mat[1,] <- c(mean,min,max,gap)
  return(mat)
}

run <- function(pop.list, gens, mu, cycles, simulations){
  mu.fit <- list()
  mu.stats <- list()
  for(i in 1:simulations){
    mu.fit[[i]] <- iterate.fitness(pop.list, gens, mu, cycles)
    mu.stats[[i]] <- map(mu.fit[[i]], pop.stats)
    mu.stats[[i]] <- bind_rows(mu.stats[[i]])
    mu.stats[[i]] <- mu.stats[[i]] %>%
      mutate("cycle" = 1:cycles)%>%
      relocate(cycle, .before = mean)
  }
  return(mu.stats)
}

#run simulations for different values of mu
gens <- 9
cycles <- 100
simulations <- 100

mu.0 <- run(pop.list, gens, mu = 0, cycles, simulations)
names(mu.0) <- paste0("sim",1:simulations)
save(mu.0, file = "final_project/run_041524/mu0.RData") ; beep(1)

mu.01 <- run(pop.list, gens, mu = 0.1, cycles, simulations)
names(mu.01) <- paste0("sim",1:simulations)
save(mu.01, file = "final_project/run_041524/mu01.RData") ; beep(1)

mu.03 <- run(pop.list, gens, mu = 0.3, cycles, simulations)
names(mu.03) <- paste0("sim",1:simulations)
save(mu.03, file = "final_project/run_041524/mu03.RData") ; beep(1)

mu.05 <- run(pop.list, gens, mu = 0.5, cycles, simulations)
names(mu.05) <- paste0("sim",1:simulations)
save(mu.05, file = "final_project/run_041524/mu05.RData") ; beep(1)

mu.07 <- run(pop.list, gens, mu = 0.7, cycles, simulations)
names(mu.07) <- paste0("sim",1:simulations)
save(mu.07, file = "final_project/run_041524/mu07.RData") ; beep(1)

mu.09 <- run(pop.list, gens, mu = 0.9, cycles, simulations)
names(mu.09) <- paste0("sim",1:simulations)
save(mu.09, file = "final_project/run_041524/mu09.RData") ; beep(1)

mu.1 <- run(pop.list, gens, mu = 1, cycles, simulations)
names(mu.1) <- paste0("sim",1:simulations)
save(mu.1, file = "final_project/run_041524/mu1.RData") ; beep(1)


#maybe i can get a variance stat for the range of fitnesses between subpopulations
#or plot the max fitness or min fitness instead of the mean
#or plot max/min  or max-min as a measure of spread?

#time to max fitness?
#maybe theres a way to approximate the slope of the timecourse before leveling out

#i might need to make the assumption that starting populations are all relatively well-balanced

#can i add a way to simulate immune memory... ie maybe better to have a C segment that is
#mismatched or something?


