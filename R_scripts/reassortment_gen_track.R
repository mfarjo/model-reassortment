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
#save(gen.list, file = "final_project/run_051324/gen_list.RData")

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
}
pop.list <- map(gen.list, build.pop)
#save(pop.list, file = "final_project/run_051324/pop_list.RData")

#choose relative weights for each gene segment (should sum to 1)
fit.vec <- c(0.3333333, 0.3333333, 0.3333333)

#this gives a vector of possible genotype combos
combo <- as.data.frame(matrix(data = c("A1","A2","B1","B2","C1","C2"), ncol = 3, nrow = 2))
colnames(combo) <- c("A","B","C")
combo <- expand(combo,A,B,C)
combo <- split(combo,1:nrow(combo))
combo <- unlist(lapply(combo, paste0, collapse = ""))

#for unit testing:
load("final_project/run_051324/pop_list.RData")
load("final_project/run_051324/gen_list.RData")
pop <- pop.list[[1]]
geno <- gen.list[[1]]

#ITERATION
#assign gene IDs to each individual
assign.geno <- function(individual){
  individual <- mutate(individual,"ID" = NA)
  individual$ID[[1]] <- paste0("A",which(geno[1,] == individual$x[[1]]))
  individual$ID[[2]] <- paste0("B",which(geno[2,] == individual$x[[2]]))
  individual$ID[[3]] <- paste0("C",which(geno[3,] == individual$x[[3]]))
  return(individual)
}
assign.pop <- function(subpop){
  assignments <- map(subpop, assign.geno)
  return(assignments)
}

pop.label <- map(pop, assign.pop)

#calculate fitness per genotype and replicate accordingly
w.calc <- function(gen){
  w <- sum(c((fit.vec[1] * (gen$x[1])* (1-abs(gen$x[1] - gen$x[2]))) , #fitness of A
             (fit.vec[2] * (gen$x[2]) * (1-abs(gen$x[1] - gen$x[2]))) , #fitness of B
             (fit.vec[3] * (gen$x[3])))) #fitness of C
  return(w)
}
rep.pop <- function(pop,generations){
  w.list <- list()
  for(i in seq_along(pop)){
    w.list[[i]] <- as.data.frame(matrix(data = NA, nrow = 3, ncol = 4))
    rownames(w.list[[i]]) <- c("A", "B", "C")
    colnames(w.list[[i]]) <- c("x", "w", "n", "id")
    w.list[[i]][,1] <- pop[[i]]$x
    w.list[[i]][,2] <- w.calc(pop[[i]])
    w.list[[i]][,3] <- round(w.list[[i]][,2] * 2^generations)
    w.list[[i]][,4] <- paste0(pop[[i]]$ID, collapse = "")
  }
  names(w.list) <- names(pop)
  w.list <- lapply(w.list, as.data.frame)
  return(w.list)
}

rep.list <- map2(pop.label, 12, rep.pop)

#break down replication by genotype
gen.breakdown <- function(pop.rep){
  gen.rep <- as.data.frame(matrix(data = NA, ncol = 8, nrow = 1))
  colnames(gen.rep) <- combo
  pop.rep <- flatten(pop.rep)
  pop.rep <- map(pop.rep, select, n, id)
  pop.rep <- map(pop.rep, function(x) x[1,])
  pop.rep <- bind_rows(pop.rep)
  pop.rep$n <- as.numeric(pop.rep$n)
  pop.split <- group_split(pop.rep, id)
  pop.sum <- map(pop.split, function(x) sum(x$n))
  pop.id <- map(pop.split, function(x) x$id[1])
  names(pop.sum) <- pop.id
  pop.sum <- as.data.frame(pop.sum)
  gen.rep <- suppressMessages(left_join(pop.sum, gen.rep))
  return(gen.rep)
}

gen.track <- gen.breakdown(rep.list)

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


##ITERATE##
iterate.fitness <- function(pops, generations, mu, cycles){
  rep.list <- list()
  gen.track <- list()
  self.pools <- list()
  other.list <- list()
  other.pools <- list()
  mixed.pools <- list()
  new.pops <- list()
  for(i in 1:cycles){
    pop.label <- map(pops,assign.pop)
    rep.list[[i]] <- map2(pop.label,generations,rep.pop)
    gen.track[[i]] <- gen.breakdown(rep.list[[i]])
    self.pools[[i]] <- map(rep.list[[i]],pool.self)
    other.list[[i]] <- list.other(rep.list[[i]])
    names(other.list[[i]]) <- names(self.pools[[i]])
    other.pools[[i]] <- map(other.list[[i]], pool.other)
    mixed.pools[[i]] <- pmap(list(self.pools[[i]], other.pools[[i]], mu), mix)
    new.pops[[i]] <- map(mixed.pools[[i]], reassort)
    pops <- new.pops[[i]]
  }
  gen.track <- bind_rows(gen.track)
  gen.track[is.na(gen.track)] <- 0
  return(gen.track)
}
iterate.fitness <- cmpfun(iterate.fitness)

generations <- 12
cycles <- 50

mu.0 <- list()
for(i in seq_along(pop.list)){
  geno <- gen.list[[i]]
  mu.0[[i]] <- iterate.fitness(pop.list[[i]], generations, mu = 0, cycles)
}
names(mu.0) <- paste0("pair",c(1:100))
mu.0 <- map(mu.0, cumsum)
save(mu.0, file = "final_project/run_051324/mu_0.RData")

mu.0.props <- map(mu.0, function(x) x/rowSums(x))
save(mu.0.props, file = "final_project/run_051324/mu_0_props.RData")

mu.01 <- list()
for(i in seq_along(pop.list)){
  geno <- gen.list[[i]]
  mu.01[[i]] <- iterate.fitness(pop.list[[i]], generations, mu = 0.1, cycles)
}
names(mu.01) <- paste0("pair",c(1:100))
mu.01 <- map(mu.01, cumsum)
save(mu.01, file = "final_project/run_051324/mu_01.RData")

mu.01.props <- map(mu.01, function(x) x/rowSums(x))
save(mu.01.props, file = "final_project/run_051324/mu_01_props.RData")

mu.03 <- list()
for(i in seq_along(pop.list)){
  geno <- gen.list[[i]]
  mu.03[[i]] <- iterate.fitness(pop.list[[i]], generations, mu = 0.3, cycles)
}
names(mu.03) <- paste0("pair",c(1:100))
mu.03 <- map(mu.03, cumsum)
save(mu.03, file = "final_project/run_051324/mu_03.RData")

mu.03.props <- map(mu.03, function(x) x/rowSums(x))
save(mu.03.props, file = "final_project/run_051324/mu_03_props.RData")

mu.05 <- list()
for(i in seq_along(pop.list)){
  geno <- gen.list[[i]]
  mu.05[[i]] <- iterate.fitness(pop.list[[i]], generations, mu = 0.5, cycles)
}
names(mu.05) <- paste0("pair",c(1:100))
mu.05 <- map(mu.05, cumsum)
save(mu.05, file = "final_project/run_051324/mu_05.RData")

mu.05.props <- map(mu.05, function(x) x/rowSums(x))
save(mu.05.props, file = "final_project/run_051324/mu_05_props.RData")

mu.07 <- list()
for(i in seq_along(pop.list)){
  geno <- gen.list[[i]]
  mu.07[[i]] <- iterate.fitness(pop.list[[i]], generations, mu = 0.7, cycles)
}
names(mu.07) <- paste0("pair",c(1:100))
mu.07 <- map(mu.07, cumsum)
save(mu.07, file = "final_project/run_051324/mu_07.RData")

mu.07.props <- map(mu.07, function(x) x/rowSums(x))
save(mu.07.props, file = "final_project/run_051324/mu_07_props.RData")

mu.09 <- list()
for(i in seq_along(pop.list)){
  geno <- gen.list[[i]]
  mu.09[[i]] <- iterate.fitness(pop.list[[i]], generations, mu = 0.9, cycles)
}
names(mu.09) <- paste0("pair",c(1:100))
mu.09 <- map(mu.09, cumsum)
save(mu.09, file = "final_project/run_051324/mu_09.RData")

mu.09.props <- map(mu.09, function(x) x/rowSums(x))
save(mu.09.props, file = "final_project/run_051324/mu_09_props.RData")

mu.1 <- list()
for(i in seq_along(pop.list)){
  geno <- gen.list[[i]]
  mu.1[[i]] <- iterate.fitness(pop.list[[i]], generations, mu = 1, cycles)
}
names(mu.1) <- paste0("pair",c(1:100))
mu.1 <- map(mu.1, cumsum)
save(mu.1, file = "final_project/run_051324/mu_1.RData")

mu.1.props <- map(mu.1, function(x) x/rowSums(x))
save(mu.1.props, file = "final_project/run_051324/mu_1_props.RData")



#SOME FITNESS CALCULATIONS#

#fitness calculations for all possible genotypes:
gen.expand <- function(geno){
  combo.mat <- as.data.frame(matrix(data = c(geno[1,1],geno[1,2],geno[2,1],geno[2,2],
                                             geno[3,1],geno[3,2]), ncol = 3, nrow = 2))
  colnames(combo.mat) <- c("A","B","C")
  combo.ex <- expand(combo.mat,A,B,C)
  combo.ex <- combo.ex %>%
    mutate("id" = combo) %>%
    relocate(id, .before = A)
  
  w <- list()
  for(i in seq_along(combo.ex$id)){
    w[[i]] <- sum((fit.vec[1]) * ((combo.ex[i,]$A) * (1-abs(combo.ex[i,]$A - combo.ex[i,]$B))),
                  (fit.vec[2]) * ((combo.ex[i,]$B) * (1-abs(combo.ex[i,]$B - combo.ex[i,]$A))),
                  (fit.vec[3]) * (combo.ex[i,]$C))
  }
  w <- unlist(w)
  combo.ex <- mutate(combo.ex, "w" = w)
  return(combo.ex)
}
all.w <- map(gen.list,gen.expand)
all.w <- map(all.w,t)
all.w <- map(all.w, as.data.frame)
for(i in seq_along(all.w)){
  colnames(all.w[[i]]) <- all.w[[i]][1,]
}
all.w <- map(all.w, function(x) x[-1,])
for(i in seq_along(all.w)){
  all.w[[i]] <- mutate_all(all.w[[i]], function(x) as.numeric(x))
}

only.w <- map(all.w, slice, 4)
#save(only.w, file = "final_project/run_051324/only_w.RData")

#calculate mean w:
mean.w.calc <- function(w,counts){
  counts <- suppressMessages(right_join(w,counts))
  mean.w <- list()
  for(i in seq_along(counts$A1B1C1)){
    mean.w[[i]] <- (rowSums(counts[i,]*w))/(rowSums(counts[i,]))
  }
  mean.w <- unlist(mean.w)
  return(mean.w)
}

mu.0.means <- map2(only.w,mu.0,mean.w.calc)
save(mu.0.means, file = "final_project/run_051324/mu_0_means.RData")
mu.01.means <- map2(only.w,mu.01,mean.w.calc)
save(mu.01.means, file = "final_project/run_051324/mu_01_means.RData")
mu.03.means <- map2(only.w,mu.03,mean.w.calc)
save(mu.03.means, file = "final_project/run_051324/mu_03_means.RData")
mu.05.means <- map2(only.w,mu.05,mean.w.calc)
save(mu.05.means, file = "final_project/run_051324/mu_05_means.RData")
mu.07.means <- map2(only.w,mu.07,mean.w.calc)
save(mu.07.means, file = "final_project/run_051324/mu_07_means.RData")
mu.09.means <- map2(only.w,mu.09,mean.w.calc)
save(mu.09.means, file = "final_project/run_051324/mu_09_means.RData")
mu.1.means <- map2(only.w,mu.1,mean.w.calc)
save(mu.1.means, file = "final_project/run_051324/mu_1_means.RData")



