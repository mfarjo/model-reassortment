library(tidyverse)

set.seed(1)
#generate random normally distributed and gamma distributed fitness values
#and generate random normally distributed small fitness offsets
fit.dist.norm <- rnorm(100, mean = 1 , sd = 0.1)
hist(fit.dist.norm)
fit.dist.plus <- rlnorm(100, meanlog = 0.1, sdlog =  0.25)
hist(fit.dist.plus)
off.dist <- rnorm(100, mean = 0.02, sd = 0.001)

#generate relatedness spectrum
r1 <- runif(25, 0, 0.1)
r2 <- runif(25, 0.2,0.3)
r3 <- runif(25, 0.7,0.8)
r4 <- runif(25, 0.9,1)
r <- c(r1,r2,r3,r4)
hist(r, breaks = seq(0,1,0.01))

#build genotypes
#want relatedness score in a separate column

build.gen.abc <- function(count,fit,off){
  gen.mat <- list()
  for(i in 1:count){
    gen.mat[[i]] <- matrix(data = NA, nrow = 3, ncol = 2)
    rownames(gen.mat[[i]]) <- c("A","B","C")
    colnames(gen.mat[[i]]) <- c("x","r")
    gen.mat[[i]][,1] <- sample(fit, 3)
    gen.mat[[i]][1,2] <- sample(r,1)
    gen.mat[[i]][2,2] <- gen.mat[[i]][1,2] + (sample(c(1,-1), 1, replace = T)) * sample(off,1)
    gen.mat[[i]][3,2] <- gen.mat[[i]][1,2] + (sample(c(1,-1), 1, replace = T)) * sample(off,1)
  }
  names(gen.mat) <- c("i1","i2")
  return(gen.mat)
}

build.gen.ab <- function(count,fit,fit.plus,off){
  gen.mat <- list()
  for(i in 1:count){
    gen.mat[[i]] <- matrix(data = NA, nrow = 3, ncol = 2)
    rownames(gen.mat[[i]]) <- c("A","B","C")
    colnames(gen.mat[[i]]) <- c("x","r")
    gen.mat[[i]][1,1] <- sample(fit.plus, 1)
    gen.mat[[i]][2,1] <- sample(fit.plus, 1)
    gen.mat[[i]][3,1] <- sample(fit, 1)
    gen.mat[[i]][1,2] <- sample(r,1)
    gen.mat[[i]][2,2] <- gen.mat[[i]][1,2] + (sample(c(1,-1), 1, replace = T)) * sample(off,1)
    gen.mat[[i]][3,2] <- gen.mat[[i]][1,2] + (sample(c(1,-1), 1, replace = T)) * sample(off,1)
  }
  names(gen.mat) <- c("i1","i2")
  return(gen.mat)
}

build.gen.c <- function(count,fit,fit.plus,off){
  gen.mat <- list()
  for(i in 1:count){
    gen.mat[[i]] <- matrix(data = NA, nrow = 3, ncol = 2)
    rownames(gen.mat[[i]]) <- c("A","B","C")
    colnames(gen.mat[[i]]) <- c("x","r")
    gen.mat[[i]][1,1] <- sample(fit, 1)
    gen.mat[[i]][2,1] <- sample(fit, 1)
    gen.mat[[i]][3,1] <- sample(fit.plus, 1)
    gen.mat[[i]][1,2] <- sample(r,1)
    gen.mat[[i]][2,2] <- gen.mat[[i]][1,2] + (sample(c(1,-1), 1, replace = T)) * sample(off,1)
    gen.mat[[i]][3,2] <- gen.mat[[i]][1,2] + (sample(c(1,-1), 1, replace = T)) * sample(off,1)
  }
  names(gen.mat) <- c("i1","i2")
  return(gen.mat)
}

gen.count <- 2

gen.list.abc <- replicate(100, build.gen.abc(gen.count, fit.dist.norm, off.dist), 
                          simplify = F)
save(gen.list.abc, file = "runs/run_091024/gen_list_abc.RData")
gen.list.ab <- replicate(100, build.gen.ab(gen.count, fit.dist.norm, fit.dist.plus, off.dist), 
                         simplify = F)
save(gen.list.ab, file = "runs/run_091024/gen_list_ab.RData")
gen.list.c <- replicate(100, build.gen.c(gen.count, fit.dist.norm, fit.dist.plus, off.dist), 
                         simplify = F)
save(gen.list.c, file = "runs/run_091024/gen_list_c.RData")

gen.lists <- list("abc" = gen.list.abc, "ab" = gen.list.ab, "c" = gen.list.c)
save(gen.lists, file = "runs/run_091024/gen_lists.RData")

#generate list of possible genotype combinations
combo <- as.data.frame(matrix(data = c("A1","A2","B1","B2","C1","C2"), ncol = 3, nrow = 2))
colnames(combo) <- c("A","B","C")
combo <- expand(combo,A,B,C)
combo <- split(combo,1:nrow(combo))
combo <- unlist(lapply(combo, paste0, collapse = ""))
save(combo,file = "runs/run_091024/combo.RData")










