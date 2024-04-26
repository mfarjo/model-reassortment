library(tidyverse)

#generate random normally distributed fitness values 
set.seed(1)
fit.dist <- rnorm(100, mean = 1 , sd = 0.1)
fit.tab <- as.data.frame(fit.dist)
fit.tab$fit.dist <- as.numeric(fit.tab$fit.dist)

ggplot(data = fit.tab, aes(x = fit.dist))+
  geom_density()+
  xlab("fitness")

#build genotype matrices
gen.count <- 20 #number of starting genotypes
initial.pop <- 50 #initial abundance of starting genotypes

gen.mat <- matrix(data = sample(fit.dist, (gen.count * 3)), nrow = 3, ncol = gen.count)
rownames(gen.mat) <- c("A","B","C")
colnames(gen.mat) <- c(1:gen.count)
save(gen.mat, file = "final_project/gen_mat.RData")

a.group <- as.data.frame(cbind("x" = gen.mat[1,], "n0" = initial.pop))
b.group <- as.data.frame(cbind("x" = gen.mat[2,], "n0" = initial.pop))
c.group <- as.data.frame(cbind("x" = gen.mat[3,], "n0" = initial.pop))

#calculate fitnesses
w.calc <- function(gen){
  w <- c()
  for(i in 1:dim(gen)[2]){
    w[i] <- mean(c(
      (gen[,i][1])* (1-abs(gen[,i][1] - gen[,i][2])) , #fitness of A
      (gen[,i][2]) * (1-abs(gen[,i][1] - gen[,i][2])) , #fitness of B
      (gen[,i][3]) #fitness of C
    ))
  }
  names(w) <- colnames(gen)
  return(w)
}

#replication
replicate <- function(a,b,c,w){
  a.rep <- cbind(a,w)
  a.rep <- cbind(a.rep, "n" = round(a.rep$n0 * a.rep$w))
  b.rep <- cbind(b,w)
  b.rep <- cbind(b.rep, "n" = round(b.rep$n0 * b.rep$w))
  c.rep <- cbind(c,w)
  c.rep <- cbind(c.rep, "n" = round(c.rep$n0 * c.rep$w))
  rep.list <- list("a.group" = a.rep, "b.group" = b.rep, "c.group" = c.rep)
  return(rep.list)
}

#pooling
pool <- function(reps){
  a.pool <- rep(reps$a.group$x, reps$a.group$n)
  b.pool <- rep(reps$b.group$x, reps$b.group$n)
  c.pool <- rep(reps$c.group$x, reps$c.group$n)
  pool.list <- list("a.group" = a.pool, "b.group" = b.pool, "c.group" = c.pool)
  return(pool.list)
}

#reassortment
reassort <- function(pools){
  newmat <- matrix(data = NA, nrow = 3, ncol = gen.count)
  newmat[1,] <- sample(pools$a.group, gen.count)
  newmat[2,] <- sample(pools$b.group, gen.count)
  newmat[3,] <- sample(pools$c.group, gen.count)
  rownames(newmat) <- c("A","B","C")
  colnames(newmat) <- c(1:gen.count)
  return(newmat)
}

#iterate over generations
generations <- 200
iterate.fitness <- function(generations,gen,a,b,c) {
  fitness <- list()
  rplct <- list()
  pools <- list()
  nextgen <- list()
  for(i in 1:generations){
    fitness[[i]] <- w.calc(gen)
    rplct[[i]] <- replicate(a, b, c, fitness[[i]])
    pools[[i]] <- pool(rplct[[i]])
    nextgen[[i]] <- reassort(pools[[i]])
    gen <- nextgen[[i]]
    a <- as.data.frame(cbind("x" = gen[1,], "n0" = initial.pop))
    b <- as.data.frame(cbind("x" = gen[2,], "n0" = initial.pop))
    c <- as.data.frame(cbind("x" = gen[3,], "n0" = initial.pop))
  }
  fit.means <- unlist(lapply(fitness, mean))
  return(fit.means)
}

runs <- list()
for(i in 1:100){
  set.seed(i)
  runs[[i]] <- iterate.fitness(generations,gen.mat, a.group, b.group,c.group)
}
run.tab <- bind_cols(runs)
colnames(run.tab) <- paste0("run",c(1:length(runs)))

run.tab <- run.tab %>%
  mutate("generation" = c(1:generations)) %>%
  relocate(generation, .before = run1)

run.gather <- gather(run.tab,key = run,value = "fitness", -generation)

ggplot(data = run.gather, aes(x = generation, y = fitness, color = run))+
  geom_line()+
  xlab("Generation")+
  ylab("Population fitness")+
  theme_bw()+
  scale_color_manual(values = rep("black", 100))+
  theme(legend.position = "none")

run.max <- unlist(lapply(runs, max))
length(which(run.max > 1)) #68 > 1
length(which(run.max <= 1)) #32 <= 1
mean(run.max) #1.015814
run.max <- as.data.frame(run.max)
ggplot(data = run.max, aes(x = run.max))+
  geom_density(bins = 30)

#now try without reassortment
generations <- 200
iterate.nomix <- function(generations,gen,a,b,c) {
  fitness <- list()
  rplct <- list()
  pools <- list()
  nextgen <- list()
  for(i in 1:generations){
    fitness[[i]] <- w.calc(gen)
    rplct[[i]] <- replicate(a, b, c, fitness[[i]])
    pools[[i]] <- pool(rplct[[i]])
    nextgen[[i]] <- reassort(pools[[i]])
    gen <- nextgen[[i]]
    a <- as.data.frame(cbind("x" = gen[1,], "n0" = initial.pop))
    b <- as.data.frame(cbind("x" = gen[2,], "n0" = initial.pop))
    c <- as.data.frame(cbind("x" = gen[3,], "n0" = initial.pop))
  }
  fit.means <- unlist(lapply(fitness, mean))
  return(fit.means)
}






#there's more to pull out from this function than just pop fitness
#this works as far as i can tell
fitness <- list()
rplct <- list()
pools <- list()
nextgen <- list()

for(i in 1:generations){
  fitness[[i]] <- w.calc(gen.mat)
  rplct[[i]] <- replicate(a.group, b.group, c.group, fitness[[i]])
  pools[[i]] <- pool(rplct[[i]])
  nextgen[[i]] <- reassort(pools[[i]])
  gen.mat <- nextgen[[i]]
  a.group <- as.data.frame(cbind("x" = gen.mat[1,], "n0" = initial.pop))
  b.group <- as.data.frame(cbind("x" = gen.mat[2,], "n0" = initial.pop))
  c.group <- as.data.frame(cbind("x" = gen.mat[3,], "n0" = initial.pop))
}
  

fit.means <- unlist(lapply(fitness, mean))
plot(1:generations,fit.means, type = "l", xlab = "generation", ylab = "population fitness")
#omg!!!!

