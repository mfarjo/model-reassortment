library(tidyverse)

#generate random normally distributed and lognormally distributed fitness values
set.seed(1)
fit.dist.norm <- as.data.frame(rnorm(100, mean = 1 , sd = 0.1))
colnames(fit.dist.norm) <- "x"
ggplot(data = fit.dist.norm, aes(x = x))+
  geom_histogram(bins = 30, fill = "grey", color = "black")+
  #xlim(c(0,1.5))+
  ylim(c(0,15))+
  ylab("Count")+
  ggtitle("No immune pressure")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 22))
ggsave("runs/run_120224/figs/nopressure_hist.png")

set.seed(1)
fit.dist.total <- as.data.frame(rgamma(90,shape = 0.5, scale = 0.3))
colnames(fit.dist.total) <- "x"
ggplot(data = fit.dist.total, aes(x = x))+
  geom_histogram(bins = 30, fill = "grey", color = "black")+
  #xlim(c(0,1.2))+
  #ylim(c(0,15))+
  ylab("Count")+
  ggtitle("Immune pressure")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 22))
ggsave("runs/run_120224/figs/pressure_hist.png")


#generate random normally distributed small fitness offsets
set.seed(1)
off.dist <- rnorm(100, mean = 0.02, sd = 0.001)

#generate relatedness spectra
#punctuated diversity
set.seed(1)
r.punc <- as.data.frame(c(runif(25, 0, 0.1), runif(25, 0.2,0.3), 
                     runif(25, 0.7,0.8), runif(25, 0.9,1)))
colnames(r.punc) <- "r"

#diffuse diversity
set.seed(1)
r.diff <- as.data.frame(runif(100,0,1))
colnames(r.diff) <- "r"

ggplot(data = r.punc, aes(x = r))+
  geom_histogram(bins = 50, fill = "grey", color = "black")+
  ylim(c(0,8))+
  ylab("Count")+
  ggtitle("Diversity distribution: punctuated")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 22))
ggsave("runs/run_120224/figs/r_punctuated.png")

ggplot(data = r.diff, aes(x = r))+
  geom_histogram(bins = 50, fill = "grey", color = "black")+
  ylim(c(0,8))+
  ylab("Count")+
  ggtitle("Diversity distribution: diffuse")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 22))
ggsave("runs/run_120224/figs/r_diffuse.png")

#build genotypes
build.abc.punc <- function(fit,off,r){
  gen.mat <- list()
  for(i in 1:2){
    gen.mat[[i]] <- matrix(data = NA, nrow = 3, ncol = 2)
    rownames(gen.mat[[i]]) <- c("A","B","C")
    colnames(gen.mat[[i]]) <- c("x","r")
    gen.mat[[i]][,1] <- sample(fit$x, 3)
    gen.mat[[i]][1,2] <- sample(r$r,1)
    gen.mat[[i]][2,2] <- gen.mat[[i]][1,2] + (sample(c(1,-1), 1, replace = T)) * sample(off,1)
    gen.mat[[i]][3,2] <- gen.mat[[i]][1,2] + (sample(c(1,-1), 1, replace = T)) * sample(off,1)
  }
  names(gen.mat) <- c("i1","i2")
  return(gen.mat)
}
build.ab.punc <- function(fit,fit.plus,off,r){
  gen.mat <- list()
  for(i in 1:2){
    gen.mat[[i]] <- matrix(data = NA, nrow = 3, ncol = 2)
    rownames(gen.mat[[i]]) <- c("A","B","C")
    colnames(gen.mat[[i]]) <- c("x","r")
    gen.mat[[i]][1,1] <- sample(fit.plus$x, 1)
    gen.mat[[i]][2,1] <- sample(fit.plus$x, 1)
    gen.mat[[i]][3,1] <- sample(fit$x, 1)
    gen.mat[[i]][1,2] <- sample(r$r,1)
    gen.mat[[i]][2,2] <- gen.mat[[i]][1,2] + (sample(c(1,-1), 1, replace = T)) * sample(off,1)
    gen.mat[[i]][3,2] <- gen.mat[[i]][1,2] + (sample(c(1,-1), 1, replace = T)) * sample(off,1)
  }
  names(gen.mat) <- c("i1","i2")
  return(gen.mat)
}
build.c.punc <- function(fit,fit.plus,off,r){
  gen.mat <- list()
  for(i in 1:2){
    gen.mat[[i]] <- matrix(data = NA, nrow = 3, ncol = 2)
    rownames(gen.mat[[i]]) <- c("A","B","C")
    colnames(gen.mat[[i]]) <- c("x","r")
    gen.mat[[i]][1,1] <- sample(fit$x, 1)
    gen.mat[[i]][2,1] <- sample(fit$x, 1)
    gen.mat[[i]][3,1] <- sample(fit.plus$x, 1)
    gen.mat[[i]][1,2] <- sample(r$r,1)
    gen.mat[[i]][2,2] <- gen.mat[[i]][1,2] + (sample(c(1,-1), 1, replace = T)) * sample(off,1)
    gen.mat[[i]][3,2] <- gen.mat[[i]][1,2] + (sample(c(1,-1), 1, replace = T)) * sample(off,1)
  }
  names(gen.mat) <- c("i1","i2")
  return(gen.mat)
}

set.seed(1)
gen.list.abc.punc <- replicate(100, build.abc.punc(fit.dist.norm,off.dist,r.punc), 
                          simplify = F)
set.seed(1)
gen.list.ab.punc <- replicate(100, build.ab.punc(fit.dist.norm,fit.dist.total,off.dist,
                                                r.punc), simplify = F)
set.seed(1)
gen.list.c.punc <- replicate(100, build.c.punc(fit.dist.norm,fit.dist.total,off.dist,
                                              r.punc), simplify = F)
gen.lists.punc <- list("abc" = gen.list.abc.punc, "ab" = gen.list.ab.punc,
                       "c" = gen.list.c.punc)
save(gen.lists.punc, file = "runs/run_120224/rdata/gen_lists_punctuated.RData")

build.abc.diff <- function(fit,r){
  gen.mat <- list()
  for(i in 1:2){
    gen.mat[[i]] <- matrix(data = NA, nrow = 3, ncol = 2)
    rownames(gen.mat[[i]]) <- c("A","B","C")
    colnames(gen.mat[[i]]) <- c("x","r")
    gen.mat[[i]][,1] <- sample(fit$x, 3)
    gen.mat[[i]][,2] <- sample(r$r,3)
  }
  names(gen.mat) <- c("i1","i2")
  return(gen.mat)
}
build.ab.diff <- function(fit,fit.plus,r){
  gen.mat <- list()
  for(i in 1:2){
    gen.mat[[i]] <- matrix(data = NA, nrow = 3, ncol = 2)
    rownames(gen.mat[[i]]) <- c("A","B","C")
    colnames(gen.mat[[i]]) <- c("x","r")
    gen.mat[[i]][1,1] <- sample(fit.plus$x, 1)
    gen.mat[[i]][2,1] <- sample(fit.plus$x, 1)
    gen.mat[[i]][3,1] <- sample(fit$x, 1)
    gen.mat[[i]][,2] <- sample(r$r,3)
  }
  names(gen.mat) <- c("i1","i2")
  return(gen.mat)
}
build.c.diff <- function(fit,fit.plus,r){
  gen.mat <- list()
  for(i in 1:2){
    gen.mat[[i]] <- matrix(data = NA, nrow = 3, ncol = 2)
    rownames(gen.mat[[i]]) <- c("A","B","C")
    colnames(gen.mat[[i]]) <- c("x","r")
    gen.mat[[i]][1,1] <- sample(fit$x, 1)
    gen.mat[[i]][2,1] <- sample(fit$x, 1)
    gen.mat[[i]][3,1] <- sample(fit.plus$x, 1)
    gen.mat[[i]][,2] <- sample(r$r,3)
  }
  names(gen.mat) <- c("i1","i2")
  return(gen.mat)
}

set.seed(1)
gen.list.abc.diff <- replicate(100, build.abc.diff(fit.dist.norm,r.diff), 
                               simplify = F)
set.seed(1)
gen.list.ab.diff <- replicate(100, build.ab.diff(fit.dist.norm,fit.dist.total,r.diff), 
                              simplify = F)
set.seed(1)
gen.list.c.diff <- replicate(100, build.c.diff(fit.dist.norm,fit.dist.total,r.diff),
                             simplify = F)
gen.lists.diff <- list("abc" = gen.list.abc.diff, "ab" = gen.list.ab.diff,
                       "c" = gen.list.c.diff)

save(gen.lists.diff, file = "runs/run_120224/rdata/gen_lists_diffuse.RData")

#generate list of possible genotype combinations
combo <- as.data.frame(matrix(data = c("A1","A2","B1","B2","C1","C2"), ncol = 3, nrow = 2))
colnames(combo) <- c("A","B","C")
combo <- expand(combo,A,B,C)
combo <- split(combo,1:nrow(combo))
combo <- unlist(lapply(combo, paste0, collapse = ""))
save(combo,file = "runs/run_120224/rdata/combo.RData")


