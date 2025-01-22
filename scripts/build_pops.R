library(sn)
library(tidyverse)

#generate fitness distribution
dms.dist <- read_csv("dms_dist.csv")
colnames(dms.dist) <- "x"
dms.dist[which(dms.dist$x <= 0),] <- 0.00001
dms.dist <- na.omit(dms.dist)

ggplot(data = dms.dist, aes(x = x))+
  geom_histogram(bins = 30, fill = "grey", color = "black")+
  scale_x_continuous(limits = c(0, 1.5), oob = scales::oob_keep)+
  ylab("Count")+
  ggtitle("No immune pressure")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 22))
ggsave("runs/run_010725/figs/nopressure_hist.png")

#generaate fitness distributions under immune pressure
viable.dist <- dms.dist[which(dms.dist$x > 0),]
set.seed(1)
low.dist <- as.data.frame(c(rep(0,950),sample(viable.dist$x,50))) #5% resistance
colnames(low.dist) <- "x"
set.seed(1)
mid.dist <- as.data.frame(c(rep(0,900),sample(viable.dist$x,100))) #10% resistance
colnames(mid.dist) <- "x"
set.seed(1)
high.dist <- as.data.frame(c(rep(0,800),sample(viable.dist$x,200)))
colnames(high.dist) <- "x"

#low resistance (5%)
ggplot(data = low.dist, aes(x = x))+
  geom_histogram(bins = 30, fill = "grey", color = "black")+
  scale_x_continuous(limits = c(0, 1.5), oob = scales::oob_keep)+
  ylab("Count")+
  ggtitle("Immune pressure, 5% resistance")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 22))
ggsave("runs/run_010725/figs/lowresistance_hist.png")

#mid resistance (10%)
ggplot(data = mid.dist, aes(x = x))+
  geom_histogram(bins = 30, fill = "grey", color = "black")+
  scale_x_continuous(limits = c(0, 1.5), oob = scales::oob_keep)+
  ylab("Count")+
  ggtitle("Immune pressure, 10% resistance")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 22))
ggsave("runs/run_010725/figs/midresistance_hist.png")

#high resistance (20%)
ggplot(data = high.dist, aes(x = x))+
  geom_histogram(bins = 30, fill = "grey", color = "black")+
  scale_x_continuous(limits = c(0, 1.5), oob = scales::oob_keep)+
  ylab("Count")+
  ggtitle("Immune pressure, 20% resistance")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 22))
ggsave("runs/run_010725/figs/highresistance_hist.png")


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
ggsave("runs/run_010725/figs/r_punctuated.png")

ggplot(data = r.diff, aes(x = r))+
  geom_histogram(bins = 50, fill = "grey", color = "black")+
  ylim(c(0,8))+
  ylab("Count")+
  ggtitle("Diversity distribution: diffuse")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        plot.title = element_text(size = 22))
ggsave("runs/run_010725/figs/r_diffuse.png")

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
build.ab.punc <- function(fit,fit.imm,off,r){
  gen.mat <- list()
  for(i in 1:2){
    gen.mat[[i]] <- matrix(data = NA, nrow = 3, ncol = 2)
    rownames(gen.mat[[i]]) <- c("A","B","C")
    colnames(gen.mat[[i]]) <- c("x","r")
    gen.mat[[i]][1,1] <- sample(fit.imm$x, 1)
    gen.mat[[i]][2,1] <- sample(fit.imm$x, 1)
    gen.mat[[i]][3,1] <- sample(fit$x, 1)
    gen.mat[[i]][1,2] <- sample(r$r,1)
    gen.mat[[i]][2,2] <- gen.mat[[i]][1,2] + (sample(c(1,-1), 1, replace = T)) * sample(off,1)
    gen.mat[[i]][3,2] <- gen.mat[[i]][1,2] + (sample(c(1,-1), 1, replace = T)) * sample(off,1)
  }
  names(gen.mat) <- c("i1","i2")
  return(gen.mat)
}
build.c.punc <- function(fit,fit.imm,off,r){
  gen.mat <- list()
  for(i in 1:2){
    gen.mat[[i]] <- matrix(data = NA, nrow = 3, ncol = 2)
    rownames(gen.mat[[i]]) <- c("A","B","C")
    colnames(gen.mat[[i]]) <- c("x","r")
    gen.mat[[i]][1,1] <- sample(fit$x, 1)
    gen.mat[[i]][2,1] <- sample(fit$x, 1)
    gen.mat[[i]][3,1] <- sample(fit.imm$x, 1)
    gen.mat[[i]][1,2] <- sample(r$r,1)
    gen.mat[[i]][2,2] <- gen.mat[[i]][1,2] + (sample(c(1,-1), 1, replace = T)) * sample(off,1)
    gen.mat[[i]][3,2] <- gen.mat[[i]][1,2] + (sample(c(1,-1), 1, replace = T)) * sample(off,1)
  }
  names(gen.mat) <- c("i1","i2")
  return(gen.mat)
}

set.seed(1)
gen.list.abc.punc <- replicate(100, build.abc.punc(dms.dist,off.dist,r.punc), 
                          simplify = F)
set.seed(1)
gen.list.ab.punc <- replicate(100, build.ab.punc(dms.dist,mid.dist,off.dist,
                                                r.punc), simplify = F)
set.seed(1)
gen.list.c.punc <- replicate(100, build.c.punc(dms.dist,mid.dist,off.dist,
                                              r.punc), simplify = F)
gen.lists.punc <- list("abc" = gen.list.abc.punc, "ab" = gen.list.ab.punc,
                       "c" = gen.list.c.punc)
save(gen.lists.punc, file = "runs/run_010725/rdata/gen_lists_punctuated.RData")

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
build.ab.diff <- function(fit,fit.imm,r){
  gen.mat <- list()
  for(i in 1:2){
    gen.mat[[i]] <- matrix(data = NA, nrow = 3, ncol = 2)
    rownames(gen.mat[[i]]) <- c("A","B","C")
    colnames(gen.mat[[i]]) <- c("x","r")
    gen.mat[[i]][1,1] <- sample(fit.imm$x, 1)
    gen.mat[[i]][2,1] <- sample(fit.imm$x, 1)
    gen.mat[[i]][3,1] <- sample(fit$x, 1)
    gen.mat[[i]][,2] <- sample(r$r,3)
  }
  names(gen.mat) <- c("i1","i2")
  return(gen.mat)
}
build.c.diff <- function(fit,fit.imm,r){
  gen.mat <- list()
  for(i in 1:2){
    gen.mat[[i]] <- matrix(data = NA, nrow = 3, ncol = 2)
    rownames(gen.mat[[i]]) <- c("A","B","C")
    colnames(gen.mat[[i]]) <- c("x","r")
    gen.mat[[i]][1,1] <- sample(fit$x, 1)
    gen.mat[[i]][2,1] <- sample(fit$x, 1)
    gen.mat[[i]][3,1] <- sample(fit.imm$x, 1)
    gen.mat[[i]][,2] <- sample(r$r,3)
  }
  names(gen.mat) <- c("i1","i2")
  return(gen.mat)
}

set.seed(1)
gen.list.abc.diff <- replicate(100, build.abc.diff(dms.dist,r.diff), 
                               simplify = F)
set.seed(1)
gen.list.ab.diff <- replicate(100, build.ab.diff(dms.dist,mid.dist,r.diff), 
                              simplify = F)
set.seed(1)
gen.list.c.diff <- replicate(100, build.c.diff(dms.dist,mid.dist,r.diff),
                             simplify = F)
gen.lists.diff <- list("abc" = gen.list.abc.diff, "ab" = gen.list.ab.diff,
                       "c" = gen.list.c.diff)

save(gen.lists.diff, file = "runs/run_010725/rdata/gen_lists_diffuse.RData")

#build punc x diff coinfection pairs
gen.lists.pd <- list("abc" = NA, "ab" = NA, "c" = NA)
for(i in seq_along(gen.lists.punc)){
  gen.lists.pd[[i]] <- map2(gen.lists.punc[[i]], gen.lists.diff[[i]], function(x,y) list("i1"= x$i1,"i2" = y$i2))
}
save(gen.lists.pd, file = "runs/run_010725/rdata/gen_lists_pd.RData")

#generate list of possible genotype combinations
combo <- as.data.frame(matrix(data = c("A1","A2","B1","B2","C1","C2"), ncol = 3, nrow = 2))
colnames(combo) <- c("A","B","C")
combo <- expand(combo,A,B,C)
combo <- split(combo,1:nrow(combo))
combo <- unlist(lapply(combo, paste0, collapse = ""))
save(combo,file = "runs/run_010725/rdata/combo.RData")


#different resistance levels
#low resistance
set.seed(1)
gen.list.abc.low <- replicate(100, build.abc.punc(dms.dist,off.dist,r.punc), 
                               simplify = F)
set.seed(1)
gen.list.ab.low <- replicate(100, build.ab.punc(dms.dist,low.dist,off.dist,
                                                 r.punc), simplify = F)
set.seed(1)
gen.list.c.low <- replicate(100, build.c.punc(dms.dist,low.dist,off.dist,
                                               r.punc), simplify = F)
gen.lists.low <- list("abc" = gen.list.abc.low, "ab" = gen.list.ab.low,
                       "c" = gen.list.c.low)
save(gen.lists.low, file = "runs/run_010725/rdata/gen_lists_lowres.RData")

#high resistance
set.seed(1)
gen.list.abc.high <- replicate(100, build.abc.punc(dms.dist,off.dist,r.punc), 
                              simplify = F)
set.seed(1)
gen.list.ab.high <- replicate(100, build.ab.punc(dms.dist,high.dist,off.dist,
                                                r.punc), simplify = F)
set.seed(1)
gen.list.c.high <- replicate(100, build.c.punc(dms.dist,high.dist,off.dist,
                                              r.punc), simplify = F)
gen.lists.high <- list("abc" = gen.list.abc.high, "ab" = gen.list.ab.high,
                      "c" = gen.list.c.high)
save(gen.lists.high, file = "runs/run_010725/rdata/gen_lists_highres.RData")

