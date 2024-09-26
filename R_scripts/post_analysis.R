library(tidyverse)

#load replication statistics
load("runs/run_091024/rep_stats_init1.RData")
load("runs/run_091024/rep_stats_init5.RData")

#final viable genome counts
count.genomes.mu0 <- function(repx){
  counts <- list()
  viable <- list()
  for(i in seq_along(repx)){
    counts[[i]] <- repx[[i]][dim(repx[[i]])[1],]
    viable[[i]] <- (min(counts[[i]]$A1, counts[[i]]$B1, counts[[i]]$C1))+
      (min(counts[[i]]$A2, counts[[i]]$B2, counts[[i]]$C2))
  }
  return(viable)
}

viable.mu0.init1 <- map(rep.stats.init1$mu0, count.genomes.mu0)
viable.mu0.init5 <- map(rep.stats.init5$mu0, count.genomes.mu0)

count.genomes.muplus <- function(repx){
  counts <- list()
  viable <- list()
  for(i in seq_along(repx)){
    counts[[i]] <- repx[[i]][dim(repx[[i]])[1],]
    viable[[i]] <- min((counts[[i]]$A1 + counts[[i]]$A2),
                       (counts[[i]]$B1 + counts[[i]]$B2),
                       (counts[[i]]$C1 + counts[[i]]$C2))
  }
  return(viable)
}

viable.mu01.init1 <- map(rep.stats.init1$mu01, count.genomes.muplus)
viable.mu01.init5 <- map(rep.stats.init5$mu01, count.genomes.muplus)

viable.mu05.init1 <- map(rep.stats.init1$mu05, count.genomes.muplus)
viable.mu05.init5 <- map(rep.stats.init5$mu05, count.genomes.muplus)

viable.mu09.init1 <- map(rep.stats.init1$mu09, count.genomes.muplus)
viable.mu09.init5 <- map(rep.stats.init5$mu09, count.genomes.muplus)

viable.mu1.init1 <- map(rep.stats.init1$mu1, count.genomes.muplus)
viable.mu1.init5 <- map(rep.stats.init5$mu1, count.genomes.muplus)

viable.init1.list <- list("mu0" = viable.mu0.init1, "mu01" = viable.mu01.init1, 
                          "mu05" = viable.mu05.init1,"mu09" = viable.mu09.init1,
                          "mu1" = viable.mu1.init1)
viable.init5.list <- list("mu0" = viable.mu0.init5, "mu01" = viable.mu01.init5, 
                          "mu05" = viable.mu05.init5,"mu09" = viable.mu09.init5,
                          "mu1" = viable.mu1.init5)

viable.lists <- list("init1" = viable.init1.list, "init5" = viable.init5.list)

build.viable.tables <- function(viables){
  viable.abc <- map(viables, function(x) x$abc)
  viable.abc <- bind_rows(map(viable.abc, unlist))
  viable.abc.gather <- gather(viable.abc, key = "mu", value = "count")
  
  viable.ab <- map(viables, function(x) x$ab)
  viable.ab <- bind_rows(map(viable.ab, unlist))
  viable.ab.gather <- gather(viable.ab, key = "mu", value = "count")
  
  viable.c <- map(viables, function(x) x$c)
  viable.c <- bind_rows(map(viable.c, unlist))
  viable.c.gather <- gather(viable.c, key = "mu", value = "count")
  
  viable.tables <- list("abc" = viable.abc.gather, "ab" = viable.ab.gather,
                        "c" = viable.c.gather)
  return(viable.tables)
  
}

viable.tables <- map(viable.lists, build.viable.tables)

viable.plots.init1 <- list()
for(i in seq_along(viable.tables$init1)){
  titles <- c("A = B = C", "AB > C", "C > AB")
  viable.plots.init1[[i]] <- 
    ggplot(data = viable.tables$init1[[i]], aes(x = count, fill = mu))+
    geom_density(alpha = 0.6)+
    xlab("Viable genomes")+
    ylab("Density")+
    ggtitle(titles[[i]])+
    scale_fill_manual(values = mu.cols, labels = c(0, 0.1, 0.5, 0.9, 1))+
    theme_bw()+
    theme(axis.text = element_text(size = 19),
          axis.title = element_text(size = 22),
          plot.title = element_text(size = 22),
          legend.text = element_text(size = 19),
          legend.title = element_text(size = 22))
  filenames <- paste0("runs/run_091024/figs/",
                      c("viable_genomes_abc_init1.png", "viable_genomes_ab_init1.png",
                        "viable_genomes_c_init1.png"))
  ggsave(filenames[[i]])
}

viable.plots.init5 <- list()
for(i in seq_along(viable.tables$init5)){
  titles <- c("A = B = C", "AB > C", "C > AB")
  viable.plots.init5[[i]] <- 
    ggplot(data = viable.tables$init5[[i]], aes(x = count, fill = mu))+
    geom_density(alpha = 0.6)+
    xlab("Viable genomes")+
    ylab("Density")+
    ggtitle(titles[[i]])+
    scale_fill_manual(values = mu.cols, labels = c(0, 0.1, 0.5, 0.9, 1))+
    theme_bw()+
    theme(axis.text = element_text(size = 19),
          axis.title = element_text(size = 22),
          plot.title = element_text(size = 22),
          legend.text = element_text(size = 19),
          legend.title = element_text(size = 22))
  filenames <- paste0("runs/run_091024/figs/",
                      c("viable_genomes_abc_init5.png", "viable_genomes_ab_init5.png",
                      "viable_genomes_c_init5.png"))
  ggsave(filenames[[i]])
}

#packaging fitness stats
load("runs/run_091024/pack_stats_init1.RData")
load("runs/run_091024/pack_stats_init5.RData")

tab.final.genotypes <- function(packx){
  final.genotypes <- list("abc" = NA, "ab" = NA, "c" = NA)
  for(i in seq_along(packx)){
    final.genotypes[[i]] <- map_dfr(packx[[i]], function(x) x[dim(x)[1],c(8:15)])
  }
  return(final.genotypes)
}

genotypes.init1 <- map(pack.stats.init1, tab.final.genotypes)
genotypes.init5 <- map(pack.stats.init5, tab.final.genotypes)

#calculate fitness for each genotype
load("runs/run_091024/gen_lists.RData")
load("runs/run_091024/combo.RData")
calc.fitness <- function(gen){
  fitness <- as.data.frame(matrix(ncol = 8, nrow = 1))
  colnames(fitness) <- combo
  fitness[1] <- mean(c(gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i1[2,2])),
                       gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i1[1,2])), 
                       gen$i1[3,1])) #A1B1C1
  fitness[2] <- mean(c(gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i1[2,2])),
                       gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i1[1,2])), 
                       gen$i2[3,1])) #A1B1C2
  fitness[3] <- mean(c(gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i2[2,2])),
                       gen$i2[2,1] * (1 - abs(gen$i2[2,2] - gen$i1[1,2])), 
                       gen$i1[3,1])) #A1B2C1
  fitness[4] <- mean(c(gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i2[2,2])),
                       gen$i2[2,1] * (1 - abs(gen$i2[2,2] - gen$i1[1,2])), 
                       gen$i2[3,1])) #A1B2C2
  fitness[5] <- mean(c(gen$i2[1,1] * (1 - abs(gen$i2[1,2] - gen$i1[2,2])),
                       gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i2[1,2])), 
                       gen$i1[3,1])) #A2B1C1
  fitness[6] <- mean(c(gen$i2[1,1] * (1 - abs(gen$i2[1,2] - gen$i1[2,2])),
                       gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i2[1,2])), 
                       gen$i2[3,1])) #A2B1C2
  fitness[7] <- mean(c(gen$i2[1,1] * (1 - abs(gen$i2[1,2] - gen$i2[2,2])),
                       gen$i2[2,1] * (1 - abs(gen$i2[2,2] - gen$i2[1,2])), 
                       gen$i1[3,1])) #A2B2C1
  fitness[8] <- mean(c(gen$i2[1,1] * (1 - abs(gen$i2[1,2] - gen$i2[2,2])),
                       gen$i2[2,1] * (1 - abs(gen$i2[2,2] - gen$i2[1,2])), 
                       gen$i2[3,1])) #A2B2C2
  return(fitness)
}

genotype.fitness <- list("abc" = NA, "ab" = NA, "c" = NA)
for(i in seq_along(gen.lists)){
  genotype.fitness[[i]] <- map_dfr(gen.lists[[i]], calc.fitness)
}

counts <- genotypes.init1$mu0$abc
fit <- genotype.fitness$abc

calc.weighted.fitness <- function(fit, counts){
  weighted.fit <- list()
  for(i in seq_along(rownames(fit))){
    weighted.fit[[i]] <- weighted.mean(fit[i,], counts[i,])
  }
  weighted.fit <- unlist(weighted.fit)
  return(weighted.fit)
}

weighted.fitness.init1 <- list("mu0" = NA, "mu01" = NA, "mu05" = NA, "mu09" = NA, "mu1" = NA)
for(i in seq_along(genotypes.init1)){
  weighted.fitness.init1[[i]] <- map2(genotypes.init1[[i]], genotype.fitness, calc.weighted.fitness)
}

weighted.fitness.init5 <- list("mu0" = NA, "mu01" = NA, "mu05" = NA, "mu09" = NA, "mu1" = NA)
for(i in seq_along(genotypes.init5)){
  weighted.fitness.init5[[i]] <- map2(genotypes.init5[[i]], genotype.fitness, calc.weighted.fitness)
}

