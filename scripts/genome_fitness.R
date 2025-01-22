library(tidyverse)
library(ggpubr)
library(dunn.test)

#packaging fitness stats
load("runs/run_010725/rdata/final_gen_punc.RData")
load("runs/run_010725/rdata/final_gen_diff.RData")
load("runs/run_010725/rdata/final_gen_diff_noepi.RData")
load("runs/run_010725/rdata/final_gen_pd_epi.RData")
load("runs/run_010725/rdata/final_gen_pd_noepi.RData")
load("runs/run_010725/rdata/final_gen_pd_midepi.RData")
load("runs/run_010725/rdata/final_gen_pd_puncepi.RData")
counts.list <- list("punc" = final.gen.punc, "diff" = final.gen.diff,
                    "diff_noepi" = final.gen.diff.noepi, "pd_epi" = final.gen.pd.epi,
                    "pd_noepi" = final.gen.pd.noepi, "pd_midepi" = final.gen.pd.midepi,
                    "pd_puncepi" = final.gen.pd.puncepi)
save(counts.list, file = "runs/run_010725/rdata/counts_list.RData")


#calculate fitness for each genotype
load("runs/run_010725/rdata/gen_lists_punctuated.RData")
load("runs/run_010725/rdata/gen_lists_diffuse.RData")
load("runs/run_010725/rdata/gen_lists_pd.RData")
load("runs/run_010725/rdata/combo.RData")

calc.fitness <- function(gen){
  fitness <- as.data.frame(matrix(ncol = 8, nrow = 1))
  colnames(fitness) <- combo
  fitness[1] <- mean(c((gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i1[2,2]))),
                       (gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i1[1,2]))), 
                       (gen$i1[3,1]))) #A1B1C1
  fitness[2] <- mean(c((gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i1[2,2]))),
                       (gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i1[1,2]))), 
                       (gen$i2[3,1]))) #A1B1C2
  fitness[3] <- mean(c((gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i2[2,2]))),
                       (gen$i2[2,1] * (1 - abs(gen$i2[2,2] - gen$i1[1,2]))), 
                       (gen$i1[3,1]))) #A1B2C1
  fitness[4] <- mean(c((gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i2[2,2]))),
                       (gen$i2[2,1] * (1 - abs(gen$i2[2,2] - gen$i1[1,2]))), 
                       (gen$i2[3,1]))) #A1B2C2
  fitness[5] <- mean(c((gen$i2[1,1] * (1 - abs(gen$i2[1,2] - gen$i1[2,2]))),
                       (gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i2[1,2]))), 
                       (gen$i1[3,1]))) #A2B1C1
  fitness[6] <- mean(c((gen$i2[1,1] * (1 - abs(gen$i2[1,2] - gen$i1[2,2]))),
                       (gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i2[1,2]))), 
                       (gen$i2[3,1]))) #A2B1C2
  fitness[7] <- mean(c((gen$i2[1,1] * (1 - abs(gen$i2[1,2] - gen$i2[2,2]))),
                       (gen$i2[2,1] * (1 - abs(gen$i2[2,2] - gen$i2[1,2]))), 
                       (gen$i1[3,1]))) #A2B2C1
  fitness[8] <- mean(c((gen$i2[1,1] * (1 - abs(gen$i2[1,2] - gen$i2[2,2]))),
                       (gen$i2[2,1] * (1 - abs(gen$i2[2,2] - gen$i2[1,2]))), 
                       (gen$i2[3,1]))) #A2B2C2
  return(fitness)
}

calc.fitness.diff.noepi <- function(gen){
  fitness <- as.data.frame(matrix(ncol = 8, nrow = 1))
  colnames(fitness) <- combo
  fitness[1] <- mean(c(gen$i1[1,1],
                       gen$i1[2,1], 
                       gen$i1[3,1])) #A1B1C1
  fitness[2] <- mean(c(gen$i1[1,1],
                       gen$i1[2,1], 
                       gen$i2[3,1])) #A1B1C2
  fitness[3] <- mean(c(gen$i1[1,1],
                       gen$i2[2,1], 
                       gen$i1[3,1])) #A1B2C1
  fitness[4] <- mean(c(gen$i1[1,1],
                       gen$i2[2,1], 
                       gen$i2[3,1])) #A1B2C2
  fitness[5] <- mean(c(gen$i2[1,1],
                       gen$i1[2,1], 
                       gen$i1[3,1])) #A2B1C1
  fitness[6] <- mean(c(gen$i2[1,1],
                       gen$i1[2,1], 
                       gen$i2[3,1])) #A2B1C2
  fitness[7] <- mean(c(gen$i2[1,1],
                       gen$i2[2,1], 
                       gen$i1[3,1])) #A2B2C1
  fitness[8] <- mean(c(gen$i2[1,1],
                       gen$i2[2,1], 
                       gen$i2[3,1])) #A2B2C2
  return(fitness)
}

calc.fitness.pd.noepi <- function(gen){
  fitness <- as.data.frame(matrix(ncol = 8, nrow = 1))
  colnames(fitness) <- combo
  fitness[1] <- mean(c((gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i1[2,2]))),
                       (gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i1[1,2]))), 
                       (gen$i1[3,1]))) #A1B1C1
  fitness[2] <- mean(c((gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i1[2,2]))),
                       (gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i1[1,2]))), 
                       (gen$i2[3,1]))) #A1B1C2
  fitness[3] <- mean(c(gen$i1[1,1],
                       gen$i2[2,1], 
                       gen$i1[3,1])) #A1B2C1
  fitness[4] <- mean(c(gen$i1[1,1],
                       gen$i2[2,1], 
                       gen$i2[3,1])) #A1B2C2
  fitness[5] <- mean(c(gen$i2[1,1],
                       gen$i1[2,1], 
                       gen$i1[3,1])) #A2B1C1
  fitness[6] <- mean(c(gen$i2[1,1],
                       gen$i1[2,1], 
                       gen$i2[3,1])) #A2B1C2
  fitness[7] <- mean(c(gen$i2[1,1],
                       gen$i2[2,1], 
                       gen$i1[3,1])) #A2B2C1
  fitness[8] <- mean(c(gen$i2[1,1],
                       gen$i2[2,1], 
                       gen$i2[3,1])) #A2B2C2
  return(fitness)
}

calc.fitness.pd.midepi <- function(gen){
  fitness <- as.data.frame(matrix(ncol = 8, nrow = 1))
  colnames(fitness) <- combo
  fitness[1] <- mean(c((gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i1[2,2]))),
                       (gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i1[1,2]))), 
                       (gen$i1[3,1]))) #A1B1C1
  fitness[2] <- mean(c((gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i1[2,2]))),
                       (gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i1[1,2]))), 
                       (gen$i2[3,1]))) #A1B1C2
  fitness[3] <- mean(c((gen$i1[1,1] * (1 - ((abs(gen$i1[1,2] - gen$i2[2,2]))/2))),
                       (gen$i2[2,1] * (1 - ((abs(gen$i2[2,2] - gen$i1[1,2]))/2))), 
                       (gen$i1[3,1]))) #A1B2C1
  fitness[4] <- mean(c((gen$i1[1,1] * (1 - ((abs(gen$i1[1,2] - gen$i2[2,2]))/2))),
                       (gen$i2[2,1] * (1 - ((abs(gen$i2[2,2] - gen$i1[1,2]))/2))), 
                       (gen$i2[3,1]))) #A1B2C2
  fitness[5] <- mean(c((gen$i2[1,1] * (1 - ((abs(gen$i2[1,2] - gen$i1[2,2]))/2))),
                       (gen$i1[2,1] * (1 - ((abs(gen$i1[2,2] - gen$i2[1,2]))/2))), 
                       (gen$i1[3,1]))) #A2B1C1
  fitness[6] <- mean(c((gen$i2[1,1] * (1 - ((abs(gen$i2[1,2] - gen$i1[2,2]))/2))),
                       (gen$i1[2,1] * (1 - ((abs(gen$i1[2,2] - gen$i2[1,2]))/2))), 
                       (gen$i2[3,1]))) #A2B1C2
  fitness[7] <- mean(c(gen$i2[1,1],
                       gen$i2[2,1], 
                       gen$i1[3,1])) #A2B2C1
  fitness[8] <- mean(c(gen$i2[1,1],
                       gen$i2[2,1], 
                       gen$i2[3,1])) #A2B2C2
  return(fitness)
}

calc.fitness.pd.puncepi <- function(gen){
  fitness <- as.data.frame(matrix(ncol = 8, nrow = 1))
  colnames(fitness) <- combo
  fitness[1] <- mean(c((gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i1[2,2]))),
                       (gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i1[1,2]))), 
                       (gen$i1[3,1]))) #A1B1C1
  fitness[2] <- mean(c((gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i1[2,2]))),
                       (gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i1[1,2]))), 
                       (gen$i2[3,1]))) #A1B1C2
  fitness[3] <- mean(c((gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i2[2,2]))),
                       (gen$i2[2,1] * (1 - abs(gen$i2[2,2] - gen$i1[1,2]))), 
                       (gen$i1[3,1]))) #A1B2C1
  fitness[4] <- mean(c((gen$i1[1,1] * (1 - abs(gen$i1[1,2] - gen$i2[2,2]))),
                       (gen$i2[2,1] * (1 - abs(gen$i2[2,2] - gen$i1[1,2]))), 
                       (gen$i2[3,1]))) #A1B2C2
  fitness[5] <- mean(c((gen$i2[1,1] * (1 - abs(gen$i2[1,2] - gen$i1[2,2]))),
                       (gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i2[1,2]))), 
                       (gen$i1[3,1]))) #A2B1C1
  fitness[6] <- mean(c((gen$i2[1,1] * (1 - abs(gen$i2[1,2] - gen$i1[2,2]))),
                       (gen$i1[2,1] * (1 - abs(gen$i1[2,2] - gen$i2[1,2]))), 
                       (gen$i2[3,1]))) #A2B1C2
  fitness[7] <- mean(c(gen$i2[1,1],
                       gen$i2[2,1], 
                       gen$i1[3,1])) #A2B2C1
  fitness[8] <- mean(c(gen$i2[1,1],
                       gen$i2[2,1], 
                       gen$i2[3,1])) #A2B2C2
  return(fitness)
}


genotype.punc.fitness <- list("abc" = NA, "ab" = NA, "c" = NA)
for(i in seq_along(gen.lists.punc)){
  genotype.punc.fitness[[i]] <- map_dfr(gen.lists.punc[[i]], calc.fitness)
}

genotype.diff.fitness <- list("abc" = NA, "ab" = NA, "c" = NA)
for(i in seq_along(gen.lists.diff)){
  genotype.diff.fitness[[i]] <- map_dfr(gen.lists.diff[[i]], calc.fitness)
}

genotype.diff.noepi.fitness <- list("abc" = NA, "ab" = NA, "c" = NA)
for(i in seq_along(gen.lists.diff)){
  genotype.diff.noepi.fitness[[i]] <- map_dfr(gen.lists.diff[[i]], calc.fitness.diff.noepi)
}

genotype.pd.fitness <- list("abc" = NA, "ab" = NA, "c" = NA)
for(i in seq_along(gen.lists.pd)){
  genotype.pd.fitness[[i]] <- map_dfr(gen.lists.pd[[i]], calc.fitness)
}

genotype.pd.noepi.fitness <- list("abc" = NA, "ab" = NA, "c" = NA)
for(i in seq_along(gen.lists.pd)){
  genotype.pd.noepi.fitness[[i]] <- map_dfr(gen.lists.pd[[i]], calc.fitness.pd.noepi)
}

genotype.pd.midepi.fitness <- list("abc" = NA, "ab" = NA, "c" = NA)
for(i in seq_along(gen.lists.pd)){
  genotype.pd.midepi.fitness[[i]] <- map_dfr(gen.lists.pd[[i]], calc.fitness.pd.midepi)
}

genotype.pd.puncepi.fitness <- list("abc" = NA, "ab" = NA, "c" = NA)
for(i in seq_along(gen.lists.pd)){
  genotype.pd.puncepi.fitness[[i]] <- map_dfr(gen.lists.pd[[i]], calc.fitness.pd.puncepi)
}

fitness.list <- list("punc" = genotype.punc.fitness, "diff" = genotype.diff.fitness,
                     "diff_noepi" = genotype.diff.fitness, "pd_epi" = genotype.pd.fitness,
                     "pd_noepi" = genotype.pd.noepi.fitness, "pd_midepi" = genotype.pd.midepi.fitness,
                     "pd_puncepi" = genotype.pd.puncepi.fitness)
save(fitness.list, file = "runs/run_010725/rdata/fitness_list.RData")

#calculate weighted fitnesses per population
calc.weighted.fitness <- function(fit, counts){
  weighted.fit <- list()
  for(i in seq_along(rownames(fit))){
    weighted.fit[[i]] <- weighted.mean(fit[i,], counts[i,])
  }
  weighted.fit <- unlist(weighted.fit)
  return(weighted.fit)
}

apply.weighted.fitness <- function(fitnesses, counts){
  weighted.list <- list("mu0" = NA, "mu01" = NA, "mu05" = NA, "mu09" = NA, "mu1" = NA)
  for(i in seq_along(counts)){
    weighted.list[[i]] <- map2(fitnesses,counts[[i]], calc.weighted.fitness)
  }
  return(weighted.list)
}

weighted.fitness.list <- map2(fitness.list, counts.list, apply.weighted.fitness)

#per run mu1/mu0 fitness
comp.fitness <- function(weighted.fit){
  ratio.abc <- (weighted.fit$mu1$abc) / (weighted.fit$mu0$abc)
  ratio.ab <- (weighted.fit$mu1$ab) / (weighted.fit$mu0$ab)
  ratio.c <- (weighted.fit$mu1$c) / (weighted.fit$mu0$c)
  ratios <- as.data.frame(cbind("A = B = C" = ratio.abc, "AB > C" = ratio.ab, 
                                "C > AB" = ratio.c))
  ratios.gather <- gather(ratios, value = "ratio", key = "scheme")
  return(ratios.gather)
}

fitness.comps <- map(weighted.fitness.list, comp.fitness)

#plot fitness differences
comp.plots <- list()
for(i in seq_along(fitness.comps)){
  titles <- c("Punctuated", "Diffuse", "Diffuse + no epistasis", 
              "Punctuated x diffuse", "Punctuated x diffuse, no epistasis",
              "Punctuated x diffuse, intermediate epistasis",
              "Punctuated x diffuse, full epistasis")
  comp.plots[[i]] <- 
    ggplot(data = fitness.comps[[i]], aes(x = scheme, y = ratio))+
    #geom_point(cex = 4, alpha = 0.2)+
    geom_violin()+
    #ylim(c(0.2, 2.2))+
    geom_hline(yintercept = 1, linetype = "dotted")+
    xlab("Scheme")+
    ylab("Fitness ratio between \n mu = 1 and mu = 0")+
    ggtitle(titles[[i]])+
    theme_bw()+
    theme(axis.text = element_text(size = 19),
          axis.title = element_text(size = 22),
          plot.title = element_text(size = 22))
  filenames <- paste0("runs/run_010725/figs/fitness_comps/", 
                      names(weighted.fitness.list), ".png")
  ggsave(filenames[[i]])
}


#comparison statistics
dunntest.fitness <- function(comp){
  dunn <- dunn.test(x = comp$ratio, g = comp$scheme, method = "bonferroni")
  dunn.mat <- matrix(dunn$P.adjusted, nrow = 3)
  rownames(dunn.mat) <- dunn$comparisons
  colnames(dunn.mat) <- "p_adj"
  return(dunn.mat)
}

comp.stats <- map(fitness.comps, dunntest.fitness)
filenames <- paste0("runs/run_010725/figs/fitness_comps/stats/dunn_",names(fitness.comps),
                    ".csv")
for(i in seq_along(comp.stats)){
  write.csv(comp.stats[[i]], file = filenames[[i]])
}



