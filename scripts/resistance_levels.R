library(tidyverse)
library(ggpubr)
library(dunn.test)

#packaging fitness stats (resistance levels)
load("runs/run_010725/rdata/final_gen_low.RData")
load("runs/run_010725/rdata/final_gen_punc.RData")
load("runs/run_010725/rdata/final_gen_high.RData")

res.counts.list <- list("low" = final.gen.low, "mid" = final.gen.punc,"high" = final.gen.high)

#calculate fitness for each genotype
load("runs/run_010725/rdata/gen_lists_lowres.RData")
load("runs/run_010725/rdata/gen_lists_punctuated.RData")
load("runs/run_010725/rdata/gen_lists_highres.RData")
load("runs/run_010725/rdata/combo.RData")
res.fitness.list <- list("low" = gen.lists.low, "mid" = gen.lists.punc, "high" = gen.lists.high)

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

apply.fitness <- function(gen){
  genotype.fitness <- list("abc" = NA, "ab" = NA, "c" = NA)
  for(i in seq_along(gen)){
    genotype.fitness[[i]] <- map_dfr(gen[[i]], calc.fitness)
  }
  return(genotype.fitness)
}

res.fitness.list <- map(res.fitness.list, apply.fitness)


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

res.weighted.fitness.list <- map2(res.fitness.list, res.counts.list, apply.weighted.fitness)


#per run mu1/mu0 fitness
comp.fitness <- function(weighted.fit){
  ratio.abc <- (weighted.fit$mu1$abc) / (weighted.fit$mu0$abc)
  ratio.ab <- (weighted.fit$mu1$ab) / (weighted.fit$mu0$ab)
  ratio.c <- (weighted.fit$mu1$c) / (weighted.fit$mu0$c)
  ratios <- as.data.frame(cbind("A = B = C" = ratio.abc, "AB > C" = ratio.ab, 
                                "C > AB" = ratio.c))
  return(ratios)
}

res.fitness.comps <- map(res.weighted.fitness.list, comp.fitness)
res.abcs <- bind_cols("0.05" = res.fitness.comps$low$`A = B = C`,
                      "0.10" = res.fitness.comps$mid$`A = B = C`,
                      "0.20" = res.fitness.comps$high$`A = B = C`) 
res.abs <- bind_cols("0.05" = res.fitness.comps$low$`AB > C`,
                      "0.10" = res.fitness.comps$mid$`AB > C`,
                      "0.20" = res.fitness.comps$high$`AB > C`) 
res.cs <- bind_cols("0.05" = res.fitness.comps$low$`C > AB`,
                     "0.10" = res.fitness.comps$mid$`C > AB`,
                     "0.20" = res.fitness.comps$high$`C > AB`)

#AB > C
res.ab.gather <- gather(res.abs, value = "fitness",key = "resistance")
ggplot(data = res.ab.gather, aes(x = resistance, y = fitness))+
  geom_violin()+
  geom_hline(yintercept = 1, linetype = "dotted")+
  xlab("Resistance level")+
  ylab("Fitness ratio between \n mu = 1 and mu = 0")+
  ggtitle("AB > C")+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 22))
ggsave("runs/run_010725/figs/resistance/ab_levels.png")

res.ab.gather$resistance <- as.numeric(res.ab.gather$resistance)
cor.test(res.ab.gather$resistance, res.ab.gather$fitness, method = "spearman")

#C > AB
res.c.gather <- gather(res.cs, value = "fitness",key = "resistance")
ggplot(data = res.c.gather, aes(x = resistance, y = fitness))+
  geom_violin()+
  geom_hline(yintercept = 1, linetype = "dotted")+
  xlab("Resistance level")+
  ylab("Fitness ratio between \n mu = 1 and mu = 0")+
  ggtitle("C > AB")+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22),
        plot.title = element_text(size = 22))
ggsave("runs/run_010725/figs/resistance/c_levels.png")
res.c.gather$resistance <- as.numeric(res.c.gather$resistance)
cor.test(res.c.gather$resistance, res.c.gather$fitness, method = "spearman")

