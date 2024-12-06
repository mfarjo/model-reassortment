library(tidyverse)
library(ggpubr)

#load replication statistics
load("runs/run_100124/rdata/rep_stats_init1_punc.RData")
load("runs/run_100124/rdata/rep_stats_init5_punc.RData")
load("runs/run_100124/rdata/rep_stats_init1_diff.RData")
load("runs/run_100124/rdata/rep_stats_init5_diff.RData")

rep.list <- list("init1_punc" = rep.stats.init1.punc, "init5_punc" = rep.stats.init5.punc,
                 "init1_diff" = rep.stats.init1.diff, "init5_diff" = rep.stats.init5.diff)

#calculate total viable genomes post-replication
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

tab.viables <- function(repz){
  viable.mu0 <- map(repz$mu0, count.genomes.mu0)
  viable.mu01 <- map(repz$mu01, count.genomes.muplus)
  viable.mu05 <- map(repz$mu05, count.genomes.muplus)
  viable.mu09 <- map(repz$mu09, count.genomes.muplus)
  viable.mu1 <- map(repz$mu1, count.genomes.muplus)
  
  viable.list <- list("mu0"= viable.mu0, "mu01" = viable.mu01, "mu05" = viable.mu05,
                      "mu09" = viable.mu09, "mu1" = viable.mu1)
  return(viable.list)
}

viable.list <- map(rep.list, tab.viables)

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

viable.tables <- map(viable.list, build.viable.tables)
viables.punc <- list("init1" = viable.tables$init1_punc, "init5" = viable.tables$init5_punc)
viables.diff <- list("init1" = viable.tables$init1_diff, "init5" = viable.tables$init5_diff)

#plot total viable genomes
plot.viables.punc <- function(viable){
  mu.cols <- c("black", RColorBrewer::brewer.pal(9,"YlOrRd")[c(3,5,7,9)])
  viable.plots.init1 <- list()
  viable.plots.init5 <- list()
  for(i in seq_along(viable$init1)){
    titles <- c("A = B = C, 1 founder", "AB > C, 1 founder", "C > AB, 1 founder")
    viable.plots.init1[[i]] <- 
      ggplot(data = viable$init1[[i]], aes(x = count, fill = mu))+
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
    filenames <- paste0("runs/run_100124/figs/viable_genomes/",
                        c("viable_genomes_punc_init1_abc.png", "viable_genomes_punc_init1_ab.png",
                          "viable_genomes_punc_init1_c.png"))
    ggsave(filenames[[i]])
  }
  for(i in seq_along(viable$init5)){
    titles <- c("A = B = C, 5 founders", "AB > C, 5 founders", "C > AB, 5 founders")
    viable.plots.init1[[i]] <- 
      ggplot(data = viable$init5[[i]], aes(x = count, fill = mu))+
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
    filenames <- paste0("runs/run_100124/figs/viable_genomes/",
                        c("viable_genomes_punc_init5_abc.png", "viable_genomes_punc_init5_ab.png",
                          "viable_genomes_punc_init5_c.png"))
    ggsave(filenames[[i]])
  }
}

plot.viables.punc(viables.punc)

plot.viables.diff <- function(viable){
  mu.cols <- c("black", RColorBrewer::brewer.pal(9,"YlOrRd")[c(3,5,7,9)])
  viable.plots.init1 <- list()
  viable.plots.init5 <- list()
  for(i in seq_along(viable$init1)){
    titles <- c("A = B = C, 1 founder", "AB > C, 1 founder", "C > AB, 1 founder")
    viable.plots.init1[[i]] <- 
      ggplot(data = viable$init1[[i]], aes(x = count, fill = mu))+
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
    filenames <- paste0("runs/run_100124/figs/viable_genomes/",
                        c("viable_genomes_diff_init1_abc.png", "viable_genomes_diff_init1_ab.png",
                          "viable_genomes_diff_init1_c.png"))
    ggsave(filenames[[i]])
  }
  for(i in seq_along(viable$init5)){
    titles <- c("A = B = C, 5 founders", "AB > C, 5 founders", "C > AB, 5 founders")
    viable.plots.init1[[i]] <- 
      ggplot(data = viable$init5[[i]], aes(x = count, fill = mu))+
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
    filenames <- paste0("runs/run_100124/figs/viable_genomes/",
                        c("viable_genomes_diff_init5_abc.png", "viable_genomes_diff_init5_ab.png",
                          "viable_genomes_diff_init5_c.png"))
    ggsave(filenames[[i]])
  }
}

plot.viables.diff(viables.diff)


#packaging fitness stats
load("runs/run_100124/rdata/pack_stats_init5_punc.RData")
load("runs/run_100124/rdata/pack_stats_init5_diff.RData")
pack.list <- list("punc" = pack.stats.init5.punc, "diff" = pack.stats.init5.diff)

tab.final.genotypes <- function(packx){
  final.genotypes <- list("abc" = NA, "ab" = NA, "c" = NA)
  for(i in seq_along(packx)){
    final.genotypes[[i]] <- map_dfr(packx[[i]], function(x) tail(x,1)[,c(8:15)])
  }
  return(final.genotypes)
}

genotype.list <- list("punc" = NA, "diff" = NA)

for(i in seq_along(pack.list)){
  genotype.list[[i]] <- map(pack.list[[i]], tab.final.genotypes)
}

#calculate fitness for each genotype
load("runs/run_100124/rdata/gen_lists_punctuated.RData")
load("runs/run_100124/rdata/gen_lists_diffuse.RData")
load("runs/run_100124/rdata/combo.RData")

#fit.vecs <- list("abc" = c(0.333,0.333,0.333),
                 #"ab" = c(0.476,0.476,0.0476),
                 #"c" = c(0.0833,0.0833,0.833))

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

genotype.punc.fitness <- list("abc" = NA, "ab" = NA, "c" = NA)
for(i in seq_along(gen.lists.punc)){
  genotype.punc.fitness[[i]] <- map_dfr(gen.lists.punc[[i]], calc.fitness)
}

genotype.diff.fitness <- list("abc" = NA, "ab" = NA, "c" = NA)
for(i in seq_along(gen.lists.diff)){
  genotype.diff.fitness[[i]] <- map_dfr(gen.lists.diff[[i]], calc.fitness)
}

calc.weighted.fitness <- function(fit, counts){
  weighted.fit <- list()
  for(i in seq_along(rownames(fit))){
    weighted.fit[[i]] <- weighted.mean(fit[i,], counts[i,])
  }
  weighted.fit <- unlist(weighted.fit)
  return(weighted.fit)
}

weighted.fitness.punc <- list("mu0" = NA, "mu01" = NA, "mu05" = NA, "mu09" = NA, "mu1" = NA)
for(i in seq_along(genotype.list$punc)){
  weighted.fitness.punc[[i]] <- map2(genotype.punc.fitness, genotype.list$punc[[i]], calc.weighted.fitness)
}

weighted.fitness.diff <- list("mu0" = NA, "mu01" = NA, "mu05" = NA, "mu09" = NA, "mu1" = NA)
for(i in seq_along(genotype.list$diff)){
  weighted.fitness.diff[[i]] <- map2(genotype.diff.fitness, genotype.list$diff[[i]], calc.weighted.fitness)
}

weighted.fitness.list <- list("punc" = weighted.fitness.punc,
                              "diff" = weighted.fitness.diff)

build.fitness.tables <- function(weighted.fit){
  fitness.abc <- map(weighted.fit, function(x) x$abc)
  fitness.abc <- bind_rows(fitness.abc)
  fitness.abc.gather <- gather(fitness.abc, key = "mu", value = "mean")
  
  fitness.ab <- map(weighted.fit, function(x) x$ab)
  fitness.ab <- bind_rows(fitness.ab)
  fitness.ab.gather <- gather(fitness.ab, key = "mu", value = "mean")
  
  fitness.c <- map(weighted.fit, function(x) x$c)
  fitness.c <- bind_rows(fitness.c)
  fitness.c.gather <- gather(fitness.c, key = "mu", value = "mean")
  
  fitness.tables <- list("abc" = fitness.abc.gather, "ab" = fitness.ab.gather,
                         "c" = fitness.c.gather)
  return(fitness.tables)
}

fitness.tables <- list("punc" = NA, "diff" = NA)
fitness.tables <- map(weighted.fitness.list, build.fitness.tables)

mu.cols <- c("black", RColorBrewer::brewer.pal(9,"YlOrRd")[c(3,5,7,9)])

fitness.plots.punc <- list()
for(i in seq_along(fitness.tables$punc)){
  titles <- c("A = B = C", "AB > C", "C > AB")
  fitness.plots.punc[[i]] <- 
    ggplot(data = fitness.tables$punc[[i]], aes(x = mean, fill = mu))+
    geom_density(alpha = 0.6)+
    xlab("Mean population fitness")+
    ylab("Density")+
    ggtitle(paste0(titles[[i]], ", punctuated diversity"))+
    xlim(c(0.2,2))+
    scale_fill_manual(values = mu.cols, labels = c(0, 0.1, 0.5, 0.9, 1))+
    theme_bw()+
    theme(axis.text = element_text(size = 19),
          axis.title = element_text(size = 22),
          plot.title = element_text(size = 22),
          legend.text = element_text(size = 19),
          legend.title = element_text(size = 22))
  filenames <- paste0("runs/run_100124/figs/genome_fitness/",
                      c("fitness_abc_punc.png", "fitness_ab_punc.png",
                        "fitness_c_punc.png"))
  ggsave(filenames[[i]])
}


fitness.plots.diff <- list()
for(i in seq_along(fitness.tables$diff)){
  titles <- c("A = B = C", "AB > C", "C > AB")
  fitness.plots.diff[[i]] <- 
    ggplot(data = fitness.tables$diff[[i]], aes(x = mean, fill = mu))+
    geom_density(alpha = 0.6)+
    xlab("Mean population fitness")+
    ylab("Density")+
    ggtitle(paste0(titles[[i]], ", diffuse diversity"))+
    xlim(c(0.2,2))+
    scale_fill_manual(values = mu.cols, labels = c(0, 0.1, 0.5, 0.9, 1))+
    theme_bw()+
    theme(axis.text = element_text(size = 19),
          axis.title = element_text(size = 22),
          plot.title = element_text(size = 22),
          legend.text = element_text(size = 19),
          legend.title = element_text(size = 22))
  filenames <- paste0("runs/run_100124/figs/genome_fitness/",
                      c("fitness_abc_diff.png", "fitness_ab_diff.png",
                        "fitness_c_diff.png"))
  ggsave(filenames[[i]])
}


#per run mu1 - mu0 fitness
comp.fitness <- function(weighted.fit){
  diff.abc <- (weighted.fit$mu1$abc) - (weighted.fit$mu0$abc)
  diff.ab <- (weighted.fit$mu1$ab) - (weighted.fit$mu0$ab)
  diff.c <- (weighted.fit$mu1$c) - (weighted.fit$mu0$c)
  diffs <- as.data.frame(cbind("A = B = C" = diff.abc, "AB > C" = diff.ab, 
                               "C > AB" = diff.c))
  diffs.gather <- gather(diffs, value = "diff", key = "scheme")
  return(diffs.gather)
}

fitness.comps <- map(weighted.fitness.list, comp.fitness)

#stats
kruskal.test(diff ~ scheme, fitness.comps$punc) #p-value = 0.0002338
dunn.test::dunn.test(fitness.comps$punc$diff, fitness.comps$punc$scheme, 
                     method = "bonferroni")
#A = B = C vs. AB > C --> p = 0.0001

kruskal.test(diff ~ scheme, fitness.comps$diff) #p-value = 0.0096
dunn.test::dunn.test(fitness.comps$diff$diff, fitness.comps$diff$scheme, 
                     method = "bonferroni")
#A = B = C vs. AB > C --> p = 0.0066

for(i in seq_along(fitness.comps)){
  titles <- c("Punctuated diversity", "Diffuse diversity")
  comp.plots <- list()
  comp.plots[[i]] <- 
    ggplot(data = fitness.comps[[i]], aes(x = scheme, y = diff))+
    geom_point(cex = 4, alpha = 0.2)+
    ylim(c(-0.5, 0.75))+
    geom_hline(yintercept = 0, linetype = "dotted")+
    xlab("Scheme")+
    ylab("Fitness difference between \n mu = 1 and mu = 0")+
    ggtitle(titles[[i]])+
    theme_bw()+
    theme(axis.text = element_text(size = 19),
          axis.title = element_text(size = 22),
          plot.title = element_text(size = 22))
  filenames <- paste0("runs/run_100124/figs/fitness_comp_", 
                      c("punctuated.png", "diffuse.png"))
  ggsave(filenames[[i]])
}


#change in per-genotype gene contribution between parental generation and offspring
freq.check <- function(offspring.count, gen.fit){
  offspring.freq <- as.data.frame(t(offspring.count/(sum(offspring.count))))
  colnames(offspring.freq) <- "freq"
  offspring.tab <- mutate(offspring.freq,
                          "parent1.id" = c(1,0.67,0.67,0.33,0.67,0.33,0.33,0),
                          "parent2.id" = c(0, 0.33,0.33,0.67,0.33,0.67,0.67,1))
  offspring.tab <- mutate(offspring.tab,
                          "freqxid.1" = freq * parent1.id,
                          "freqxid.2" = freq * parent2.id)
  
  offspring.fx <- colSums(offspring.tab)[c(4,5)] - 0.5
  names(offspring.fx) <- c("parent1", "parent2")
  parent.diff <- c("parent1" = unname(gen.fit[1] - gen.fit[8]),
                   "parent2" = unname(gen.fit[8] - gen.fit[1]))
  
  freq.mat <- as.data.frame(bind_rows(offspring.fx, parent.diff))
  rownames(freq.mat) <- c("offspring_fx","parental_diff")
  return(freq.mat)
}

freq.table <- function(offspring.count, gen.fit){
  comp <- list()
  for(i in seq_along(rownames(offspring.count))){
    comp[[i]] <- freq.check(offspring.count[i,], gen.fit[i,])
  }
  comp <- map(comp, function(x) t(x)[1,])
  comp <- bind_rows(comp)
  comp <- mutate(comp,"type" = "mismatch")
  for(i in seq_along(comp$offspring_fx)){
    if((comp$offspring_fx[[i]] < 0) && (comp$parental_diff[[i]] < 0)){
      comp$type[[i]] <- "match"
    }
    if((comp$offspring_fx[[i]] >= 0) && (comp$parental_diff[[i]] >= 0)){
      comp$type[[i]] <- "match"
    }
  }
  return(comp)
}

abc.freqs.punc <- list("mu0" = freq.table(genotype.list$punc$mu0$abc, genotype.punc.fitness$abc),
                       "mu1" = freq.table(genotype.list$punc$mu1$abc, genotype.punc.fitness$abc))
ab.freqs.punc <- list("mu0" = freq.table(genotype.list$punc$mu0$ab, genotype.punc.fitness$ab),
                      "mu1" = freq.table(genotype.list$punc$mu1$ab, genotype.punc.fitness$ab))
c.freqs.punc <- list("mu0" = freq.table(genotype.list$punc$mu0$c, genotype.punc.fitness$c),
                     "mu1" = freq.table(genotype.list$punc$mu1$c, genotype.punc.fitness$c))
freqs.punc <- list("A = B = C" = abc.freqs.punc,"AB > C" = ab.freqs.punc,
                   "C > AB" = c.freqs.punc)


abc.freqs.diff <- list("mu0" = freq.table(genotype.list$diff$mu0$abc, genotype.diff.fitness$abc),
                       "mu1" = freq.table(genotype.list$diff$mu1$abc, genotype.diff.fitness$abc))
ab.freqs.diff <- list("mu0" = freq.table(genotype.list$diff$mu0$ab, genotype.diff.fitness$ab),
                      "mu1" = freq.table(genotype.list$diff$mu1$ab, genotype.diff.fitness$ab))
c.freqs.diff <- list("mu0" = freq.table(genotype.list$diff$mu0$c, genotype.diff.fitness$c),
                     "mu1" = freq.table(genotype.list$diff$mu1$c, genotype.diff.fitness$c))
freqs.diff <- list("A = B = C" = abc.freqs.diff,"AB > C" = ab.freqs.diff,
                   "C > AB" = c.freqs.diff)


plot.freqs <- function(freqs){
  titles <- c("mu = 0", "mu = 1")
  freq.plots <- list()
  cor.list <- list()
  lm.list <- list()
  for(i in seq_along(freqs)){
    cor.list[[i]] <- (cor.test(freqs[[i]]$parental_diff, freqs[[i]]$offspring_fx))$estimate
    lm.list[[i]] <- (lm(freqs[[i]]$offspring_fx ~ freqs[[i]]$parental_diff))$coefficients
    freq.plots[[i]] <- ggplot(data = freqs[[i]], 
                              aes(x = parental_diff, y = offspring_fx, color = type))+
      geom_point(cex = 3)+
      geom_abline(intercept = lm.list[[i]][1], slope = lm.list[[i]][2])+
      xlim(c(-0.8,0.8))+
      ylim(c(-0.5,0.5))+
      xlab("Fitness(genotype i) - fitness(genotype j)")+
      ylab("∆ Population frequency of i")+
      ggtitle(paste(titles[[i]], ", r =",round(cor.list[[i]], 2)))+
      scale_color_manual(values = c("black","red"))+
      theme_bw()+
      theme(axis.title = element_text(size = 22),
            axis.text = element_text(size = 19),
            plot.title = element_text(size = 22),
            legend.position = "none")
  }
  return(freq.plots)
}

freq.plots.punc <- map(freqs.punc, plot.freqs)
freq.plots.punc <- flatten(freq.plots.punc)
names(freq.plots.punc) <- c("abc_mu0","abc_mu1","ab_mu0","ab_mu1","c_mu0","c_mu1")
filenames.punc <- paste0("runs/run_100124/figs/freqs/freqs_punc_",
                    names(freq.plots.punc),".png")

for(i in seq_along(freq.plots.punc)){
  plot(freq.plots.punc[[i]])
  ggsave(filenames.punc[[i]])
}

freq.plots.diff <- map(freqs.diff, plot.freqs)
freq.plots.diff <- flatten(freq.plots.diff)
names(freq.plots.diff) <- c("abc_mu0","abc_mu1","ab_mu0","ab_mu1","c_mu0","c_mu1")
filenames.diff <- paste0("runs/run_100124/figs/freqs/freqs_diff_",
                         names(freq.plots.diff),".png")

for(i in seq_along(freq.plots.diff)){
  plot(freq.plots.diff[[i]])
  ggsave(filenames.diff[[i]])
}










#older
abc.freqs <- list("mu0" = freq.table(genotypes.init5$mu0$abc, genotype.fitness$abc),
                   "mu1" = freq.table(genotypes.init5$mu1$abc, genotype.fitness$abc))
ab.freqs <- list("mu0" = freq.table(genotypes.init5$mu0$ab, genotype.fitness$ab),
                  "mu1" = freq.table(genotypes.init5$mu1$ab, genotype.fitness$ab))
c.freqs <- list("mu0" = freq.table(genotypes.init5$mu0$c, genotype.fitness$c),
                 "mu1" = freq.table(genotypes.init5$mu1$c, genotype.fitness$c))





plot.freqs <- function(freqs){
  titles <- c("mu = 0", "mu = 1")
  freq.plots <- list()
  cor.list <- list()
  lm.list <- list()
  for(i in seq_along(freqss)){
    cor.list[[i]] <- (cor.test(freqs[[i]]$parental_diff, freqs[[i]]$offspring_fx))$estimate
    lm.list[[i]] <- (lm(freqs[[i]]$offspring_fx ~ freqs[[i]]$parental_diff))$coefficients
    freq.plots[[i]] <- ggplot(data = freqs[[i]], 
                               aes(x = parental_diff, y = offspring_fx, color = type))+
      geom_point(cex = 3)+
      geom_abline(intercept = lm.list[[i]][1], slope = lm.list[[i]][2])+
      xlim(c(-0.5,0.5))+
      ylim(c(-0.5,0.5))+
      xlab("Fitness(genotype i) - fitness(genotype j)")+
      ylab("∆ Population frequency of i")+
      ggtitle(paste(titles[[i]], ", r =",round(cor.list[[i]], 2)))+
      scale_color_manual(values = c("black","red"))+
      theme_bw()+
      theme(axis.title = element_text(size = 22),
            axis.text = element_text(size = 19),
            plot.title = element_text(size = 22),
            legend.position = "none")
    
  }
  return(freq.plots)
}

abc.freq.plots <- plot.freqs(abc.freqs)
abc.freq.plots[[2]]
ggsave("runs/run_091024/figs/abc_price_mu1.png")
ab.freq.plots <- plot.freqs(ab.freqs)
ab.freq.plots[[2]]
ggsave("runs/run_091024/figs/ab_price_mu1.png")
c.freq.plots <- plot.freqs(c.freqs)
c.freq.plots[[2]]
ggsave("runs/run_091024/figs/c_price_mu1.png")


