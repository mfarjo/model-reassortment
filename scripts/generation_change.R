library(tidyverse)

#change in per-genotype gene contribution between parental generation and offspring
load("runs/run_010725/rdata/counts_list.RData")
load("runs/run_010725/rdata/fitness_list.RData")

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

freq.check.2 <- function(offspring.count, gen.fit){
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
  parent.diff <- c("parent1" = unname(gen.fit[1]),
                   "parent2" = unname(gen.fit[8]))
  
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

calc.fx <- function(offspring.count, gen.fit){
  fx <- list("mu0"=NA,"mu01"=NA,"mu05"=NA,"mu09"=NA,"mu1"=NA)
  for(i in seq_along(offspring.count)){
    fx[[i]] <- map2(offspring.count[[i]],gen.fit,freq.table)
  }
  return(fx)
}
fx.list <- map2(counts.list, fitness.list, calc.fx)


mu0.list <- map(fx.list, function(x) x$mu0)
mu1.list <- map(fx.list, function(x) x$mu1)

cor.mu <- function(mu.list){
  cor <- list()
  cor.p <- list()
  cor.r <- list()
  cor.results <- list()
  for(i in seq_along(mu.list)){
    cor[[i]] <- suppressWarnings(cor.test(mu.list[[i]]$parental_diff, mu.list[[i]]$offspring_fx, method = "spearman"))
    cor.p[[i]] <- cor[[i]]$p.value
    cor.r[[i]] <- cor[[i]]$estimate
    cor.results[[i]] <- as.data.frame(matrix(nrow = 2, ncol = 1))
    rownames(cor.results[[i]]) <- c("correlation","p_value")
    cor.results[[i]][,1] <- c(cor.r[[i]], cor.p[[i]])
  }
  names(cor.results) <- c("abc","ab","c")
  return(cor.results)
}

cor.mu0 <- map(mu0.list, cor.mu)
cor.mu1 <- map(mu1.list, cor.mu)

combine.cors <- function(mu0,mu1){
  for(i in seq_along(mu0)){
    colnames(mu0[[i]]) <- "mu0"
  }
  for(i in seq_along(mu1)){
    colnames(mu1[[i]]) <- "mu1"
  }
  
  combined <- map2(mu0, mu1, function(x,y) bind_cols(x,y))
  return(combined)
}

cor.list <- map2(cor.mu0, cor.mu1, combine.cors)

filenames <- list()
for(i in seq_along(cor.list)){
  filenames[[i]] <- paste0("runs/run_010725/figs/generation_comps/stats/",
    names(cor.list)[i], c("_abc.csv","_ab.csv","_c.csv"))
}

save.cors <- function(cor, filename){
  for(i in seq_along(cor)){
    write.csv(cor[[i]], file = filename[[i]])
  }
}

map2(cor.list, filenames, save.cors)

plot.change <- function(mu){
  mu.plots <- list()
  for(i in seq_along(mu)){
    mu.plots[[i]] <- ggplot(data = mu[[i]], 
                            aes(x = parental_diff, y = offspring_fx, color = type))+
      geom_point(cex = 3)+
      xlab("Fitness(genotype i) - fitness(genotype j)")+
      ylab("∆ Population frequency of i")+
      scale_color_manual(values = c("black","red"))+
      theme_bw()+
      theme(axis.title = element_text(size = 22),
            axis.text = element_text(size = 19),
            plot.title = element_text(size = 22),
            legend.position = "none")
  }
  return(mu.plots)
}
mu0.plots <- map(mu0.list, plot.change)
mu1.plots <- map(mu1.list, plot.change)

header.list <- c("Punctuated", "Diffuse", "Diffuse + no epistasis", 
                 "Punctuated x diffuse", "Punctuated x diffuse, no epistasis",
                 "Punctuated x diffuse, intermediate epistasis",
                 "Punctuated x diffuse, full epistasis")
titles <- c("A = B = C", "AB > C", "C > AB")

title.change <- function(mu, header){
  for(i in seq_along(mu)){
    mu[[i]] <- mu[[i]] + ggtitle(paste0(header, "\n", titles[[i]]))
  }
  return(mu)
}

mu0.plots <- map2(mu0.plots, header.list, title.change)
mu1.plots <- map2(mu1.plots, header.list, title.change)


prefixes.mu0 <- paste0("runs/run_010725/figs/generation_comps/","comp_mu0_",names(mu0.plots),"_")
prefixes.mu1 <- paste0("runs/run_010725/figs/generation_comps/","comp_mu1_",names(mu0.plots),"_")

save.change <- function(mu, prefix){
  suffix <- c("abc.png","ab.png","c.png")
  filenames <- paste0(prefix,suffix)
  for(i in seq_along(mu)){
    ggsave(plot = mu[[i]], filename = filenames[[i]])
  }
}

map2(mu0.plots, prefixes.mu0, save.change)
map2(mu1.plots, prefixes.mu1, save.change)



