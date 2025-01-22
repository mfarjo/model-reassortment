library(tidyverse)
library(ggpubr)

#load viable genome stats
load("runs/run_010725/rdata/viables_punc.RData")
load("runs/run_010725/rdata/viables_diff.RData")
load("runs/run_010725/rdata/viables_diff_noepi.RData")
load("runs/run_010725/rdata/viables_pd_epi.RData")
load("runs/run_010725/rdata/viables_pd_noepi.RData")
load("runs/run_010725/rdata/viables_pd_midepi.RData")
load("runs/run_010725/rdata/viables_pd_puncepi.RData")

viables.list <- list("punc" = viables.punc, "diff" = viables.diff, 
                     "diff_noepi" = viables.diff.noepi, "pd_epi" = viables.pd.epi,
                     "pd_noepi" = viables.pd.noepi, "pd_midepi" = viables.pd.midepi,
                     "pd_puncepi" = viables.pd.puncepi)

gather.viables <- function(viable){
  via.abc <- gather(bind_cols(map(viable, function(x) x$abc)), key = "mu", value = "viables")
  via.ab <- gather(bind_cols(map(viable, function(x) x$ab)), key = "mu", value = "viables")
  via.c <- gather(bind_cols(map(viable, function(x) x$c)), key = "mu", value = "viables")
  via.list <- list("abc" = via.abc,"ab" = via.ab, "c" = via.c)
  return(via.list)
}
viables.list <- map(viables.list, gather.viables)

#viable genome violin plots
plot.violin <- function(viable){
  vio.plots <- list()
  for(i in seq_along(viable)){
    mu.cols <- c("black", RColorBrewer::brewer.pal(9,"YlOrRd")[c(3,5,7,9)])
    vio.plots[[i]] <- ggplot(data = viable[[i]], aes(x = mu, y = viables, fill = mu))+
      geom_violin()+
      ylab("Viable genomes")+
      xlab("\U003BC")+
      scale_x_discrete(labels=c("mu0" = "0.0", "mu01" = "0.1","mu05" = "0.5", 
                                "mu09" = "0.9","mu1" = "1.0"))+
      scale_fill_manual(values = mu.cols)+
      theme_bw()+
      theme(axis.text = element_text(size = 19),
            axis.title = element_text(size = 22),
            plot.title = element_text(size = 22),
            legend.position = "none")
  }
  return(vio.plots)
}

viable.violins <- map(viables.list, plot.violin)

header.list <- c("Punctuated", "Diffuse", "Diffuse + no epistasis", 
                 "Punctuated x diffuse", "Punctuated x diffuse, no epistasis",
                 "Punctuated x diffuse, intermediate epistasis",
                 "Punctuated x diffuse, full epistasis")

title.violin <- function(viable, header){
  titles <- c("A = B = C", "AB > C", "C > AB")
  titled.plots <- list()
  for(i in seq_along(viable)){
    titled.plots[[i]] <- viable[[i]] + ggtitle(paste0(header,"\n",titles[[i]]))
  }
  return(titled.plots)
}

viable.violins <- map2(viable.violins, header.list,title.violin)

prefixes <- paste0("runs/run_010725/figs/violins/",names(viables.list),"_")

save.violin <- function(viable, prefix){
  suffix <- c("abc.png","ab.png","c.png")
  filenames <- paste0(prefix,suffix)
  for(i in seq_along(viable)){
    ggsave(plot = viable[[i]], filename = filenames[[i]])
  }
}

map2(viable.violins, prefixes, save.violin)

#calculate mu vs viability correlations
calc.viable.cor <- function(viable){
 cor.mat <- matrix(data = NA, nrow = 3, ncol = 2)
 rownames(cor.mat) <- c("abc","ab","c")
 colnames(cor.mat) <- c("rho","p")
 for(i in seq_along(viable)){
   viable[[i]]$mu[which(viable[[i]]$mu == "mu0")] <- 0
   viable[[i]]$mu[which(viable[[i]]$mu == "mu01")] <- 0.1
   viable[[i]]$mu[which(viable[[i]]$mu == "mu05")] <- 0.5
   viable[[i]]$mu[which(viable[[i]]$mu == "mu09")] <- 0.9
   viable[[i]]$mu[which(viable[[i]]$mu == "mu1")] <- 1
   viable[[i]]$mu <- as.numeric(viable[[i]]$mu)
   cor.mat[i,] <- c((cor.test(viable[[i]]$mu, viable[[i]]$viables, method = "spearman"))$estimate,
                    (cor.test(viable[[i]]$mu, viable[[i]]$viables, method = "spearman"))$p.value)
 }
 return(cor.mat)
}

viable.cors <- map(viables.list, calc.viable.cor)

filenames <- paste0("runs/run_010725/figs/violins/stats/cor_",names(viable.cors),".csv")

for(i in seq_along(viable.cors)){
  write.csv(viable.cors[[i]], filenames[[i]])
}

