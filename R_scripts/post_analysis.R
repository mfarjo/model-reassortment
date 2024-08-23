library(tidyverse)

load("runs/run_080224/mu0.RData")
load("runs/run_080224/mu01.RData")
load("runs/run_080224/mu05.RData")
load("runs/run_080224/mu09.RData")
load("runs/run_080224/mu1.RData")


#example trajectories
mu0.test <- results.mu0[[1]][,c(1,8:15)]
mu0.gather <- gather(mu0.test, key = "genotype", value = "count", -t)

mu1.test <- results.mu1[[1]][,c(1,8:15)]
mu1.gather <- gather(mu1.test, key = "genotype", value = "count", -t)

gen.cols <- RColorBrewer::brewer.pal(11,"RdBu")[c(1:4,8:11)]
ggplot(data = mu0.gather, aes(x = t, y = count, color = genotype))+
  geom_line()+
  ylab("N")+
  scale_color_manual(values = gen.cols)
ggsave("runs/run_080224/figs/exampletraj_mu0.png")
ggplot(data = mu1.gather, aes(x = t, y = count, color = genotype))+
  geom_line()+
  ylab("N")+
  scale_color_manual(values = gen.cols)
ggsave("runs/run_080224/figs/exampletraj_mu1.png")

#final population sizes
results.mu0 <- map(results.mu0, function(x) x[,c(1,8:15)])
results.mu01 <- map(results.mu01, function(x) x[,c(1,8:15)])
results.mu05 <- map(results.mu05, function(x) x[,c(1,8:15)])
results.mu09 <- map(results.mu09, function(x) x[,c(1,8:15)])
results.mu1 <- map(results.mu1, function(x) x[,c(1,8:15)])

final.mu0 <- map(results.mu0, function(x) x[dim(x)[1],])
final.mu0.total <- map_dbl(final.mu0, function(x) sum(x[,-1]))
final.mu01 <- map(results.mu01, function(x) x[dim(x)[1],])
final.mu01.total <- map_dbl(final.mu01, function(x) sum(x[,-1]))
final.mu05 <- map(results.mu05, function(x) x[dim(x)[1],])
final.mu05.total <- map_dbl(final.mu05, function(x) sum(x[,-1]))
final.mu09 <- map(results.mu09, function(x) x[dim(x)[1],])
final.mu09.total <- map_dbl(final.mu09, function(x) sum(x[,-1]))
final.mu1 <- map(results.mu1, function(x) x[dim(x)[1],])
final.mu1.total <- map_dbl(final.mu1, function(x) sum(x[,-1]))

final.sums <- bind_cols("mu0"=final.mu0.total,"mu01"=final.mu01.total,"mu05"=final.mu05.total,
                      "mu09"=final.mu09.total, "mu1"=final.mu1.total)

final.sums.gather <- gather(final.sums, key = "mu", value = "N")

mu.cols <- c("black", RColorBrewer::brewer.pal(9,"YlOrRd")[c(3,5,7,9)])
ggplot(final.sums.gather, aes(x = N, fill = mu))+
  geom_histogram(bins = 60)+
  scale_fill_manual(values = mu.cols)+
  ylab("Count")+
  ggtitle("Total population size at t = 10")+
  theme_bw()
ggsave("runs/run_080224/figs/pop_hist.png")


#calculate fitnesses
load("runs/run_080224/gen_list.RData")
load("runs/run_080224/combo.RData")
#define fitness weighting vector
#choose between:
#c(0.3333333, 0.3333333, 0.3333333) [A = B = C]
#c(0.476, 0.476, 0.0476) [AB > C]
#c(0.0833, 0.0833, 0.833) [C > AB]
fit.vec <- c(0.3333333, 0.3333333, 0.3333333)
fit.calc <- function(geno){
  combo.mat <- as.data.frame(matrix(data = NA, ncol = 3, nrow = 8))
  colnames(combo.mat) <- c("A","B","C")
  combo.mat[1,] <- c(geno[1,1],geno[2,1],geno[3,1])
  combo.mat[2,] <- c(geno[1,1],geno[2,1],geno[3,2])
  combo.mat[3,] <- c(geno[1,1],geno[2,2],geno[3,1])
  combo.mat[4,] <- c(geno[1,1],geno[2,2],geno[3,2])
  combo.mat[5,] <- c(geno[1,2],geno[2,1],geno[3,1])
  combo.mat[6,] <- c(geno[1,2],geno[2,1],geno[3,2])
  combo.mat[7,] <- c(geno[1,2],geno[2,2],geno[3,1])
  combo.mat[8,] <- c(geno[1,2],geno[2,2],geno[3,2])
  combo.mat <- mutate(combo.mat, "id" = combo, .before = A)
  
  w <- list()
  for(i in seq_along(combo.mat$id)){
    w[[i]] <- sum((fit.vec[1]) * ((combo.mat[i,]$A) * (1-abs(combo.mat[i,]$A - combo.mat[i,]$B))),
                  (fit.vec[2]) * ((combo.mat[i,]$B) * (1-abs(combo.mat[i,]$B - combo.mat[i,]$A))),
                  (fit.vec[3]) * (combo.mat[i,]$C))
  }
  w <- unlist(w)
  combo.mat <- mutate(combo.mat, "w" = w)
  combo.mat <- as.data.frame(t(combo.mat))
  colnames(combo.mat) <- combo
  combo.mat <- combo.mat[-c(1:4),]
  return(combo.mat)
}

all.w <- map(gen.list,fit.calc)

fit.mu0 <- map2(all.w,final.mu0, function(x,y) rbind(x,y[,-1]))
fit.mu01 <- map2(all.w,final.mu01, function(x,y) rbind(x,y[,-1]))
fit.mu05 <- map2(all.w,final.mu05, function(x,y) rbind(x,y[,-1]))
fit.mu09 <- map2(all.w,final.mu09, function(x,y) rbind(x,y[,-1]))
fit.mu1 <- map2(all.w,final.mu1, function(x,y) rbind(x,y[,-1]))

weight.mean <- function(fit){
  a <- list()
  for(i in seq_along(fit[1,])){
    a[[i]] <- sum(as.numeric(fit[1,i]) * as.numeric(fit[2,i]))
    b <- sum(unlist(a))
    c <- b/sum(as.numeric(fit[2,]))
  }
  return(c)
}
means.mu0 <- map_dbl(fit.mu0,weight.mean)
means.mu01 <- map_dbl(fit.mu01,weight.mean)
means.mu05 <- map_dbl(fit.mu05,weight.mean)
means.mu09 <- map_dbl(fit.mu09,weight.mean)
means.mu1 <- map_dbl(fit.mu1,weight.mean)

all.means <- bind_cols("mu0" = means.mu0, "mu01" = means.mu01, "mu05" = means.mu05,
                       "mu09" = means.mu09,"mu1" = means.mu1)
means.gather <- gather(all.means, value = "mean", key = "mu")

ggplot(means.gather, aes(x = mean, color = mu))+
  #geom_histogram(bins = 50)+
  geom_density()+
  scale_color_manual(values = mu.cols)+
  xlab("Fitness")+
  ylab("Count")+
  ggtitle("Weighted mean of fitness at t = 10")+
  theme_bw()
ggsave("runs/run_080224/figs/w_density.png")



