library(tidyverse)

load("final_project/run_040324/mu01.RData")
load("final_project/run_040324/mu03.RData")
load("final_project/run_040324/mu05.RData")
load("final_project/run_040324/mu07.RData")
load("final_project/run_040324/mu09.RData")
load("final_project/run_040324/mu1.RData")

mu.list <- list("mu_01" = mu.01, "mu_03" = mu.03, "mu_05" = mu.05, "mu_07" = mu.07, 
                "mu_09" = mu.09, "mu_1" = mu.1)
titles <- c("mu = 0.1", "mu = 0.3", "mu = 0.5", "mu = 0.7", "mu = 0.9", "mu = 1")

#plot mean trajectories
mean.traj <- function(mu, title){
  mu.mean <- map(mu, select, mean)
  mu.mean <- bind_cols(mu.mean, .name_repair = "minimal")
  colnames(mu.mean) <- names(mu)
  mu.mean <- mu.mean %>%
    mutate("cycle" = c(1:100)) %>%
    relocate(cycle, .before = sim1) %>%
    gather(key = "sim", value = "w", -cycle)
  mu.plot <- ggplot(data = mu.mean, aes(x = cycle, y = w, color = sim))+
    geom_line()+
    xlab("Generation")+
    ylab("Mean fitness")+
    ggtitle(paste("Mean population fitness,",title))+
    theme_bw()+
    scale_color_manual(values = rep("black", 100))+
    theme(legend.position = "none")
  return(mu.plot)
}
mean.traj.plots <- map2(mu.list, titles, mean.traj)

#plot max trajectories
max.traj <- function(mu, title){
  mu.max <- map(mu, select, max)
  mu.max <- bind_cols(mu.max, .name_repair = "minimal")
  colnames(mu.max) <- names(mu)
  mu.max <- mu.max %>%
    mutate("cycle" = c(1:100)) %>%
    relocate(cycle, .before = sim1) %>%
    gather(key = "sim", value = "w", -cycle)
  mu.plot <- ggplot(data = mu.max, aes(x = cycle, y = w, color = sim))+
    geom_line()+
    xlab("Generation")+
    ylab("Max fitness")+
    ggtitle(paste("Max population fitness,",title))+
    theme_bw()+
    scale_color_manual(values = rep("black", 100))+
    theme(legend.position = "none")
  return(mu.plot)
}
max.traj.plots <- map2(mu.list, titles, max.traj)

