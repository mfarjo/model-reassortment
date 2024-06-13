library(tidyverse)
library(gtools)

test.props <- as.data.frame(rdirichlet(20, c(5,1,1,1,1,1,1,5)))
colnames(test.props) <- c("A","B","C","D","E","F","G","H")

prop.gather <- gather(test.props, key = "rep", value = "freq")
ggplot(data = prop.gather, aes(x = rep, y = freq))+
  geom_point(cex = 3)+
  xlab("Genotype")+
  ylab("Frequency")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19))
ggsave("runs/run_061224/dirichlet_dist.png")

test.props <- test.props %>%
  mutate("rep" = c(1:20)) %>%
  group_split(rep)

test.props <- map(test.props, select, -rep)
test.props <- lapply(test.props, t)
test.props <- map(test.props, as.data.frame)
test.props <- map(test.props, rownames_to_column)
for(i in seq_along(test.props)){
  colnames(test.props[[i]]) <- c("ID","freq")
}

prop.sample <- function(prop, trials){
  prop.draw <- sample(prop$ID,trials,replace = T,prob = prop$freq)
  post.freq <- as.data.frame(table(prop.draw)/trials)
  colnames(post.freq) <- c("ID","post_freq")
  freq.tab <- left_join(prop, post.freq, by = "ID")
  freq.tab <- relocate(freq.tab, ID, .before = freq)
  freq.tab[is.na(freq.tab)] <- 0
  return(freq.tab)
}

re.samples <- list()
seq.vec <- seq(5,100,5)
for(i in seq_along(seq.vec)){
  re.samples[[i]] <- map2(test.props, seq.vec[[i]], prop.sample)
}
names(re.samples) <- seq(5,100,5)

calc.see <- function(re){
  num <- list()
  see <- list()
  for(i in seq_along(re)){
    num[[i]] <- sum((re[[i]]$freq - re[[i]]$post_freq)^2)
    see[[i]] <- sqrt(num[[i]]/8)
  }
  see <- unlist(see)
  return(see)
}
comp <- map(re.samples, calc.see)

comp <- bind_cols(comp)
comp.gather <- gather(comp, value = "off", key = "sample_count")
comp.gather$sample_count <- as.numeric(comp.gather$sample_count)

ggplot(comp.gather, aes(x= sample_count, y= off))+
  geom_point(cex = 3)+
  xlab("N")+
  ylab("Standard error of estimate")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19))
ggsave("runs/run_061224/dirichlet_see.png")

#or calculate standard error for a range of population frequencies
calc.se <- function(p,n){
  num <- p * (1-p)
  inner <- num/n
  se <- sqrt(inner)
  return(se)
}

test01 <- unlist(map2(0.1, c(1:1000), calc.se))
test03 <- unlist(map2(0.3, c(1:1000), calc.se))
test05 <- unlist(map2(0.5, c(1:1000), calc.se))

test.se <- bind_cols("n" = c(1:1000),"p = 0.1" = test01,"p = 0.3" = test03, 
                     "p = 0.5" = test05)
test.se.gather <- gather(test.se, key = "p", value = "se", -n)

ggplot(test.se.gather, aes(x = n, y = se, color = p))+
  geom_line()+
  xlab("N")+
  scale_x_continuous(breaks = c(0,200,400,600,800,1000))+
  ylab("Standard error")+
  scale_color_manual(values = c("black","blue","red"))+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19))
ggsave("runs/run_061224/stderror.png")




