library(tidyverse)
library(gtools)

test.props <- as.data.frame(rdirichlet(20, c(5,1,1,1,1,1,1,5)))
colnames(test.props) <- c("A","B","C","D","E","F","G","H")

prop.gather <- gather(test.props, key = "rep", value = "freq")
ggplot(data = prop.gather, aes(x = rep, y = freq))+
  geom_point()

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

compare.trials <- function(re.samples){
  compare <- map(re.samples, function(x) abs(x$freq - x$post_freq))
  summarize <- unlist(map(compare, mean))
  return(summarize)
}

test <- map(re.samples, compare.trials)
a <- bind_cols(test)
colnames(a) <- seq(5,100,5)
a.gather <- gather(a, value = "off", key = "sample_count")
a.gather$sample_count <- as.numeric(a.gather$sample_count)

ggplot(a.gather, aes(x= sample_count, y= off))+
  geom_point()+
  geom_smooth()


