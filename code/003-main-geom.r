library(tidyverse)
library(palmerpenguins)
library(cowplot)

# Variable continua
# eje x --> variable continua
# eje y --> Nada - geom_histogram calcula número de ocurrencias
p1 <- penguins %>% ggplot() +
  geom_histogram(aes(bill_length_mm))

# Variable discreta
# eje x --> variable discreta
# eje y --> Nada - geom_bar calcula número de ocurrencias
p2 <- penguins %>% ggplot() +
  geom_bar(aes(island))

# 2 variables numéricas
p3 <- penguins %>% ggplot() +
  geom_point(aes(bill_length_mm,bill_depth_mm))

# 2 variables: numérica + categórica
p4 <- penguins %>% ggplot() +
  geom_violin(aes(y = bill_length_mm,x = island))

# 3 variables
p5 <- penguins %>% group_by(year, island) %>% count() %>% ggplot() +
  geom_tile(aes(y = year, x = island, fill = n))

p <- plot_grid(p1,p2,p3,p4,p5, nrow = 3)

ggsave(file = "../assets/003_001_examples.png",p,width = 6, height = 4, units = "in", dpi = 100)
ggsave(file = "../assets/003_002_example_hist.png",p1,width = 4, height = 4, units = "in", dpi = 100)

p1 <- penguins %>% ggplot() +
  geom_histogram(aes(bill_length_mm), bins = 10) # set number of bins

p2 <- penguins %>% ggplot() +
  geom_histogram(aes(bill_length_mm), binwidth = 5) # set width of bins

p <- plot_grid(p1,p2, nrow = 1)
ggsave(file = "../assets/003_003_bins.png",p,width = 6, height = 4, units = "in", dpi = 100)


penguins %>% ggplot()+
  geom_point(aes(bill_length_mm, island))
penguins %>% ggplot() +
  geom_histogram(aes(bill_length_mm, fill = island), alpha = 0.5, color = "gray40")+
  geom_vline(xintercept = 45)+
  geom_vline(xintercept = 47)


