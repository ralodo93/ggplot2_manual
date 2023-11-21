library(tidyverse)
library(palmerpenguins)

p <- penguins %>% ggplot() + # Otra forma de introducir los datos
  geom_point(aes(x = bill_length_mm, y = bill_depth_mm)) + # Capa del objeto geométrico, indicando la capa estética
  scale_x_continuous(trans = "log10") + # capa escala para modificar la escala del eje x
  theme(axis.text = element_text(color = "blue")) + # capa estilo para modificar el color del texto de los ejes
  facet_wrap(~island) # capa multipanel

p

ggsave(file = "../assets/002_001_scatterplot.png",p,width = 4, height = 4, units = "in", dpi = 100)


p <- penguins %>% ggplot() + 
  geom_boxplot(aes(x = island, y = bill_depth_mm), outlier.shape = NA) + # usamos geom_boxplot para visualizar el boxplot
  geom_jitter(aes(x = island, y = bill_depth_mm), color = "darkred") # A continuación se introducen los puntos coloreados de rojo (el orden es importante)

p
ggsave(file = "../assets/002_002_boxpoints.png",p,width = 4, height = 4, units = "in", dpi = 100)

p <- penguins %>% ggplot() + 
  geom_boxplot(aes(x = island, y = bill_depth_mm), outlier.shape = NA) +
  geom_jitter(aes(x = island, y = bill_depth_mm, color = island)) # introducimos color en la función aes() y la asignamos a variable categórica

p

p1 <- penguins %>% ggplot() + 
  geom_jitter(aes(x = island, y = bill_depth_mm, color = body_mass_g)) # introducimos color en la función aes() y la asignamos a variable numérica

p1


library(cowplot)

p <- plot_grid(p,p1, nrow = 1)

ggsave(file = "../assets/002_003_multi.png",p,width = 8, height = 4, units = "in", dpi = 100)
