library(tidyverse)
library(palmerpenguins)

p <- penguins %>% ggplot() + # Otra forma de introducir los datos
  
  # Capa del objeto geométrico, indicando la capa estética
  geom_point(aes(x = bill_length_mm, y = bill_depth_mm)) +
  # capa escala para modificar la escala del eje x
  scale_x_continuous(trans = "log10") +
  # capa estilo para modificar el color del texto de los ejes
  theme(axis.text = element_text(color = "blue")) +
  # capa multipanel
  facet_wrap( ~ island)

p

ggsave(
  file = "../assets/002_001_scatterplot.png",
  p,
  width = 4,
  height = 4,
  units = "in",
  dpi = 100
)


p <- penguins %>% ggplot() +
  # usamos geom_boxplot para visualizar el boxplot
  geom_boxplot(aes(x = island, y = bill_depth_mm), outlier.shape = NA) +
  # A continuación se introducen los puntos coloreados de rojo
  # El orden de los objetos geométricos es importante
  # En este caso los puntos se superponen al boxplot
  geom_jitter(aes(x = island, y = bill_depth_mm), color = "darkred")

p

ggsave(
  file = "../assets/002_002_boxpoints.png",
  p,
  width = 4,
  height = 4,
  units = "in",
  dpi = 100
)

p <- penguins %>% ggplot() +
  geom_boxplot(aes(x = island, y = bill_depth_mm), outlier.shape = NA) +
  # introducimos color en la función aes() y la asignamos a variable categórica
  geom_jitter(aes(x = island, y = bill_depth_mm, color = island))

p

p1 <- penguins %>% ggplot() +
  # introducimos color en la función aes() y la asignamos a variable numérica
  geom_jitter(aes(x = island, y = bill_depth_mm, color = body_mass_g))

p1


library(cowplot)

p <- plot_grid(p, p1, nrow = 1)

ggsave(
  file = "../assets/002_003_multi.png",
  p,
  width = 8,
  height = 4,
  units = "in",
  dpi = 100
)
