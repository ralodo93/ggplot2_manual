library(tidyverse)
# This library contains collected data from penguin population
library(palmerpenguins)

p <- ggplot(data = penguins)
print(p) # Genera una figura vacía
class(p) # [1] "gg"     "ggplot"

p <- ggplot(data = penguins, aes(x = bill_length_mm, y = bill_depth_mm)) # Genera figura con ejes x e y pero vacía
print(p)

p <- ggplot(data = penguins, aes(x = bill_length_mm, y = bill_depth_mm)) + # Usamos el símbolo + para añadir componentes
  geom_point() # Al introducir el objeto geométrico, genera el scatter plot
print(p)

p <- ggplot(data = penguins, aes(x = bill_length_mm, y = bill_depth_mm))
p <- p + geom_point() # También podemos añadir componentes a una variable de clase ggplot

# Usamos print(p) o simplemente p para visualizar el gráfico
print(p)
p


ggsave(file = "../assets/001_001_scatterplot.png",p,width = 4, height = 4, units = "in", dpi = 100)
