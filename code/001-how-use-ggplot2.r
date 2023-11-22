library(tidyverse)
library(palmerpenguins) # This library contains collected data from penguin population

p <- ggplot(data = penguins)
print(p) # Genera una figura vacía
class(p) # [1] "gg"     "ggplot"

# Genera figura con ejes x e y pero vacía
p <-
  ggplot(data = penguins, aes(x = bill_length_mm, y = bill_depth_mm))
print(p)

# Usamos el símbolo + para añadir componentes
p <-
  ggplot(data = penguins, aes(x = bill_length_mm, y = bill_depth_mm)) +
  # Al introducir el objeto geométrico, genera el scatter plot
  geom_point()

print(p)

p <-
  ggplot(data = penguins, aes(x = bill_length_mm, y = bill_depth_mm))

# También podemos añadir componentes a una variable de clase ggplot
p <- p + geom_point()

# Usamos print(p) o simplemente p para visualizar el gráfico
print(p)
p


ggsave(
  file = "../assets/001_001_scatterplot.png",
  p,
  width = 4,
  height = 4,
  units = "in",
  dpi = 100
)
