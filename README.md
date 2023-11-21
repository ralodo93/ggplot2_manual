# Manual ggplot2

## Importancia de la visualización de datos

Por regla general es muy complicado analizar un tabla de datos a partir de los valores de la misma, sean numéricos o carácteres.  Por este motivo es muy importante conocer y saber aplicar las diferentes técnicas de visualización de datos. A continuación se muestra un ejemplo de un conjunto pequeño de datos y la visualización de los mismos.

Valores de una variable numérica distribuidos entre 30 muestras de dos condiciones diferentes en formato tabla:

| Condition | Values |  Condition | Values |
| -------- | ------- | -------- | ------- |
| A |  3.879049 | B | 10.573826 |
| A |  4.539645 | B |  7.995701 |
| A |  8.1174174 | B |  3.066766 |
| A |  5.141017 | B |  8.402712 |
| A |  5.258575 | B |  6.054417 |
| A |  8.430130 | B |  4.864353 |
| A |  5.921832 | B |  6.564050 |
| A |  2.469878 | B |  4.947991 |
| A |  3.626294 | B |  5.542218 |
| A |  4.108676 | B |  5.749921 |
| A |  7.448164 | B |  3.626613 |
| A |  5.719628 | B |  8.675574 |
| A |  5.801543 | B |  7.306746 |
| A |  5.221365 | B |  4.723726 |
| A |  3.888318 | B |  9.507630 |

Y mediante la visualización de dichos datos

![000_ex1](assets/000_01_example.png)

A simple vista, a partir de los datos de la tabla es muy complicado sacar ninguna conclusión de los mismos, sin embargo, mediante el uso de visualizaciones sencillas podemos interpretar de forma mucho más eficiente los datos de los que disponemos.

## ¿Cómo hacer un gráfico con ggplot2?

`ggplot2` es una librería de R que se engloba dentro de un conjunto de paquetes llamado `tidyverse`, que incluye algunas librerías muy útiles para la manipulación de datos como `dplyr` o `stringr`. Aunque en R hay numerosas formas de diseñar gráficos, `ggplot2` es quizá la más utilizada ya que es fácil de aprender y permite un elevado grado de personalización. Se basa en una sintaxis conocida como gramática de gráficos (grammar of graphics) mediante la cual, con solo conocer una serie de funciones podemos generar miles de gráficos diferentes. Además, `ggplot2` resulta sencillo ya que se trata de componer la figura mediante bloques. Por el contrario presenta una gran desventaja y es que el formato de entrada de los datos es muy estático, conocido como datos `tidy` y que consiste en tablas en las que cada registro se almacena en una fila de la tabla. Por lo general, los datos que se recogen suelen ser de este estilo salvo que nos encontremos con matrices, las cuales deben ser transformadas a formato `tidy`.

Antes de comenzar a generar gráficos, debemos tener en cuenta un aspecto específico de los gráficos de `ggplot2` y es que están construidos a base de capas. Para generar un gráfico es necesario utilizar al menos 3 capas:

- Datos. La tabla con los datos que vamos a visualizar, siempre en formato tidy.
- Objeto geométrico. La representación en la que vamos a visualizar los datos (puntos, boxplot, líneas, texto, etc)
- Estética: La forma en la que vamos a visualizar los datos (coordenadas x e y, etiquetas, colores etc.)

El resto de capas, que son las que nos permiten hacer de nuestro gráfico algo personal las iremos viendo a lo largo de este tutorial. Algunas de ellas son: escala, etiquetas, legenda, tema o paneles.

Ahora si, vamos a comenzar a generar nuestros primeros gráficos. Para ello usaremos 3 funciones fundamentales. La primera de ellas es `ggplot()` que se utiliza para indicar que se va a crear un objeto de R de tipo ggplot. Esta función se utiliza, generalmente para cargar los datos y definir la estética. Esta estética se determina con la función `aes()`. Finalmente usaremos una función para el objeto geométrico. Todas ellas siguen la misma gramática `geom_X()`, de este modo tenemos `geom_point()` para puntos, `geom_line()` para líneas o `geom_text()` para texto. A lo largo del tutorial veremos muchísimas.

```r
library(tidyverse)
library(palmerpenguins) # This library contains collected data from penguin population

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
```

![001_001_scatterplot](assets/001_001_scatterplot.png)

El código completo de este apartado está en [code/001-how-use-ggplot2.r](code/001-how-use-ggplot2.r)

## Personalizar un gráfico en ggplot2

En este apartado abordaremos las principales características de la visualización de datos con `ggplot2`. Para ello, se ha elaborado una serie de secciones en las que se trata de dar relevancia a cada una de ellas. Estas secciones son: cómo añadir y modificar capas, los principales objetos geométricos, escalas y colores, paneles múltiples y estilos.

### Añadir y modificar capas

Como hemos comentado, en `ggplot2` creamos los gráficos a partir de ir añadiendo capas. Cada capa o componente puede definir elementos muy diferentes: geometrías, estadística, escalas o estilos. Un ejemplo típico de comando de `ggplot2` sería este:

```r
library(tidyverse)
library(palmerpenguins)

p <- penguins %>% ggplot() + # Otra forma de introducir los datos
  geom_point(aes(x = bill_length_mm, y = bill_depth_mm)) + # Capa del objeto geométrico, indicando la capa estética
  scale_x_continuous(trans = "log10") + # capa escala para modificar la escala del eje x
  theme(axis.text = element_text(color = "blue")) + # capa estilo para modificar el color del texto de los ejes
  facet_wrap(~island) # capa multipanel

p
```

![002_001_scatterplot](assets/002_001_scatterplot.png)

También es posible concatenar varios objetos geométricos en una misma figura, por ejemplo, vamos a hacer un boxplot sobre el cual vamos a poner los puntos. No te preocupes si hay algo que no entiendas, en el siguiente apartado veremos las principales opciones de los objetos geométricos más usados.

```r
p <- penguins %>% ggplot() + 
  geom_boxplot(aes(x = island, y = bill_depth_mm), outlier.shape = NA) + # usamos geom_boxplot para visualizar el boxplot
  geom_jitter(aes(x = island, y = bill_depth_mm), color = "darkred") # A continuación se introducen los puntos coloreados de rojo (el orden es importante)

p
```

![002_002_boxplot](assets/002_002_boxpoints.png)

En la figura anterior hemos visto como colorear los puntos de un color determinado. Sin embargo, en ocasiones nos interesará colorearlos en base a alguna variable de nuestra tabla. Para ello debemos introducir el nombre de dicha variable asignándola al parámetro color, pero dentro de la estética, es decir, en la función `aes()`. La variable que se puede asignar a un color (o forma, tamaño) puede ser categórica (factor o carácter) o numérica. Vemos un ejemplo de cada una de ellas.

```r
p <- penguins %>% ggplot() + 
  geom_boxplot(aes(x = island, y = bill_depth_mm), outlier.shape = NA) +
  geom_jitter(aes(x = island, y = bill_depth_mm, color = island)) # introducimos color en la función aes() y la asignamos a variable categórica

p

p1 <- penguins %>% ggplot() + 
  geom_jitter(aes(x = island, y = bill_depth_mm, color = body_mass_g)) # introducimos color en la función aes() y la asignamos a variable numérica

p1
```

![002_003_multi](assets/002_003_multi.png)

El código completo de este apartado está en [code/002-add-modify-layers.r](code/002-add-modify-layers.r)


# capas

## stats
## position
## scale
## facet
## labels
## estilo



### Principales objetos geométricos

A continuación veremos los principales objetos geométricos que se utilizan a la hora de visualizar datos así como los parámetros que se usan y se pueden modificar en cada uno de ellos. Como hemos visto anteriormente, los objetos geométricos se llaman con la función `geom_X`. Si vemos la ayuda de cualquiera de estas funciones podremos hacernos una idea de los parámetros que se pueden modificar. Aquí tenéis el ejemplo de la función `geom_point`:

```r
geom_point(
  mapping = NULL,
  data = NULL,
  stat = "identity", # indica la transformación estadística, en este caso identity
  position = "identity", # indica el ajuste de posición, en este caso identity
  ...,
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE
)

# geom_point() tiene los siguientes "aesthetics"

# x coordenada x
# y coordenada y
# alpha transparencia
# colour/color color
# fill color de relleno
# group grupos
# shape forma
# size tamaño
# stroke tamaño borde
```

Para explicar los distintos objetos geométricos los vamos a clasificar según su uso:

| 1 Variable (continua) | 1 Variable (discreta) | 2 Variables Numéricas | 2 Variables (numérica y categórica) | 3 Variables |
| ---------------------- | ---------------------- | ---------------------- | ---------------------------------- | ------------ |
| geom_area              | geom_bar               | geom_point             | geom_boxplot                      | geom_tile    |
| geom_density           |               | geom_line              | geom_jitter                       |   |
| geom_histogram         |          | geom_label              | geom_violin                       |   |
| geom_dotplot       |           | geom_text            | geom_col                    |   |
|       |           |            |         geom_segment           |   |

Hay que tener en cuenta que las funciones se han agrupado según se suelen usar, lo cual no implica que algunas puedan usarse con otro tipo de variables distinto al referenciado en esta tabla. Antes de comentar cada uno de los objetos geométricos vamos a ver algunos ejemplos:

```r
library(tidyverse)
library(palmerpenguins)

# Variable continua
# eje x --> variable continua
# eje y --> Nada - geom_histogram calcula número de ocurrencias
p <- penguins %>% ggplot() +
  geom_histogram(aes(bill_length_mm))

# Variable discreta
# eje x --> variable discreta
# eje y --> Nada - geom_bar calcula número de ocurrencias
p <- penguins %>% ggplot() +
  geom_bar(aes(island))

# 2 variables numéricas
p <- penguins %>% ggplot() +
  geom_point(aes(bill_length_mm,bill_depth_mm))

# 2 variables: numérica + categórica
p <- penguins %>% ggplot() +
  geom_violin(aes(y = bill_length_mm,x = island))

# 3 variables
p <- penguins %>% group_by(year, island) %>% count() %>% ggplot() +
  geom_tile(aes(y = year, x = island, fill = n))
```

![003_001_multi](assets/003_001_examples.png)

#### Representar una variable continua

Generalmente para visualizar una variable continua nos fijamos en la distribución de los valores de dicha variable. Por ejemplo, si queremos visualizar la distribución de la longitud del pico (`bill_lenght_mm`) lo haremos principalmente mediante histogramas y/o perfiles de densidad. Otras alternativas son las áreas y los gráficos de puntos (`dotplot`).

##### geom_histogram

Un histograma es un gráfico de barras que muestra una distribución de frecuencias. Si vemos la ayuda de esta función:

```r
geom_histogram(
  mapping = NULL,
  data = NULL,
  stat = "bin",
  position = "stack",
  ...,
  binwidth = NULL,
  bins = NULL,
  na.rm = FALSE,
  orientation = NA,
  show.legend = NA,
  inherit.aes = TRUE
)

# geom_histogram() usa los siguientes "aesthetic"
# x (variable continua)
# y (eje y, si no se indica nada usa la función stat_bin() --> count)
# alpha (transparencia)
# colour/color (color de los bordes)
# fill (color de relleno)
# group (variable para agrupar)
# linetype (tipo de línea para borde)
# linewidth (tamaño del borde)

# Ejemplo

p <- penguins %>% ggplot() +
  geom_histogram(aes(bill_length_mm))
```

![003_002_hist](assets/003_002_example_hist.png)

Lo que se muestra en el histograma es el número de pingüinos que presentan un rango de longitud del pico. Ya que se trata de una variable continua, los rangos de longitud pueden variar y ser más pequeños o más grandes. Por defecto el número de barras es 30, tal y como se indica en el mensaje que se muestra cuando se lanza el comando. Podemos modificar o bien el número de barras o bien el tamaño de rango de las mismas con los parámetros `bins` y `binwidth` respectivamente.

```r
p <- penguins %>% ggplot() +
  geom_histogram(aes(bill_length_mm), bins = 10) # set number of bins

p <- penguins %>% ggplot() +
  geom_histogram(aes(bill_length_mm), binwidth = 5) # set width of bins
```

![003_003_bins](assets/003_003_bins.png)

A la izquierda se establece la distribución usando 10 *bins* y en la derecha se define un tamaño de 5 mm para cada *bin*.

