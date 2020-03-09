library(ggplot2)
library(volesti)
P <- GenCube(3, 'H')
N <- 1000
val <- 1:N
for (i in 1:N) {
  val[i] <- print(8 - volume(P))
}
print(mean(val))