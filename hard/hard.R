# ---- GeomScale (GSOC) ------
# --- Randomized LP solver ---
# --  Iasonas Nikolaou  ------
rm(list=ls())
library(volesti)

#dot product
dot <- function(a, b) {
  res = a * b
  return(sum(res))
}

#distance
dist = function(x1, x2) {
  return(sqrt(sum((x1 - x2) ^ 2)))
}

simAnnealing = function(A, b, c) {
  P = Hpolytope$new(A, b) # construct polytope
  n = ncol(A) # dimensions
  x = sample_points(P, N=1) # initial point
  T = 1000 # temperature
  temperature_factor = (1 - 1/sqrt(n)) #factor to decrease variance.
  V = 1000 # variance
  variance_factor = 0.9 #factor to decrease variance.
  while(T > 1e-10) {
    E_c = dot(c, x)
    x_next = c(sample_points(P, N=1, distribution="gaussian",
                             Parameters=list("variance" = V), InnerPoint = x))
    E_n = dot(c, x_next)
    DE = E_c - E_n # energy difference
    if (DE > 0) {
      #x_next is a better solution
      x = x_next
    }
    else if (exp(DE/T) > runif(1, 0, 1)) {
      #x_next is not a better solution
      #but it is selected
      x = x_next
    }
    #decrease temperature
    T = T * temperature_factor
    #decrease variance
    V = V * variance_factor
  }
  return (x)
}