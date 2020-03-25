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

#calculate number of samples
numOfSamples = function(k, epsilon) {
  return (ceiling(2.2 * log(1/epsilon) + 1.1 + 0.0505*log(k)))
}

# find min point
find_min = function(points, c) {
  min_dot = dot(c, points[,1])
  min = points[,1]
  for (i in 1:ncol(points)) {
    tmp = dot(c, points[,i])
    if (tmp < min_dot) {
      min = points[,i]
      min_dot = tmp
    }
  }
  return(min)
}

#randomized cutting plane algorithm
# A, b: Polytope in H-representation, i.e. Ax <= b
# c: cost vector
# walk: WalkType for sampling
# epsilon: probability used for convergence
# precision: required precision of the answer
RCP = function(A, b, c, walk="CDHR", epsilon=0.001, precision=0.0001) {
  P = Hpolytope$new(A, b)
  k = 1;
  N_k = numOfSamples(k, epsilon)
  points = sample_points(P, WalkType=walk, N=N_k)
  P$A = rbind(A, c); #comment
  x = 0
  x_new = 0
  repeat {
    x = x_new
    x_new = find_min(points, c)
    P$b = c(b, dot(c, x_new)) #update polytope
    k = k + 1
    N_k = numOfSamples(k, epsilon)
    points = sample_points(P, WalkType=walk, N=N_k)
    
    if (dist(x, x_new) < precision) {
      return (k)
      #return (x_new)
      break;
    }
  }
}

#benchmarks
i = 0;

res = 3:99
time = 3:99
for (d in 3: 100) {
  t1 = Sys.time()
  P = GenCube(d, 'H')
  A = P$A
  b = P$b
  c = runif(d, -1, 1)
  res[i] = RCP(A, b, c)
  t2 = Sys.time()
  time[i] = (t2 - t1)/60
  i = i + 1;
}
time
plot(time, type="o", ann=FALSE)
title(main="time - dimensions")
title(xlab="dimensions")
title(ylab="time")
legend("topleft", legend=c("ε = 10E-3", "α = 10E-4"))



