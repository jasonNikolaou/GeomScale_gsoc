# ---- GeomScale (GSOC) ------
# --- Randomized LP solver ---
# --  Iasonas Nikolaou  ------

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
  P$A = rbind(A, c);
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
      return(x_new)
      break;
    }
  }
}
#=== example 1 =====
#expected output = (0, 0, 0)
#cost:
c = c(1,1,1)
#A, b
A = rbind(c(1,0,0), c(0,1,0), c(0,0,1),
           c(-1,0,0), c(0,-1,0), c(0,0,-1))
b = c(1,1,1,0,0,0)
#parameters
epsilon = 0.0001
precision = 0.000000001
walk = "RDHR"

print(RCP(A, b, c, walk, epsilon, precision))

#=== test 2 ====
#expected output = (0, 200)
#cost:
c = c(-1, 1)
#A, b
A = rbind(c(1,0), c(0,1), c(1,1),
          c(-1,0), c(0,-1))
b = c(200, 300, 400, 0, 0)
#parameters
epsilon = 0.00001
precision = 0.0000001
walk = "BW"

print(RCP(A, b, c, walk, epsilon, precision))
