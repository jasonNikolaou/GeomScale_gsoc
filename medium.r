rm(list=ls())
library(volesti)

#dot product
dot <- function(a, b) {
  res = a * b
  return(sum(res))
}

#distance
dist <- function(x1, x2) {
  return(sqrt(sum((x1 - x2) ^ 2)))
}

#find min point
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
# P: Polytope in H-representation
# walk: WalkType for sampling
# N: number of samples
# c: cost vector
# err: calculate an err-optimal solution

RCP = function(A, b, walk, N, c, err) {
  P = Hpolytope$new(A, b)
  points = sample_points(P, WalkType=walk, N=d)
  A = rbind(A, c)
  b = c(b, 0)
  NROWS = length(b)
  x = 0
  x_new = 0
  repeat {
    x = x_new
    x_new = find_min(points, c)
    b[NROWS] = dot(c, x_new)
    P = Hpolytope$new(A, b)
    points = sample_points(P, WalkType=walk, N=d)
    if (dist(x, x_new) < err) {
      return(x_new)
      break;
    }
  }
}
#=== example 1 =====
#expected output = (0, 0, 0)

#dimensions:
d = 3
#cost:
c = c(1,1,1)
# Ax <= b
A = rbind(c(1,0,0), c(0,1,0), c(0,0,1),
           c(-1,0,0), c(0,-1,0), c(0,0,-1))
b = c(1,1,1,0,0,0)
#error
err = 0.000000001
#unit cube
walk = "RDHR"
RCP(A, b, walk, d, c, err)

#=== example 2 ====
#expected output = (0, 200)
d = 2
err = 0.000000001
c = c(-1, 1)
A = rbind(c(1,0), c(0,1), c(1,1),
          c(-1,0), c(0,-1))
b = c(200, 300, 400, 0, 0)
walk = "BW"
RCP(A, b, walk, d, c, err)
