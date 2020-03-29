#Code from https://vissarion.github.io/tutorials/volesti_R_tutorial
#Sampling via random walks.
library(ggplot2)
library(volesti)
step <- 100
for (walk in c("BW")) {
  P <- GenCube(100, 'H')
  points1 <- sample_points(P, WalkType = walk, walk_step = step, N=1000)
  g<-plot(ggplot(data.frame( x=points1[1,], y=points1[2,] )) +
            geom_point( aes(x=x, y=y, color=walk)) + coord_fixed(xlim = c(-1,1),
                                                                 ylim = c(-1,1)) + ggtitle(sprintf("walk length=%s", step, walk)))
}
