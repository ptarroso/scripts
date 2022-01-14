# This code will generate random samples forming two groups that are 
# diametrically opposed in a 3 dimensional volume. It uses  PHYLIN to with a 
# custom 3 dimensinal Euclidiean distance to interpolate over a 3D grid.

library(phylin)
library(rgl)

set.seed(23)

# make an example 3 dimensions grid
grd <- expand.grid(x=seq(-1, 1, 0.1), y=seq(-1, 1, 0.1), z=seq(-1, 1, 0.1))

# just use a cross-section of the 3d grid
angle <- atan2(grd$y, grd$x)
mask <- angle > pi/4 | angle < -pi*3/4
grd <- grd[mask, ]

# Create some random species presence data with 2 lineages
n <- 25
lin1 <- data.frame(x=rnorm(n, 0.5, 0.5), y=rnorm(n, 0.5, 0.5), 
                   z=rnorm(n, 0.5, 0.5), lin=1)
lin2 <- data.frame(x=rnorm(n, -0.5, 0.5), y=rnorm(n, -0.5, 0.5), 
                   z=rnorm(n, -0.5, 0.5), lin=2)
lin <- rbind(lin1, lin2)

# User defined distance function
geo.3d.dist <- function(from, to) {
	# Calculate 3D Eucledian distances
    dst <- matrix(NA, nrow = nrow(from), ncol = nrow(to))
    dimnames(dst) <- list(rownames(from), rownames(to))
    for (i in 1:nrow(from)) {
		 dst[i, ] <- sqrt((to[, 1] - from[i, 1])^2 + 
                          (to[, 2] - from[i, 2])^2 + 
                          (to[, 3] - from[i, 3])^2)
	}
    return(dst)
}

# Geographical distance is calculated with 3D coordinates. 
# The following command is equivalent to "geo.dis <- dist(lin[,1:3])" but phylin 
# needs a distance function that accepts 'from' and 'to' arguments.
geo.dist <- geo.3d.dist(lin[,1:3], lin[,1:3]) 

# Genetic distances are derived linearly from geographical distances with an 
# added random error. To simplify the variogram, genetic distances higher than
# 2 units of geographical distance are kept constant
gen.dist <- geo.dist+rnorm(length(geo.dist))
gen.dist[geo.dist > 2] <- mean(gen.dist[geo.dist > 2])

# Build and fit model to semivariogram
gv <- gen.variogram(geo.dist, gen.dist)
gv <- gv.model(gv, range=3, sill=3.2)
plot(gv)

# interpolate lineage 1
l1 <- krig(lin$lin == 1, lin[,1:3], grd, gv, distFUN = geo.3d.dist, 
           neg.weights=FALSE)

l2 <- krig(lin$lin == 2, lin[,1:3], grd, gv, distFUN = geo.3d.dist, 
           neg.weights=FALSE)


# 3D interactive scatter plot with rgl
col.l1 <- hcl.colors(25, "RdYlGn")[as.integer(l1$Z*24)+1]

cubit <- cube3d(color="blue", alpha=0.3)
cubit$vb[cubit$vb == -1] <- 0
temp <- cubit
open3d()

for (i in 1:nrow(grd)) {
	temp$vb[1,] <- cubit$vb[1,]+grd$x[i]*10
	temp$vb[2,] <- cubit$vb[2,]+grd$y[i]*10
	temp$vb[3,] <- cubit$vb[3,]+grd$z[i]*10
	shade3d(temp, add=TRUE, color=col.l1[i], alpha=l1$Z[i])  
}


