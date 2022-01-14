# This code generates a study area grid and 100 random samples. The samples are
# divided into two groups: Lineage 1 on the left and Lineage 2 on the right. A
# climate layer is used as connectivity. A second limate layer is derived by
# adding 4 degrees with some random noise, simulating a warmer climate. Phylin
# is used with resistance distances with both layers to illustrate how it could
# be applied under a climate warming scenario

library(raster)
library(phylin)
library(gdistance)

set.seed(33)

# read Annual Mean Temperature (AMT) raster
rst <- raster("AnnualMeanTemperature_CHELSA.tif")

# Create interpolation grid based on raster
grd <- coordinates(rst)

# Get 100 species locations by randomly selecting pixels with AMT < 10
sp.px <- sample(which(rst[] < 10), 100)
sp <- data.frame(coordinates(rst)[sp.px,], lin = 1)

# second lineage if x > 1.3 
sp$lin[sp$x > -1.3] <- 2

# plot created data
plot(rst)
points(sp$x, sp$y, pch=sp$lin)
legend("topright", pch=1:2, legend=c("Lineage 1", "Lineage 2"))

# Create a simple conductance surface to calculate distances 
# (check phylin vignettes and papers for more information)
# this creates a surface that is more pearmeable with lower temperatures
conductance <- rst
conductance[] <- 1/(1+exp((rst[]-12))) #Check with 'plot(rst[], conductance[])'

# The package 'gdistance' allows to calculate a resistance distance with climate
# It needs a transition matrix and a geocorretion (refer 'gdistance' vignetes 
# and publication to more help)
tr <- transition(conductance, mean, 8)
tr <- geoCorrection(tr, type="r")

# We build our custom distance function to use 'gdistance' transition matrix.
# Phylin needs the arguments 'from' and 'to' and we also provide 'tr' to pass
# the transition matrix to the distance calculation
res.dist <- function (from, to, tr) {
	nf <- nrow(from)
	allcoords <- as.matrix(rbind(from, to))
	dist <- as.matrix(commuteDistance(tr, allcoords))
	my.dist <- dist[1:nf, (nf+1):ncol(dist)]
	return(my.dist/10000) 
}

## now we can use the function to calculate distances between samples
r.dist <- res.dist(sp[,1:2], sp[,1:2], tr)

# Genetic distances are derived linearly from climate resistance  distances with
# a small random value. To simplify the variogram, genetic distances higher than
# 4 units of resistance distance are kept constant
gen.dist <- r.dist
gen.dist[r.dist > 4] <- min(gen.dist[r.dist > 4])
gen.dist <- gen.dist + rnorm(length(r.dist), sd=0.1)


# Build and fit model to semivariogram
gv <- gen.variogram(r.dist, gen.dist, lag=0.2)
gv <- gv.model(gv, range=5, sill=8)
plot(gv)

# interpolate lineages (note: 'gdistance' is having trouble calculating the last
# pixel, so we avoid it, giving NA to this corner). Create a raster with lineage
# area defined as prob greater than 0.75
lin.rst <- rst * 0

for (l in 1:2) {
	lin <- krig(sp$lin == l, sp[,1:2], grd[1:8711,], gv, distFUN = res.dist, 
			    tr = tr, neg.weights=FALSE)

	lin.rst[1:8711] <- lin.rst[1:8711] + (lin$Z > 0.9) * l
}

plot(lin.rst)


# Simulate warming climate by adding a value around 4 to each pixel
rst.warm <- rst
rst.warm[] <- rst[] + rnorm(length(rst), 4, 0.1)

cond.warm <- rst
cond.warm[] <- 1/(1+exp((rst.warm[]-12))) #Check with 'plot(rst[], conductance[])'

# Recalculate conductance with warmer climate and use it interpolation
tr.warm <- transition(cond.warm, mean, 8)
tr.warm <- geoCorrection(tr.warm, type="r")

lin.warm.rst <- rst * 0

for (l in 1:2) {
	lin <- krig(sp$lin == l, sp[,1:2], grd[1:8711,], gv, distFUN = res.dist, 
			    tr = tr.warm, neg.weights=FALSE)

	lin.warm.rst[1:8711] <- lin.warm.rst[1:8711] + (lin$Z > 0.9) * l
}
plot(lin.warm.rst)

