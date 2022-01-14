# The following code generates two species distributed longitudinally with a
# contact zone where both species can be found. This code illustrates the use
# of a Jaccard distance matrix (instead of a genetic distance) to interpolate
# the distribution of the species (without modelling environmental niche) and to
# map a contact area between them.

library(raster)
library(phylin)
library(philentropy)

set.seed(12)

# Create two species distributed on the left and right side of the grid
sp1 <- data.frame(x=rnorm(50, -0.5, 0.5), y=rnorm(50, 0, 0.5), sp="sp1")
sp2 <- data.frame(x=rnorm(50, 0.5, 0.5), y=rnorm(50, 0, 0.5), sp="sp2")
sp <- rbind(sp1, sp2) 

# create a raster defining study area extent (each pixel with a different int)
res <- 0.1
rst <- raster(xmn=min(sp$x)-res, xmx=max(sp$x)+res, 
              ymn=min(sp$y)-res, ymx=max(sp$y)+res, resolution=res)
rst[] <- 1:length(rst)

# Plot everything to check distributions
plot(rst)
points(sp$x, sp$y, pch=ifelse(sp$sp == "sp1", 1, 2))
legend("topright", pch=1:2, legend=c("sp1", "sp2"))

# calculate Jaccard distance between sampled locations (use pixel ID)
sp$sites <- extract(rst, sp[,1:2])
occ.sp <- xtabs( ~ sites + sp, sp) > 0
jac.dist <- distance(occ.sp, method="jaccard")
dimnames(jac.dist) <- list(rownames(occ.sp), rownames(occ.sp))

# geographical distance (from sites, not sample locations)
# NOTE: using indices as dim.names is not the best coding practice...
sites.crd <- coordinates(rst)[as.integer(rownames(jac.dist)),]
rownames(sites.crd) <- rownames(jac.dist)
geo.dist <- as.matrix(dist(sites.crd))

# Build and fit model to semivariogram
gv <- gen.variogram(geo.dist, jac.dist, lag=0.25)
gv <- gv.model(gv, range=3, sill=0.5)
plot(gv)

# interpolate lineages (note: 'gdistance' is having trouble calculating the last
# pixel, so we avoid it, giving NA to this corner). Create a raster with lineage
# area defined as prob greater than 0.75
grd <- coordinates(rst)

i1 <- krig(sp$sp == "sp1", sp[,1:2], grd, gv, neg.weights=FALSE)
i2 <- krig(sp$sp == "sp2", sp[,1:2], grd, gv, neg.weights=FALSE)

rst1 <- rst2 <- rst
rst1[] <- i1$Z
rst2[] <- i2$Z

# plot simpatric area
plot(rst1*rst2)
points(sp$x, sp$y, pch=ifelse(sp$sp == "sp1", 1, 2))
legend("topright", pch=1:2, legend=c("sp1", "sp2"))

