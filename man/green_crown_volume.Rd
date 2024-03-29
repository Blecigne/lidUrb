% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/green_crown_volume.R
\name{green_crown_volume}
\alias{green_crown_volume}
\title{Compute the green crown volume, area and return a mesh}
\usage{
green_crown_volume(las, npts_in_clust = 200)
}
\arguments{
\item{las}{a LAS file with a classification field for the wood class.}

\item{npts_in_clust}{numeric. The average number of points in a cluster.
Defines the size of each cluster and therefore the size
of the convex hulls.}
}
\value{
a mesh object of the green crown and the green crown volume and area
}
\description{
Compute the green crown volume, area and return a mesh
}
\examples{
\donttest{
# import data
file = system.file("extdata", "tree_leaves.las", package="lidUrb")
las = lidR::readLAS(file,select = "xyz")

# filter noise
original = las
original@data[,original_index := 1:nrow(las@data)]

# filter noise and reduce point density
las = lidUrb::reduce_point_density(original,0.02)
las = lidUrb::filter_noise(las)

# segment leaves
las = lidUrb::LW_segmentation_dbscan(las)
las@data[,wood := as.numeric(p_wood >= 0.95)] # binary class (0 for leaves, 1 for wood)

# compute green crown volume
GCV = lidUrb::green_crown_volume(las)

# plot the green crown mesh over the original LAS
lidR::plot(las,color="wood",size=2,colorPalette = c("chartreuse4","cornsilk2"),clear_artifacts = FALSE)
rgl::shade3d(GCV$mesh,col = "chartreuse4",add=T)

# green crown volume and area
GCV$Green_crow_volume
GCV$Green_crown_area
}
}
\references{
Zhu, Z., Kleinn, C., & Nölke, N. (2020). Towards tree green crown
            volume: a methodological approach using terrestrial laser scanning.
            Remote Sensing, 12(11), 1841.
}
