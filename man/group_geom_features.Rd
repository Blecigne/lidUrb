% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_geometric_features.R
\name{group_geom_features}
\alias{group_geom_features}
\title{Compute geometric features for a group on points}
\usage{
group_geom_features(las, group, feat_list)
}
\arguments{
\item{las}{a LAS file.}

\item{group}{character. The name of the grouping variable.}

\item{feat_list}{(optional) character. A vector containing the list of the geometric features to compute.
Can be: "Eigenvalue_sum","Omnivariance", "Eigenentropy","Anisotropy","Planarity",
"Linearity","PCA1","PCA2","Surface_variation","Sphericity","Curvature".}
}
\value{
the LAS with the geometric features.
}
\description{
Compute geometric features for a group on points
}
\examples{
\donttest{
# import data
file = system.file("extdata", "tree_no_leaves.las", package="lidUrb")
las = lidR::readLAS(file)

# cluster using DBSCAN to create a group
las@data[,cluster := dbscan::dbscan(las@data[,1:3],0.03,1)$cluster]

# compute Planarity and Linearity for each cluster
las = lidUrb::group_geom_features(las,group = "cluster",feat_list = c("Planarity","Linearity"))

# plot the result
lidR::plot(las,color="Linearity_G")
lidR::plot(las,color="Planarity_G")
}
}
