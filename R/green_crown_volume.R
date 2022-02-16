#' Compute the green crown volume, area and return a mesh
#'
#' @param las a LAS file with a classification field for the wood class.
#' @param npts_in_clust numeric. The average number of points in a cluster.
#'                      Defines the size of each cluster and therefore the size
#'                      of the convex hulls.
#'
#' @return a mesh object of the green crown and the green crown volume and area
#'
#' @references Zhu, Z., Kleinn, C., & Nölke, N. (2020). Towards tree green crown
#'             volume: a methodological approach using terrestrial laser scanning.
#'             Remote Sensing, 12(11), 1841.
#' @export
#'
#' @examples
#' \donttest{
#' # import data
#' file = system.file("extdata", "tree_leaves.las", package="lidUrb")
#' las = lidR::readLAS(file,select = "xyz")
#'
#' # filter noise
#' original = las
#' original@data[,original_index := 1:nrow(las@data)]
#'
#' # filter noise and reduce point density
#' las = lidUrb::reduce_point_density(original,0.02)
#' las = lidUrb::filter_noise(las)
#'
#' # segment leaves
#' las = lidUrb::LW_segmentation_dbscan(las)
#' las@data[,wood := as.numeric(p_wood >= 0.95)] # binary class (0 for leaves, 1 for wood)
#'
#' # compute green crown volume
#' GCV = lidUrb::green_crown_volume(las)
#'
#' # plot the green crown mesh over the original LAS
#' lidR::plot(las,color="wood",size=2,colorPalette = c("chartreuse4","cornsilk2"),clear_artifacts = FALSE)
#' rgl::shade3d(GCV$mesh,col = "chartreuse4",add=T)
#'
#' # green crown volume and area
#' GCV$Green_crow_volume
#' GCV$Green_crown_area
#' }
green_crown_volume = function(las,npts_in_clust = 200){

  . = .GRP = .N = N = X = Y = Z = axis_ID = clusters = index = is.low.density =
  pts.grp = volume = wood = NULL

  if(class(las)[1] != "LAS") stop("las must be of type LAS")
  if(!is.numeric(npts_in_clust)) stop("npts_in_clust must be numeric")
  if(is.null(las@data$wood)){
    stop("las must include a leaf/wood segmentation field named 'wood' which is
         binary: 0 stands for leaves points and 1 for wood points")
  }

  leaves = las@data[wood == 0]

  # clustering the leaves point cliud using fastKmeans and defines the number of clusters
  # as the number of points devided by npts_in_clust
  leaves[,clusters := Rvcg::vcgKmeans(as.matrix(leaves[,1:3]),k=round(nrow(leaves)/npts_in_clust))$class]

  # number of points per cluster
  leaves[,N := .N, by=clusters]

  # keep only thoose with more than 4 points (needed to fit a convex hull)
  leaves = leaves[N >= 4]

  # add an index to the point cloud
  leaves[,index := 1:nrow(leaves)]

  hull_out = matrix(ncol=3) # matrix to store the hull
  volume = 0 # total green crown volume
  area = 0 # total green crown area
  for(i in unique(leaves$clusters)){
    # compute convex hull for the cluster
    hull = geometry::convhulln(leaves[clusters == i,.(X,Y,Z)],output.options = TRUE)

    # compute total volume and area
    volume = volume + hull$vol
    area = area + hull$area

    # store hull for mesh
    hull = matrix(c(leaves$index[leaves$clusters == i][hull$hull]),ncol=3)
    hull_out = rbind(hull_out,hull)
  }
  mesh = rgl::mesh3d(vertices = as.matrix(t(leaves[,1:3])),triangles = as.matrix(t(hull_out)),meshColor = "faces")

  return(list(mesh = mesh, Green_crow_volume = volume, Green_crown_area = area))
}
