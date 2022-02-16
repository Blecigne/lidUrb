
#' Computes indices for leaf / wood segmentation
#'
#' @param las a LAS file.
#' @param k integer. The number of nearest neighbors to use in the geometric features computation
#' @param search_radius numeric. The searching distance to compute geometric features.
#'                      Note that \code{2*search_radius} is used to compute linearity
#'                      that which provides better estimates.
#' @param dbscan_eps numeric. Sets the esp parameter to pass to the
#'                   \code{\link[dbscan]{dbscan}} function.
#' @param min_cluster_size integer. The minimal size of a cluster required to compute the indices. Points located
#'                         in small clusters recieve a p_wood and SoD of their nearest well classified point.
#' @param reclass_trunk_th numeric. The minimal Threshold to reclass the lower clusters (i.e. the clusters with at least
#'                         one point with \code{Z <= reclass_trunk_th}) as wood, i.e. a \code{p_wood = SoD = 1}. A negative
#'                         value desable the reclassification.
#'
#' @return a LAS with two additionnal fields, \code{p_wood} and \code{SoD}, corresponding to a wood classification index. \code{p_wood}
#'         was introduced in Wang et al. 2019 and range from 0 to 1 with 1 being very high probability of wood. \code{SoD} was introduced
#'         in Wan et al. 2020 and range between -1 and 1, higher probability of wood class close to 1.
#'
#' @description This function implements the method presented in Wan et al. 2020 but replace the connected component algorithm by a
#'              dbscan clustering.
#'
#' @references
#' \itemize{
#'   \item Wang, D., Momo Takoudjou, S., & Casella, E. (2020). LeWoS: A universal leaf‐wood classification method to facilitate the 3D
#'         modelling of large tropical trees using terrestrial LiDAR. Methods in Ecology and Evolution, 11(3), 376-389.
#'   \item Wan, P., Zhang, W., Jin, S., Wang, T., Yang, S., & Yan, G. (2020). Plot-level wood-leaf separation of trees using terrestrial
#'         LiDAR data based on a segmentwise geometric feature classification method. Methods in Ecology and Evolution, 12, 2473–2486.
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # import the data
#' file = system.file("extdata", "tree_leaves.las", package="lidUrb")
#' las = lidR::readLAS(file)
#'
#' # clean and reduce point density
#' las = lidUrb::filter_noise(las,k = 6L, sigma = 0.8)
#' las = lidUrb::reduce_point_density(las,0.02)
#'
#' # compute two criterias for wood and leaves classification
#' las = LW_segmentation_dbscan(las)
#'
#' # plot the two wood indexes criterias
#' lidR::plot(las,color = "p_wood",legend = TRUE)
#' lidR::plot(las,color = "SoD",legend = TRUE)
#'
#' # assign a class base on p_wood
#' las@data[,wood := as.numeric(p_wood >= 0.9)]
#' # plot the hard classification
#' lidR::plot(las,color="wood",size=2,colorPalette = c("chartreuse4","cornsilk2"))
#'
#' # assign a class base on p_wood
#' las@data[,wood := as.numeric(SoD >= 0.99)]
#' # plot the hard classification
#' lidR::plot(las,color="wood",size=2,colorPalette = c("chartreuse4","cornsilk2"))
#' }

LW_segmentation_dbscan = function(las,k=10L,search_radius = 0.05, dbscan_eps = 0.03, min_cluster_size = 10L, reclass_trunk_th = 0.5){

  . = .GRP = .N = Linearity_G = N = Planarity_G = SoD = Sphericity_G =
  Surface_variation = X = Y = Z = cluster = label = p_wood = NULL

  if(class(las)[1] != "LAS") stop("las must be of type LAS")
  if(!is.integer(k) | k <= 0) stop("k must be an integer >= 1")
  if(!is.integer(min_cluster_size) | min_cluster_size <= 0) stop("min_cluster_size must be an integer >= 1")
  if(!is.numeric(search_radius) | search_radius <= 0) stop("search_radius must be numeric >= 0")
  if(!is.numeric(dbscan_eps) | dbscan_eps <= 0) stop("dbscan_eps must be numeric >= 0")
  if(!is.numeric(reclass_trunk_th)) stop("reclass_trunk_th must be numeric")

  # compute surface variation if does already exists
  if(is.null(las@data$Surface_variation)){
    las = TreeLS::fastPointMetrics(las,method = TreeLS::ptm.knn(k=k,r=search_radius),which_metrics = c("Curvature"))
    data.table::setnames(las@data,old="Curvature",new = "Surface_variation")
  }

  # label the point cloud into three categories based on Wan et al. 2020 criteria
  las@data[,label := 1]
  las@data[Surface_variation > 0.1, label := 2]
  las@data[Surface_variation > 0.2, label := 3]

  # points in third label are automatically classified as leaves
  las@data[label == 3, ':='(p_wood = 0, SoD = -1)]

  # cluster within labels 1 and 2
  las@data[label == 1,cluster := dbscan::dbscan(las@data[label == 1,1:3],dbscan_eps,minPts = 1)$cluster]
  las@data[label == 2,cluster := dbscan::dbscan(las@data[label == 2,1:3],dbscan_eps,minPts = 1)$cluster]
  las@data[label == 3,cluster := dbscan::dbscan(las@data[label == 3,1:3],dbscan_eps,minPts = 1)$cluster]
  las@data[,cluster := .GRP, by = .(cluster,label)] # global cluster index

  # cluster geometric features
  las = lidUrb::group_geom_features(las,group = "cluster", feat_list = c("Linearity","Planarity","Sphericity"))

  ##### compute the probability of wood classification as in Wang et al. 2019
  las@data[,p_wood := 0]
  for(i in seq(0.01,1,0.01)){
    las@data[Linearity_G >= i & label != 3,p_wood := p_wood+0.01]
  }

  ##### compute the Significance Of Difference index as in Wan et al. 2020
  las@data[label != 3,SoD := Linearity_G+(1-Linearity_G)*(Linearity_G-max(c(Planarity_G,Sphericity_G))),by = cluster]

  # reclass if reclass_tree_base is TRUE the tree base recieve a p_wood of 1 ans a SoD of 1
  if(reclass_trunk_th > 0){
    las@data[cluster %in% unique(las@data[Z<=min(Z)+reclass_trunk_th,cluster]),':='(p_wood = 1,SoD = 1)]
  }

  # class points in small clusters
  las@data[,N := .N, by = cluster]
  if(min(las@data$N) <= min_cluster_size){
    neib = FNN::knnx.index(data = las@data[N > min_cluster_size,.(X,Y,Z)],query = las@data[N <= min_cluster_size,.(X,Y,Z)],1)
    las@data[N <= min_cluster_size,':='(
      p_wood = las@data$p_wood[las@data$N > min_cluster_size][neib],
      SoD = las@data$SoD[las@data$N > min_cluster_size][neib])]
  }

  las@data[,':='(Surface_variation = NULL, label = NULL, cluster = NULL, Linearity_G = NULL,
                 Planarity_G = NULL, Sphericity_G = NULL, N = NULL)]

  return(las)
}


