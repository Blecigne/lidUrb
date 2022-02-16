
#' Compute Leaf Area Index, Total Tree Area and Leaf Area Density profiles from TLS point cloud
#'
#' @param leaves_pc a data.table containing the XYZ coordinates of a tree leaves
#' @param layer_thickness numeric (optional). The layer tickness to compute LAD profiles.
#'                       If not provided the voxel size if use as layer tickness.
#' @description This function implements the "Voxel method" first described by Hosoi and Omasa (2006) to estimate
#'              tree LAI, Total Leaf Area and LAD profile. The voxel resolution is automatically computed from
#'              the point cloud as the average shortest distance of a point to its nearest neighbor as suggested
#'              in Li et al. 2016.
#'
#' @references
#' \itemize{
#'   \item Hosoi, F., & Omasa, K. (2006). Voxel-based 3-D modeling of individual trees for estimating leaf area density using
#'         high-resolution portable scanning lidar. IEEE transactions on geoscience and remote sensing, 44(12), 3610-3618.
#'   \item Li, Y., Guo, Q., Tao, S., Zheng, G., Zhao, K., Xue, B., & Su, Y. (2016). Derivation, validation, and sensitivity
#'         analysis of terrestrial laser scanning-based leaf area index. Canadian Journal of Remote Sensing, 42(6), 719-729.
#' }
#'
#' @return A list containing the LAI, TLA as numeric values and the LAD profile returned as a data.table.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # import the data
#' file = system.file("extdata", "tree_leaves.las", package="lidUrb")
#' las = lidR::readLAS(file)
#'
#' # filter noise
#' las = lidUrb::filter_noise(las,k = 6L,sigma = 0.8)
#'
#' # segment foliage
#' las = lidUrb::LW_segmentation_dbscan(las)
#' las@data[,wood := as.numeric(p_wood >= 0.9)]
#'
#' # compute leaves traits
#' leaves = lidUrb::leaves_traits(las@data[wood == 0])
#'
#' # estimated LAI
#' leaves$LAI
#' # estimated Total Leaf Area
#' leaves$TLA
#'
#' # plot the LAD profile
#' library(ggplot2)
#' ggplot(leaves$LAD_profile, aes(x = layer, y= LAD)) +
#'   geom_line(size=1.1) +
#'   coord_flip() +
#'   theme(legend.position = "none") +
#'   xlab("Layer elevation") +
#'   ylab("Leaf Area Density")
#' }

leaves_traits = function(leaves_pc,layer_thickness = 0.2){

  . = .N = Contact_frequency = LAD = N = X = Xvox = Y = Yvox = Z = Zvox = layer = NULL

  if(!"data.table" %in% class(leaves_pc)) stop("leaves_pc must de a data.table object.")
  if(ncol(leaves_pc) < 3) stop("leaves_pc must have at least three columns.")
  if(!all(names(leaves_pc)[1:3] == c("X","Y","Z"))) stop("leaves_pc coordinates must be named X Y Z")

  min_Z = min(leaves_pc$Z)

  # normalize leaves position
  leaves_pc[,':='(
    X = X - min(X),
    Y = Y - min(Y),
    Z = Z - min_Z
  )]

  # compute voxel size based of mean distance to nearest neihbor
  voxel_size = mean(FNN::knn.dist(leaves_pc[,.(X,Y,Z)],1))

  # voxelize the leaves point cloud
  leaves_pc[,':='(
    Xvox = round(X/voxel_size)*voxel_size,
    Yvox = round(Y/voxel_size)*voxel_size,
    Zvox = round(Z/voxel_size)*voxel_size
  )]

  # vox
  vox = unique(leaves_pc[,.(Xvox,Yvox,Zvox)])

  # compute the potential number of voxels based on a convexhull area
  N_pot_voxel  = geometry::convhulln(unique(leaves_pc[,1:2]),output.options = "FA")$area/(voxel_size^2)

  #create layers
  if(missing(layer_thickness)){
    vox[,layer := Zvox]
    layer_thickness = voxel_size
  }else{
    vox[,layer := round(Zvox/layer_thickness)*layer_thickness]
  }

  # compute number of voxels in each layer
  vox[,N:=.N, by = Zvox]

  # compute contact frequency in each layer
  vox[,Contact_frequency := N/N_pot_voxel]

  # compute tree LAI
  LAI = 1.1*sum(unique(vox[,.(Zvox,Contact_frequency)])$Contact_frequency)

  # compute tree total leaf area
  TLA = LAI*N_pot_voxel*voxel_size^2

  # compute LAD in each layer
  Lad_tab = unique(vox[,.(layer,Zvox,Contact_frequency)])
  Lad_tab[,LAD := 1.1*(1/layer_thickness)*sum(Contact_frequency),by = layer]

  return(list(
    LAI = LAI,
    TLA = TLA,
    LAD_profile = unique(Lad_tab[order(layer),.(layer,LAD)])
  ))
}

