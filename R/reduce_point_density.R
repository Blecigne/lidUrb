
#' Reduce the point density by keeping the point the closest to the center of a
#' voxel of given resolution
#'
#' @param las a LAS class object.
#' @param res numeric. The voxel resolution.
#' @param keep_points logical. Allow to keep all the points in the point cloud
#'                    and label the reduced density points as \code{is.low.density == 1}
#'                    as well as the grpou of points each point belongs to (allow easy matching
#'                    after running a function on the downsampled point cloud).
#'
#' @return the LAS file with lower reduced density.
#'
#' @export
#'
#' @examples
#' \donttest{
#' file = system.file("extdata", "urban.las", package="lidUrb")
#' las = lidR::readLAS(file)
#'
#' # reduce point density with large voxels
#' las = lidUrb::reduce_point_density(las,0.2)
#'
#' # plot
#' lidR::plot(las)
#' }

reduce_point_density = function(las,res, keep_points = FALSE){

  # to pass check
  X = Y = Z = D = Xvox = Yvox = Zvox = is.min = . = pts.grp = .GRP =
  is.low.density = NULL

  if(class(las)[1] != "LAS") stop("las must be of type LAS")
  if(missing(res)) stop("res must be specified")
  if(!is.numeric(res)) stop("res must be numeric")
  if(class(las)[1] != "LAS") stop("las must be of type LAS")

  # voxel
  las@data[,':='(
    Xvox = round(X/res)*res,
    Yvox = round(Y/res)*res,
    Zvox = round(Z/res)*res
  )]

  # distance to voxel center
  las@data[,D := sqrt( (X-Xvox)^2 + (Y-Yvox)^2 + (Z-Zvox)^2)]

  # find the point the closest from the voxel center
  las@data[,is.low.density := as.numeric(D == min(D)), by= .(Xvox,Yvox,Zvox)]

  if(keep_points){
    las@data[,pts.grp := .GRP, by =.(Xvox,Yvox,Zvox)]
    las@data[,':='(Xvox = NULL, Yvox = NULL, Zvox = NULL, D = NULL)]
  }else{
    # keep the nearest point
    las@data = las@data[las@data$is.low.density == 1]
    las@data[,':='(Xvox = NULL, Yvox = NULL, Zvox = NULL, D = NULL, is.low.density = NULL)]
  }

  return(las)
}
