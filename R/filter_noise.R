
#' Filter noise from TLS scene
#'
#' @param las a LAS class object.
#' @param k numeric. the number of nearest neighbors to compute distance.
#' @param sigma numeric. The multiplier of the distance standard deviation to
#'              class points as noise.
#' @param keep_noise logical. If TRUE the points classified as noise are kept in
#'                   the point cloud and are labeled in the is.noise fields,
#'                   if FALSE points classified as noise are removed.
#'
#' @return the LAS with a new field specifying if the points are noise (if keep_nois is TRUE) or
#'         the LAS without noise points (if keep_noise is FALSE).
#'
#' @export
#'
#' @examples
#' \donttest{
#' # import data
#' file = system.file("extdata", "urban.las", package="lidUrb")
#' las = lidR::readLAS(file)
#'
#' # filter noise and store it for ploting
#' las = lidUrb::filter_noise(las,keep_noise = TRUE)
#'
#' # plot noise in red
#' lidR::plot(las,color="is.noise")
#' }

filter_noise = function(las,k = 6L,sigma = 2,keep_noise = FALSE){

  # to pass check
  D = . = X = Y = Z = is.noise = sd = NULL

  if(class(las)[1] != "LAS") stop("las must be of type LAS")
  if(!is.integer(k) | k <= 0) stop("k must be an integer >= 1")
  if(!is.numeric(sigma) | sigma <= 0) stop("sigma must be numeric > 0")
  if(!is.logical(keep_noise)) stop("keep_noise must be logical")

  # compute points average distance to their k nearest neighbors
  las@data[,D := rowMeans(FNN::knn.dist(data = as.matrix(las@data[,.(X,Y,Z)]), k = k))]

  # points are classified as noise if their distance is greater than
  # sigma*sdandard deiavtion of distance
  las@data[,is.noise := D > mean(D)+sigma*sd(D)]

  if(keep_noise){
    # if keep_noise is TRUE keep the is.noise field
    las@data[,D := NULL]
    return(las)
  }else{
    # if keep_noise is FALSE remove the noise from the point cloud as well as
    # the is.noise field
    las@data = las@data[!las@data$is.noise]
    las@data[,':='(D = NULL, is.noise = NULL)]
    return(las)
  }
}
