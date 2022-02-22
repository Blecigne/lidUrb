#' Computes local geometrics features using a constant research distance
#'
#' @description Compute local geometric features of the point cloud based on the eigenvalues of
#'              the points covariance matrix. This function uses the python implementation provided
#'              in \href{https://jakteristics.readthedocs.io/en/latest/}{jakteristics}.
#'
#' @param las a LAS file.
#' @param search_radius numeric. The search distance to retrieve point neighborhood.
#' @param features_list (optional) character. A vector containing the list of the geometric features to compute.
#'                      Can be: "Eigenvalue_sum","Omnivariance", "Eigenentropy","Anisotropy","Planarity",
#'                      "Linearity","PCA1","PCA2","Surface_variation","Sphericity","Verticality",
#'                      "Nx","Ny","Nz".
#'
#' @return The LAS file with the geometric features added as new columns in the slot data. NOTE: the names
#'         of the features is a contraction of the feature name and search_radius.
#' @export
#'
#' @examples
#' \donttest{
#' # import data
#' file = system.file("extdata", "urban.las", package="lidUrb")
#' las = lidR::readLAS(file)
#'
#' # compute all features
#' las = lidUrb::geom_features(las,search_radius = 0.1)
#'
#' # plot Planarity and Linearity
#' lidR::plot(las,color="Planarity")
#' lidR::plot(las,color="Linearity")
#' }

geom_features = function(las, search_radius, features_list){

  # test the parameters
  if(class(las)[1] != "LAS") stop("las must be of type LAS")
  if(missing(search_radius)) stop("res must be specified")
  if(!is.numeric(search_radius)) stop("res must be numeric")
  if(!missing(features_list)){
    if(!is.character(features_list)) stop("features_list must be character")
  }

  if(!reticulate::py_module_available("jakteristics")){
    stop(
      "The geom_features function relies on the jakteristics python package
      to compute the geometric features of the point cloud. Please refer to the
      https://github.com/Blecigne/lidUrb to find some help about how to install
      it."
    )
  }

  jak = reticulate::import("jakteristics")
  features = data.table::data.table(jak$compute_features(as.matrix(las@data[,1:3]),search_radius))

  # name the geometric features
  data.table::setnames(features,c(
    "Eigenvalue_sum","Omnivariance", "Eigenentropy","Anisotropy","Planarity",
    "Linearity","PCA1","PCA2","Surface_variation","Sphericity","Verticality",
    "Nx","Ny","Nz"))

  # return the data with the needed geometric features
  if(missing(features_list)){
    las@data = cbind(las@data,features)
  }else{
    las@data = cbind(las@data,features[,features_list, with=FALSE])
  }
  return(las)
}

