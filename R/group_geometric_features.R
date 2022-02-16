#' Compute geometric features for a group on points
#'
#' @param las a LAS file.
#' @param group character. The name of the grouping variable.
#' @param feat_list (optional) character. A vector containing the list of the geometric features to compute.
#'                      Can be: "Eigenvalue_sum","Omnivariance", "Eigenentropy","Anisotropy","Planarity",
#'                      "Linearity","PCA1","PCA2","Surface_variation","Sphericity","Curvature".
#'
#' @return the LAS with the geometric features.
#' @export
#'
#' @examples
#' \donttest{
#' # import data
#' file = system.file("extdata", "tree_no_leaves.las", package="lidUrb")
#' las = lidR::readLAS(file)
#'
#' # cluster using DBSCAN to create a group
#' las@data[,cluster := dbscan::dbscan(las@data[,1:3],0.03,1)$cluster]
#'
#' # compute Planarity and Linearity for each cluster
#' las = lidUrb::group_geom_features(las,group = "cluster",feat_list = c("Planarity","Linearity"))
#'
#' # plot the result
#' lidR::plot(las,color="Linearity_G")
#' lidR::plot(las,color="Planarity_G")
#' }

group_geom_features = function(las,group,feat_list){

  . = .N = N = X = Y = Z = e1 = e2 = e3 = eigen_values = index = NULL

  # c++ function for faster computation
  Rcpp::sourceCpp(code = "
    #include <RcppArmadillo.h>
    // [[Rcpp::depends(RcppArmadillo)]]

    // [[Rcpp::export]]
    SEXP eigen_values(arma::mat A) {
    arma::mat coeff;
    arma::mat score;
    arma::vec latent;
    arma::princomp(coeff, score, latent, A);
    return(Rcpp::wrap(latent));
  }")

  ####### Build three matrix for faster computation
  # sort data by group
  las@data = las@data[order(c(las@data[,..group]))]
  las@data[,N := .N, by = group]

  # store xyz in a matrix
  temp_mat = as.matrix(las@data[,.(X,Y,Z)])

  # for each group, find the minimum and maximum index
  las@data[,index := 1:nrow(las@data)]
  las@data[,':='(min = min(index),max=max(index)),by=group]
  # store in a matrix
  ind_mat = as.matrix(unique(las@data[,.(min,max)]))

  # vector to test the number of points in a group
  # must be >= 2 to allow eigen values computation
  Ntest = c(ind_mat[,2]-ind_mat[,1] >= 2)

  # build an empty matrix to store the outputs
  out_mat = matrix(ncol=3,nrow=nrow(ind_mat))

  # compute eigen values
  for(i in 1:nrow(ind_mat)){
    if(Ntest[i]){
      out_mat[i,1:3] = eigen_values(temp_mat[ind_mat[i,1]:ind_mat[i,2],])
    }
  }

  # build a data.table to cimpute geometric features
  features = data.table::data.table(out_mat[rep(1:nrow(out_mat),times = unique(las@data[,.(get(group),N)])$N),])
  data.table::setnames(features,c("e1","e2","e3"))

  # remove unnecesarry material
  rm(temp_mat,ind_mat,out_mat)
  las@data[,':='(index = NULL, N = NULL, min = NULL, max = NULL)]

  # compute geometric features
  features[,':='(
    Eigenvalue_sum_G = e1+e2+e3,
    Omnivariance_G = (e1*e2*e3)^(1/3),
    Eigenentropy_G = -((e1*log(e1)) + (e2*log(e2)) + (e3*log(e3))),
    Anisotropy_G = (e1-e3)/e1,
    Planarity_G = (e2-e3)/e1,
    Linearity_G = (e1-e2)/e1,
    PCA1_G = e1/(e1+e2+e3),
    PCA2_G = e2/(e1+e2+e3),
    Surface_variation = e3/(e1+e2+e3),
    Sphericity_G = e3/e1,
    Curvature_G = e3/(e1+e2+e3)
  )]

  if(missing(feat_list)){
    # remove eigen values and add geometric features to original point cloud
    features[,':='(e1 = NULL, e2 = NULL, e3 = NULL)]
    las@data = cbind(las@data,features)
  }else{
    # select the required geomterin features and add it to the point cloud
    las@data = cbind(las@data,features[,paste(feat_list,"_G",sep=""), with=FALSE])
  }
  return(las)
}
