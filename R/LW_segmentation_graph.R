#' Computes indices for leaf / wood segmentation
#'
#' @param las a LAS file.
#' @param k integer. The number of nearest neighbors to build the knn graph.
#' @param search_radius numeric. The searching distance to compute geometric features.
#'                      Note that \code{2*search_radius} is used to compute linearity
#'                      that which provides better estimates.
#' @param features_similarity numeric. Sets the thresholds in geometric features used to trim the knn graph.
#'                            It sets how similar the points in a cluster must be. Typically between 0.1 and 0.2.
#' @param min_cluster_size integer. The minimal size of a cluster required to compute the indices. Points located
#'                         in small clusters receive a p_wood and SoD of their nearest well classified point.
#' @param reclass_trunk_th numeric. The minimal Threshold to reclass the lower clusters (i.e. the clusters with at least
#'                         one point with \code{Z <= reclass_trunk_th}) as wood, i.e. a \code{p_wood = SoD = 1}. A negative
#'                         value desable the reclassification.
#'
#' @return a LAS with two additionnal fields, \code{p_wood} and \code{SoD}, corresponding to a wood classification index. \code{p_wood}
#'         was introduced in Wang et al. 2019 and range from 0 to 1 with 1 being very high probability of wood. \code{SoD} was introduced
#'         in Wan et al. 2020 and range between -1 and 1, higher probability of wood class close to 1.
#'
#' @description This algorithm uses the graph clustering framework presented in Wang et al. 2019 to compute the geometric features of clusters from
#'              which the p_wood and SoD indices are computed. However, unlike the Wang et al. 2019 method, this algorithm does not use
#'              post-clustering partitioning, which is a long process, and replace it by the inclusion of two additional geometric
#'              features to remove inconstant edges from the graph. However, the three geometric features are all controlled by one parameter
#'              so it does not add any additional parameter to the original method. The class regularization step was also removed and replaced
#'              by the reclassification of points in small clusters.
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
#' las = LW_segmentation_graph(las)
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

LW_segmentation_graph = function(las, k = 10L,search_radius = 0.05, features_similarity = 0.15, min_cluster_size = 10L, reclass_trunk_th = 0.5){

    . = .N = Linearity_G = Linearity_n1 = Linearity_n2 = N = Planarity_G = SoD =
    Sphericity_G = Surf_var_n1 = Surf_var_n2 = Verticality_n1 = Verticality_n2 = X =
    Y = Z = cluster_graph = max_length = node_1 = node_2 = p_wood = sd = NULL

    if(class(las)[1] != "LAS") stop("las must be of type LAS")
    if(!is.integer(k) | k <= 0) stop("k must be an integer >= 1")
    if(!is.integer(min_cluster_size) | k <= 0) stop("min_cluster_size must be an integer >= 1")
    if(!is.numeric(search_radius) | search_radius <= 0) stop("search_radius must be numeric >= 0")
    if(!is.numeric(features_similarity) | features_similarity <= 0) stop("features_similarity must be numeric >= 0")
    if(!is.numeric(reclass_trunk_th)) stop("reclass_trunk_th must be numeric")

    # compute geometric features
    las = lidUrb::geom_features(las,search_radius,features_list = c("Verticality","Surface_variation"))
    las = lidUrb::geom_features(las,search_radius*2,features_list = "Linearity")

    ################# KNN graph
    ##### buid KNN graph with distance
    graph = lidUrb::knn_graph(las,k,local_filter = 1)

    # add geometric features to each node of the graph for graph trimming
    graph[,':='(Verticality_n1 = las@data$Verticality[node_1],
                Verticality_n2 = las@data$Verticality[node_2],
                Surf_var_n1 = las@data$Surface_variation[node_1],
                Surf_var_n2 = las@data$Surface_variation[node_2],
                Linearity_n1 = las@data$Linearity[node_1],
                Linearity_n2 = las@data$Linearity[node_2]
                )]

    # add the maximum distance for each point (needed for graph trimming)
    graph[,max_length := max(length),by=node_1]

    ##### trim the graph according
    graph = graph[
      abs(Verticality_n1-Verticality_n2) < features_similarity &
      abs(Linearity_n1-Linearity_n2) < features_similarity/2 &
      abs(Surf_var_n1-Surf_var_n2) < features_similarity/2 &
      max_length < mean(max_length) + sd(max_length)
    ]

    ##### cluster LAS based of the trimmed graph
    las = lidUrb::graph_clustering(las,graph[,.(node_1,node_2,length)])

    ##### compute cluster geometric features
    # keep clusters large enough for confidence
    las@data[,N := .N, by=cluster_graph]
    las@data[N <= min_cluster_size | is.na(cluster_graph), cluster_graph := 0]
    # compute features
    las = lidUrb::group_geom_features(las,group = "cluster_graph", feat_list = c("Linearity","Planarity","Sphericity"))

    ##### compute the probability of wood classification as in Wang et al. 2019
    las@data[,p_wood := 0]
    for(i in seq(0.01,1,0.01)){
      las@data[Linearity_G >= i & cluster_graph > 0,p_wood := p_wood+0.01]
    }

    ##### compute the Significance Of Difference index as in Wan et al. 2020
    las@data[,SoD := Linearity_G+(1-Linearity_G)*(Linearity_G-max(c(Planarity_G,Sphericity_G))),by = cluster_graph]

    # reclass if reclass_tree_base is TRUE the tree base recieve a p_wood of 1 ans a SoD of 1
    if(reclass_trunk_th > 0){
      las@data[cluster_graph %in% unique(las@data[Z<=min(Z)+reclass_trunk_th,cluster_graph]),':='(p_wood = 1,SoD = 1)]
    }

    # class points in small clusters
    if(min(las@data$cluster_graph) == 0){
      neib = FNN::knnx.index(data = las@data[cluster_graph > 0,.(X,Y,Z)],query = las@data[cluster_graph == 0,.(X,Y,Z)],1)
      las@data[cluster_graph == 0,':='(
        p_wood = las@data$p_wood[las@data$cluster_graph > 0][neib],
        SoD = las@data$SoD[las@data$cluster_graph > 0][neib])]
    }

    las@data[,':='(Verticality = NULL, Surface_variation = NULL, Linearity = NULL, cluster_graph = NULL, Linearity_G = NULL,
              Planarity_G = NULL, Sphericity_G = NULL)]

    return(las)
}


