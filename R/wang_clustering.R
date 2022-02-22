#' Implement the graph based segmentation presented in Wang et al.
#'
#' @param las a LAS file. NOTE: the data must be ordered similarly than when the graph was created.
#' @param graph a graph produced with the \code{\link{knn_graph}} or \code{\link{delaunay_graph}} functions.
#' @param heigth_th_merge_roots numeric. Defines the height of a root relative to the ground to be considered as valid.
#' @param distance_th_merge_close_roots numeric. Defines the distance between two roots to be merged. It defines the
#'                                      euclidean distance threshold, the geodesic distance threshold = 3*distance_th_merge_close_roots
#'                                      as suggested by Wang et al.
#' @param correct_elevation character. If \code{correct_elevation = "local"} a local correction of point elevation is performed by locating
#'                          the tree bases and using it as reference to correct the points elevation. If \code{correct_elevation = "global"} the
#'                          elevation is corrected relative to the minmum Z value of the point cloud that is set to 0.
#'                          If \code{correct_elevation = "none"} the Z value is used without transformation (recommended if the point cloud
#'                          elevation was previously normalized).
#' @param merge_non_connected logical. If TRUE objects that are disconnected in the graph are merged to the nearest cluster.
#'
#' @references Wang, D., Liang, X., Mofack, G. I., & Martin-Ducup, O. (2021). Individual tree extraction from terrestrial laser
#'             scanning data via graph pathing. Forest Ecosystems, 8(1), 1-11.
#'
#' @return the las file
#' @export
#'
#' @examples
#' \donttest{
#' ######################################
#' # Wang et al. tree segmentation method
#' ######################################
#' # import data
#' file = system.file("extdata", "four_trees.las", package="lidUrb")
#' las = lidR::readLAS(file)
#'
#' # reduce point density to 0.1 instead of voxelisation as in Wang et al.
#' las = lidUrb::reduce_point_density(las,0.1)
#'
#' # build an hybrid graph mixing a knn graph and a dlaunay graph
#' # with Wang et al parameters
#' hybrid_graph = data.table::rbindlist(list(
#'   lidUrb::knn_graph(las,k=10L,local_filter = 1),
#'   lidUrb::delaunay_graph(las,global_filter = 0.8)
#' ))
#'
#' # run the downward clustering following the hybrid graph, parameters are
#' # guessed after Wang et al. paper
#' las_sub=lidUrb::wang_clustering(las=las,
#'                                 graph=hybrid_graph,
#'                                 heigth_th_merge_roots = 1,
#'                                 distance_th_merge_close_roots = 1)
#'
#' # plot the result
#' lidR::plot(las_sub,color="cluster_wang")
#'
#' ########################
#' # An alternative example
#' ########################
#' # import data
#' file = system.file("extdata", "four_trees.las", package="lidUrb")
#' las = lidR::readLAS(file)
#'
#' # reduce point density to 0.1 instead of voxelisation as in Wang et al.
#' las = lidUrb::reduce_point_density(las,0.1)
#'
#' # build a highly connected knn graph
#' KNN = lidUrb::knn_graph(las,local_filter = 0.5)
#' KNN = lidUrb::highly_connected_graph(las,KNN)
#'
#' # run the downward clustering following the hybrid graph, parameters are
#' # guessed after Wang et al. paper
#' las_sub=lidUrb::wang_clustering(las=las,
#'                                 graph=KNN,
#'                                 heigth_th_merge_roots = 1,
#'                                 distance_th_merge_close_roots = 1)
#'
#' # plot the result
#' lidR::plot(las_sub,color="cluster_wang")
#' }
wang_clustering = function(las,graph,heigth_th_merge_roots = 0.5,distance_th_merge_close_roots = 1,correct_elevation = "global",merge_non_connected = TRUE){

  . = .GRP = Dij_dist = Euc_dist = X = Y = Z = Z_node1 = Z_node2 = as.dist = cutree = dist =
  index = is_min = na.omit = nearest_cl = node = node_1 = node_2 = potential = root = root_orig =
    ..group = cl = cluster = is.min = norm_Z = NULL

  if(class(las)[1] != "LAS") stop("las must be of type LAS")
  if(ncol(graph) != 3) stop("graph must have three columns. Not likely a graph.")
  if(!all(names(graph) == c("node_1","node_2","length"))) stop("graph names does not match expected names. Not likely a graph.")
  if(!is.numeric(heigth_th_merge_roots)) stop("heigth_th_merge_roots must be numeric")
  if(!is.numeric(distance_th_merge_close_roots)) stop("distance_th_merge_close_roots must be numeric")
  if(!is.logical(merge_non_connected)) stop("merge_non_connected must be logical")
  if(!is.character(correct_elevation)) stop("correct_elevation must be a character string with 'local','global' or 'none'.")
  if(!correct_elevation %in% c("local","global", "none")) stop("correct_elevation must be either 'local','global' or 'none'.")

  # optionally perform an elevation correction to maximize the chance that roots are low
  if(correct_elevation == "local"){
    # subset the lower 2 meters of the point cloud
    sub = las@data[Z <= min(Z)+2]
    # cluster the subset point cloud
    sub[,cl := dbscan::dbscan(sub[,1:2],eps=5,1)$cluster]
    # find the lower point for each cluster
    sub[,is.min := Z == min(Z),by = cl]
    sub = sub[is.min == TRUE]
    # attach each point of the original data to the nearest lower point
    las@data[,cluster := sub$cl[FNN::knnx.index(data = sub[,1:2],query = las@data[,1:2],1)]]
    rm(sub)
    # define a new field with corrected elevation
    las@data[,norm_Z := Z-min(Z), by = cluster]
  }
  if(correct_elevation == "global"){
    las@data[,norm_Z := Z-min(Z)]
  }
  if(correct_elevation == "none"){
    las@data[,norm_Z := Z]
  }

  # keep index to sort after using key matching
  las@data[,index := 1:nrow(las@data)]

  ##### trim the graph to keep downward relations only
  # add Z values for each node
  graph[,':='(Z_node1 = las@data$Z[graph$node_1] , Z_node2 = las@data$Z[graph$node_2])]
  # node 1 is the one with higher Z
  graph = data.table::rbindlist(list(graph[Z_node1 > Z_node2,],
                                     graph[Z_node1 < Z_node2,.(node_2,node_1,length,Z_node2,Z_node1)]),
                                use.names = FALSE)

  # find and keep the edge that link the lower neighbor
  graph[,is_min := as.numeric(Z_node2 == min(Z_node2)),by = node_1]
  graph_trim = graph[is_min == 1,.(node_1,node_2)]

  ##### add non connected points, i.e. the points that does not connect to neighbors with smaller Z
  graph_trim = data.table::rbindlist(list(graph_trim,data.table::data.table(node_1 = las@data$index,node_2=0)))
  graph_trim = graph_trim[,max(node_2),by=node_1]
  data.table::setnames(graph_trim,c("node_1","node_2"))

  ##### remove loops
  graph_trim = graph_trim[node_1 != node_2]

  ##### sort the graph so node_1 correspond to the points index in original data
  graph_trim = graph_trim[order(graph_trim$node_1)]

  ##### find the root for each point by following the paths downward until the path ends
  las@data[,':='(root = 1:nrow(las@data),potential = 1)]
  while(max(las@data$potential)>0){
    las@data[,potential := graph_trim$node_2[root]]
    las@data[potential > 0, root := potential]
  }

  las@data[,root_orig := root] # store root for next step

  rm(graph_trim)
  graph[,':='(Z_node1 = NULL, Z_node2 = NULL, is_min = NULL)]

  ############################################
  # reconnect non connected roots to the graph
  ############################################

  disconected = unique(las@data$root)[which(!unique(las@data$root) %in% c(graph$node_1,graph$node_2))]

  # if there are disconectes root -> attach it to the nearest node
  if(length(disconected) > 0){
    temp_graph = data.table::data.table(
      node_1 = disconected,
      node_2 = las@data$index[!las@data$root %in% disconected][FNN::knnx.index(data = las@data[!root %in% disconected, 1:3], query = las@data[disconected,1:3],1)]
    )
    temp_graph[,length := sqrt( (las@data$X[node_1] - las@data$X[node_2])^2 +
                                  (las@data$Y[node_1] - las@data$Y[node_2])^2 +
                                  (las@data$Z[node_1] - las@data$Z[node_2])^2) ]

    graph = data.table::rbindlist(list(graph,temp_graph))

    rm(temp_graph)
  }

  ##########################
  ##### trim and merge roots
  ##########################

  ################ merge each root to its nearest valid root

  ##### identify correct and not correct root nodes based on their Z value
  valid = data.table::data.table(node = unique(las@data[norm_Z <= heigth_th_merge_roots,root]))
  not_valid = unique(las@data[,root])
  not_valid = data.table::data.table(node = not_valid[which(!not_valid %in% valid$node)])

  ##### merge invalid roots to valid roots by shortest pathing
  # distance matrix based on Dijkstra
  dist_mat = cppRouting::get_distance_matrix(cppRouting::makegraph(graph[,.(node_1,node_2,length)],directed = F),from = not_valid$node, to = valid$node)
  dist_mat[is.na(dist_mat)] <- 9999 # replace NAs by a very large value
  # find the valid root with shortest path but only for thoose that are connected to the graph
  not_valid[Rfast::rowMins(dist_mat,value = TRUE)<9999,valid := valid$node[Rfast::rowMins(dist_mat)[Rfast::rowMins(dist_mat,value = TRUE)<9999]]]

  # modify invalid roots in the original data
  data.table::setkey(las@data,root)
  data.table::setkey(not_valid,node)
  las@data[not_valid, root := valid]

  # reorder the data to keep original order
  las@data = las@data[order(las@data$index)]

  rm(valid,not_valid)

  ################ merge roots that are close to each other in terms of both euclidean distance and geodesic distance
  roots = data.table::data.table(node = na.omit(unique(las@data$root)))
  # merge roots that are close (<= merge_d) in terms of euclidean distance
  roots[,Euc_dist := cutree(fastcluster::hclust(dist(las@data[roots$node,1:3]),method = "single"), h = distance_th_merge_close_roots)]

  # clustering using dijkstra
  dist_mat = cppRouting::get_distance_matrix(cppRouting::makegraph(graph[,.(node_1,node_2,length)],directed = F),from = roots$node, to = roots$node)
  dist_mat[is.na(dist_mat)] <- 9999 # replace NAs by a very large value
  # merge roots that are close (<=3*merge_d) in terms of geodesic distance
  roots[,Dij_dist := cutree(fastcluster::hclust(as.dist(dist_mat),method = "single"), h = 3*distance_th_merge_close_roots)]

  # combine the two merging type to produce the final clustering
  roots[,merge := .GRP, by=.(Euc_dist,Dij_dist)]

  # add the final clusters to the original data
  data.table::setkey(las@data,root)
  data.table::setkey(roots,node)
  las@data[roots, root := merge]

  # reorder the data to keep original order
  las@data = las@data[order(las@data$index)]

  rm(roots,dist_mat)

  ################### merge parts of the trimmed graph that were not connected to the global graph
  if(merge_non_connected){
    to_connect = las@data[unique(las@data[is.na(root),root_orig])] # non connected roots

    if(nrow(to_connect) > 0){
      # cluster of the nearest point
      to_connect[,nearest_cl := las@data$root[!is.na(las@data$root)][FNN::knnx.index(data = las@data[!is.na(root),.(X,Y,Z)], query = to_connect[,.(X,Y,Z)], 1)]]

      # add the cluster in the data
      data.table::setkey(las@data,root_orig)
      data.table::setkey(to_connect,root_orig)
      las@data[to_connect, root := nearest_cl]

      # reorder the data to keep original order
      las@data = las@data[order(las@data$index)]

      rm(to_connect)
    }
  }

  las@data[,':='(index = NULL,potential = NULL,root_orig = NULL)]
  data.table::setnames(las@data,old="root",new="cluster_wang")

  return(las)
}
