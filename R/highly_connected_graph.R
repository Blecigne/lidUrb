#' Transform a graph into a highly connected graph, i.e. a graph with nearly all
#' points belonging to a cluster that have a minimum Z value lower than a given
#' threshold
#'
#' @param las a LAS file.
#' @param graph a graph as created by the \code{\link{knn_graph}} or \code{\link{delaunay_graph}}
#' @param min_Z_th numeric. The value of Z to decide if a cluster is valid or not.
#' @param max_edge_length numeric. The maximum edge length that can be added to the graph.
#'
#' @return a highly connected grapg
#' @export
#'
#' @examples
#' \donttest{
#' # import data
#' file = system.file("extdata", "four_trees.las", package="lidUrb")
#' las = lidR::readLAS(file,select = "xyz")
#'
#' # reduce point density to 0.1 instead of voxelisation as in Wang et al.
#' las = lidUrb::reduce_point_density(las,0.1)
#'
#' # build a knn  graph
#' KNN = lidUrb::knn_graph(las,local_filter = 0.5)
#'
#' # graph clustering without highly connected graph
#' las = lidUrb::graph_clustering(las,KNN)
#' # random colors to plot not highly connected clusters
#' las@data[,plot_not_high := sample(1:10,1),by = cluster_graph]
#'
#' # transform graph into a highly connected graph
#' KNN = lidUrb::highly_connected_graph(las=las,graph = KNN)
#' # graph clustering based on the highly connected graph
#' las = lidUrb::graph_clustering(las,KNN)
#' las@data[,plot_high := sample(1:10,1),by = cluster_graph]
#'
#' # plot the result
#' lidR::plot(las,color="plot_not_high")
#' lidR::plot(las,color="plot_high")
#' }

highly_connected_graph = function(las,graph,min_Z_th = 1,max_edge_length = 0.5){

  . = D = Z = cluster = index = is.closest = minZ = nearest = node_1 = node_2 = NULL

  if(class(las)[1] != "LAS") stop("las must be of type LAS")
  if(ncol(graph) != 3) stop("graph must have three columns. Not likely a graph.")
  if(!all(names(graph) == c("node_1","node_2","length"))) stop("graph names does not match expected names. Not likely a graph.")
  if(!is.numeric(min_Z_th)) stop("min_Z_th must be numeric")
  if(!is.numeric(max_edge_length)) stop("max_edge_length must be numeric")

  # keep only nodes of the graph
  graph = graph[,.(node_1,node_2)]

  # add index to data
  las@data[,index := 1:nrow(las@data)]

  # first graph clustering
  las@data[1:max(graph[,c(node_1,node_2)]),cluster :=
             igraph::components( # graph clustering
               igraph::make_graph( # build graph
                 c(t(graph[,.(node_1,node_2)]))
               ),mode = "weak"
             )$membership
  ]

  # find the valid vs. invalid clusters and store invalid cluster points in temp
  las@data[,minZ := min(Z),by=cluster]
  las@data[,minZ := minZ-min(Z)]
  temp = las@data[minZ > min_Z_th]

  run = TRUE
  while(run){

    # for points in temp, find the nearest point of the valid clusters
    neib_tab = FNN::get.knnx(data=las@data[minZ <= min_Z_th,1:3],query = temp[,1:3],1)
    temp[,':='(nearest = neib_tab$nn.index, D = neib_tab$nn.dist)]

    # find the the invalid point the closest of a valid point in each invalid cluster
    temp[,is.closest := as.numeric(D == min(D)),by=cluster]

    # add an edge that links the closest point of an invalid cluster to in nearest valid point
    # but only if the valid point is less than max_edge_length
    graph = data.table::rbindlist(list(graph,temp[is.closest == 1 & D <= max_edge_length,.(index,las@data$index[las@data$minZ <= min_Z_th][nearest])]),use.names = FALSE)

    # new graph clustering
    las@data[1:max(graph[,c(node_1,node_2)]),cluster :=
               igraph::components( # graph clustering
                 igraph::make_graph( # build graph
                   c(t(graph[,.(node_1,node_2)]))
                 ),mode = "weak"
               )$membership
    ]

    # stop the loop if there are no invalid clusters close enough of any valid points
    # to link them
    if(min(temp$D) > max_edge_length){
      run = FALSE
    }

    # find the valid vs. invalid clusters and store invalid cluster points in temp
    las@data[,minZ := min(Z),by=cluster]
    las@data[,minZ := minZ-min(Z)]
    temp = las@data[minZ > min_Z_th]
  }

  # add edge length to the graph
  graph[,length := sqrt( (las@data$X[node_1] - las@data$X[node_2])^2 +
                         (las@data$Y[node_1] - las@data$Y[node_2])^2 +
                         (las@data$Z[node_1] - las@data$Z[node_2])^2)]

  # remove non needed columns from the las file
  las@data[,':='(minZ = NULL, cluster = NULL, index = NULL)]

  return(graph)
}
