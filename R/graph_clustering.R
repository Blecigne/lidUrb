#' Cluster data base on k nearest neighbor graph. Non connected components are
#' parts of different objects.
#'
#' @param las a LAS file. NOTE: the data must be ordered similarly than when the graph was created.
#' @param graph a graph produced with the \code{\link{knn_graph}} or \code{\link{delaunay_graph}} functions.
#'
#' @return a vector with the points cluster
#'
#' @importFrom data.table :=
#'
#' @export
#'
#' @examples
#' \donttest{
#' # import data
#' file = system.file("extdata", "urban.las", package="lidUrb")
#' las = lidR::readLAS(file)
#'
#' # remove the ground
#' las@data = las@data[Classification == 0]
#'
#' # build the knn graph with no filter
#' KNN = lidUrb::knn_graph(las, k = 5L)
#'
#' # clustering
#' las = lidUrb::graph_clustering(las = las, graph = KNN)
#'
#' # add random color to clusters
#' las@data[,cluster_plot := sample(1:20,1),by=cluster_graph]
#'
#' # plot
#' lidR::plot(las,color="cluster_plot")
#' }

graph_clustering = function(las,graph){

  . = node_1 = node_2 = cluster_graph = NULL

  if(class(las)[1] != "LAS") stop("las must be of type LAS")
  if(ncol(graph) != 3) stop("graph must have three columns. Not likely a graph.")
  if(!all(names(graph) == c("node_1","node_2","length"))) stop("graph names does not match expected names. Not likely a graph.")

  las@data[1:max(graph[,c(node_1,node_2)]),cluster_graph :=
             igraph::components( # graph clustering
               igraph::make_graph( # build graph
                 c(t(graph[,.(node_1,node_2)]))
               ),mode = "weak"
             )$membership
  ]

  return(las)
}
