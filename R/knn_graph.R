#' Build a nearest neighbor graph and filter the edges.
#'
#' @param las a LAS file.
#' @param k integer. The number of neighbors to use to construct the graph.
#' @param local_filter numeric (optional). If declared, applies a local edge filter that,
#'                     allow to remove, for each node, the edges that are longer than
#'                     the mean edge \emph{distance+local_filter*sd(mean(distance))}.
#' @param global_filter numeric between 0 and 1 (optional). If declared, applies a global
#'                      edge filter by removing all edges longer than the \eqn{global_filter^{th}}
#'                      quantile of all edges length.
#' @param max_edge_length numeric (optional). If declared, removes all edges longer than \emph{max_edge_length}.<
#'
#' @return an data.table containing the nodes ID (their position within the LAS data) and the length of the edge.
#' @export
#'
#' @examples
#' \donttest{
#' # import data
#' file = system.file("extdata", "tree_no_leaves.las", package="lidUrb")
#' las = lidR::readLAS(file)
#'
#' # build the knn graph with no filter
#' KNN = lidUrb::knn_graph(las, k = 5L)
#'
#' # plot the graph
#' lidUrb::plot_graph(las = las, graph = KNN)
#'
#' # build the knn graph with local filter and all edges longer than 5cm removed
#' KNN = lidUrb::knn_graph(las, k = 5L,local_filter = 0.5,max_edge_length = 0.05)
#'
#' # plot the new graph
#' lidUrb::plot_graph(las = las, graph = KNN)
#' }
knn_graph = function(las,k = 10L,local_filter = NULL,global_filter = NULL,max_edge_length = NULL){

  . = node_1 = node_2 = quantile = sd = th = NULL

  if(class(las)[1] != "LAS") stop("las must be of type LAS")
  if(!is.integer(k) | k <= 0) stop("k must be an integer >= 1")

  # find the k nearest neibors for each point
  graph = FNN::get.knn(las@data[,1:3],k=k)

  # make a data.table with origin point, the end point and length for each segment of the graph
  graph = data.table::data.table(node_1 = rep(1:nrow(las@data),each=k),node_2 = c(t(graph$nn.index)), length = c(t(graph$nn.dist)))

  # local statistical edge filter
  if(!missing(local_filter)){
    if(!is.numeric(local_filter)) stop("local_filter must be numeric")
    if(local_filter <= 0) stop("local_filter must be > 0")
    graph[,th := mean(length) + local_filter*sd(length),by=node_1]
    graph = graph[length <= th]
    graph[,th := NULL]
  }

  # global statistical edge filter
  if(!missing(global_filter)){
    if(!is.numeric(global_filter)) stop("global_filter must be numeric")
    if(global_filter < 0 | global_filter>1) stop("global_filter must be between 0 and 1")
    graph = graph[length <= quantile(length,global_filter) ]
  }

  # removes edges that are longer than a given length
  if(!missing(max_edge_length)){
    if(!is.numeric(max_edge_length)) stop("max_edge_length must be numeric")
    if(max_edge_length <= 0) stop("max_edge_length must be > 0")
    graph = graph[length <= max_edge_length]
  }

  ########## as global_filter and max_edge_length can remove edges, add it as non conected points
  # look for missing nodes
  to_add = which(! c(1:nrow(las@data)) %in% unique(graph[,c(node_1,node_2)]))
  # if needed add missing nodes connected to themselves
  if(length(to_add) > 0){
    graph = data.table::rbindlist(list(graph,data.table::data.table(node_1 = to_add, node_2 = to_add, length = 0)))
  }

  return(graph[order(node_1,node_2)])
}
