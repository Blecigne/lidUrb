#' Build and filter a delaunay graph
#'
#' @param las aLAS file.
#' @param global_filter numeric between 0 and 1 (optional). If declared, applies a global
#'                      edge filter by removing all edges longer than the \eqn{global_filter^{th}}
#'                      quantile of all edges length.
#' @param max_edge_length numeric (optional). If declared, removes all edges longer than \emph{max_edge_length}.
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
#' delaunay = lidUrb::delaunay_graph(las)
#'
#' # plot the graph
#' lidUrb::plot_graph(las = las, graph = delaunay)
#'}

delaunay_graph = function(las,global_filter = 0.8,max_edge_length = NULL){

  . = node_1 = node_2 = quantile = NULL

  if(class(las)[1] != "LAS") stop("las must be of type LAS")
  if(!is.numeric(global_filter)) stop("global_filter must be numeric")

  # delaunay triangulation using the geometry package
  graph = geometry::delaunayn(las@data[,1:3],output.options=TRUE)$tri

  # transform in the right form and add edge length
  graph = unique(data.table::data.table(node_1 = rep(graph[,1],each = 3), node_2 = c(t(graph[,2:4]))))
  graph[,length := sqrt( (las@data$X[node_1] - las@data$X[node_2])^2 +
                         (las@data$Y[node_1] - las@data$Y[node_2])^2 +
                         (las@data$Z[node_1] - las@data$Z[node_2])^2) ]

  # global statistical edge filter
  graph = graph[length <= quantile(length,global_filter)]

  # removes edges that are longer than a given length
  if(!missing(max_edge_length)){
    if(!is.numeric(max_edge_length)) stop("max_edge_length must be numeric")
    if(max_edge_length <= 0) stop("max_edge_length must be > 0")
    graph = graph[length <= max_edge_length]
  }

  # as global_filter and max_edge_length can remove edges, add it as non conected points
  # look for missing nodes
  to_add = which(! c(1:nrow(las@data)) %in% unique(graph[,c(node_1,node_2)]))
  # if needed add missing nodes connected to themselves
  if(length(to_add) > 0){
    graph = data.table::rbindlist(list(graph,data.table::data.table(node_1 = to_add, node_2 = to_add, length = 0)))
  }

  return(graph)
}
