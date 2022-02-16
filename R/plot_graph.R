
#' Plot a graph.
#'
#' @param las a LAS file.
#' @param graph a graph produced with the \code{\link{knn_graph}} or \code{\link{delaunay_graph}} functions.
#' @param graph_color the color of the graph edges to plot.
#' @param ... additionnal parameters to pass to the lidR \code{\link[lidR]{plot}} function.
#'
#' @return plot the graph.
#' @export
#'
#' @examples
#' \donttest{
#' # import data
#' file = system.file("extdata", "tree_no_leaves.las", package="lidUrb")
#' las = lidR::readLAS(file)
#'
#' # build a knn graph
#' KNN = lidUrb::knn_graph(las, k = 5L)
#'
#' # plot the graph
#' lidUrb::plot_graph(las = las, graph = KNN)
#'
#' # change the size of the nodes
#' lidUrb::plot_graph(las = las, graph = KNN, size = 5)
#' }


plot_graph = function(las,graph,graph_color = "white",...){

  . = node_1 = node_2 = NULL

  if(class(las)[1] != "LAS") stop("las must be of type LAS")
  if(ncol(graph) != 3) stop("graph must have three columns. Not likely a graph.")
  if(!all(names(graph) == c("node_1","node_2","length"))) stop("graph names does not match expected names. Not likely a graph.")

  # plot the nodes
  lidR::plot(las,clear_artifacts = FALSE,...)
  # plot the edges
  rgl::segments3d(las@data[c(t(graph[,.(node_1,node_2)]))],col=graph_color)
}



