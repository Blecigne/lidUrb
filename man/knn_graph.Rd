% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/knn_graph.R
\name{knn_graph}
\alias{knn_graph}
\title{Build a nearest neighbor graph and filter the edges.}
\usage{
knn_graph(
  las,
  k = 10L,
  local_filter = NULL,
  global_filter = NULL,
  max_edge_length = NULL
)
}
\arguments{
\item{las}{a LAS file.}

\item{k}{integer. The number of neighbors to use to construct the graph.}

\item{local_filter}{numeric (optional). If declared, applies a local edge filter that,
allow to remove, for each node, the edges that are longer than
the mean edge \emph{distance+local_filter*sd(mean(distance))}.}

\item{global_filter}{numeric between 0 and 1 (optional). If declared, applies a global
edge filter by removing all edges longer than the \eqn{global_filter^{th}}
quantile of all edges length.}

\item{max_edge_length}{numeric (optional). If declared, removes all edges longer than \emph{max_edge_length}.<}
}
\value{
an data.table containing the nodes ID (their position within the LAS data) and the length of the edge.
}
\description{
Build a nearest neighbor graph and filter the edges.
}
\examples{
\donttest{
# import data
file = system.file("extdata", "tree_no_leaves.las", package="lidUrb")
las = lidR::readLAS(file)

# build the knn graph with no filter
KNN = lidUrb::knn_graph(las, k = 5L)

# plot the graph
lidUrb::plot_graph(las = las, graph = KNN)

# build the knn graph with local filter and all edges longer than 5cm removed
KNN = lidUrb::knn_graph(las, k = 5L,local_filter = 0.5,max_edge_length = 0.05)

# plot the new graph
lidUrb::plot_graph(las = las, graph = KNN)
}
}
