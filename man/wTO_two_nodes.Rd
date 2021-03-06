% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bowtie_functions.R
\name{wTO_two_nodes}
\alias{wTO_two_nodes}
\title{wTO between two nodes in a network}
\usage{
wTO_two_nodes(igraph_obj, node1_id, node2_id, node_vector, neighbour_list)
}
\arguments{
\item{igraph_obj}{An igraph network object}

\item{node1_id}{The name of the first node}

\item{node2_id}{The name of the second node}

\item{node_vector}{Vector with the names of all the nodes in the network}

\item{neighbour_list}{A list containing the neighbor nodes of each node}
}
\description{
Calculate the weighted topological overlap (wTO) between two nodes in an iGraph network
}
\examples{

# Create graph
random_graph           = sample_pa(n=200, power = 1.2, directed=F)
E(random_graph)$weight = runif(n=length(E(random_graph)))

# Calculate wTO
wTO_list   = bowtie::wTO.network(node_vector = as.vector(V(random_graph)), igraph_object = random_graph, thread_numb = 4)
wTO_matrix = bowtie::from.list.to.df(wTO_list)

# Ordering of nodes
hclust_object = hclust(as.dist(1-wTO_matrix), method = 'average')
node_order    = hclust_object$labels[hclust_object$order]

# Visualize
image(wTO_matrix[node_order,node_order], useRaster = T, col = colorRampPalette(brewer.pal(9,"YlGnBu"))(49))

}
\keyword{overlap,}
\keyword{pairwise}
\keyword{topological}
\keyword{wTO,}
\keyword{weighted}
