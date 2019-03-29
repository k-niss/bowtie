#' wTO between two nodes in a network
#'
#' Calculate the weighted topological overlap (wTO) between two nodes in an iGraph network
#' @param igraph_obj An igraph network object
#' @param node1_id The name of the first node
#' @param node2_id The name of the second node
#' @param node_vector Vector with the names of all the nodes in the network
#' @param neighbour_list A list containing the neighbor nodes of each node
#' @keywords wTO, weighted topological overlap, pairwise
#' @export
#' @examples
#' 
#' # Create graph
#' random_graph           = sample_pa(n=200, power = 1.2, directed=F)
#' E(random_graph)$weight = runif(n=length(E(random_graph)))
#' 
#' # Calculate wTO
#' wTO_list   = bowtie::wTO.network(node_vector = as.vector(V(random_graph)), igraph_object = random_graph, thread_numb = 4)
#' wTO_matrix = bowtie::from.list.to.df(wTO_list)
#' 
#' # Ordering of nodes
#' hclust_object = hclust(as.dist(1-wTO_matrix), method = 'average')
#' node_order    = hclust_object$labels[hclust_object$order]
#' 
#' # Visualize
#' image(wTO_matrix[node_order,node_order], useRaster = T, col = colorRampPalette(brewer.pal(9,"YlGnBu"))(49))
#' 

wTO.two.nodes = function(igraph_obj, node1_id, node2_id, node_vector, neighbour_list){
  if (node1_id == node2_id){
    w_TO = NA
  }else{
    # calculate denominator
    nbs_of_1 = node_vector[as.vector(neighbour_list[[node1_id]])]
    nbs_of_2 = node_vector[as.vector(neighbour_list[[node2_id]])]
    
    connectivity_of_1 = sum(igraph_obj[node1_id, nbs_of_1]) # get weight
    connectivity_of_2 = sum(igraph_obj[node2_id, nbs_of_2]) # get weight
    
    denominator = sort(c(connectivity_of_1, connectivity_of_2), decreasing=F)
    denominator = denominator[1]
    
    # connection status
    if (are.connected(igraph_obj, node1_id, node2_id) == TRUE){
      direct_connect = igraph_obj[node1_id, node2_id] # get weight
    }else{
      direct_connect = 0
    }
    
    denominator = (denominator+1)-direct_connect
    
    shared_nbs = intersect(nbs_of_1, nbs_of_2)
    numerator  = 0
    
    for (shared_NB in shared_nbs){
      NB_edge_w_of_1 = igraph_obj[node1_id, shared_NB] # get weight
      NB_edge_w_of_2 = igraph_obj[node2_id, shared_NB] # get weight
      product        = NB_edge_w_of_1*NB_edge_w_of_2
      numerator      = numerator + product
    }
    
    numerator = numerator + direct_connect
    w_TO      = numerator/denominator
  }
  return(w_TO)
}

#' wTO between all nodes in a network
#'
#' Calculate the weighted topological overlap for all nodes in a network
#' @param igraph_object An igraph network object
#' @param node_vector Vector with the names of all the nodes in the network
#' @param thread_numb The number of threads to use
#' @keywords wTO, weighted topological overlap, pairwise
#' @export
#' @examples
#' 
#' # Create graph
#' random_graph           = sample_pa(n=200, power = 1.2, directed=F)
#' E(random_graph)$weight = runif(n=length(E(random_graph)))
#' 
#' # Calculate wTO
#' wTO_list   = bowtie::wTO.network(node_vector = as.vector(V(random_graph)), igraph_object = random_graph, thread_numb = 4)
#' wTO_matrix = bowtie::from.list.to.df(wTO_list)
#' 
#' # Ordering of nodes
#' hclust_object = hclust(as.dist(1-wTO_matrix), method = 'average')
#' node_order    = hclust_object$labels[hclust_object$order]
#' 
#' # Visualize
#' image(wTO_matrix[node_order,node_order], useRaster = T, col = colorRampPalette(brewer.pal(9,"YlGnBu"))(49))
#' 

wTO.network = function(igraph_object, node_vector, thread_numb=1){
  
  # multithreading
  clust_obj = makeCluster(thread_numb, type="FORK")
  
  # get all combinations
  cat('# Getting all pairwise combinations\n')
  all_pairwise                = combn(node_vector, 2, simplify = T)
  pairs_length                = dim(all_pairwise)[2]
  all_pairwise_list           = list()
  all_pairwise_list[['node1']] = all_pairwise[1,]
  all_pairwise_list[['node2']] = all_pairwise[2,]
  cat('Found:', pairs_length, '\n')
  
  # all neighbor calculation needed
  cat('# Getting neighbors of each node\n')
  ngb_list          = pblapply(node_vector, function(x) neighbors(igraph_object, v=x), cl=clust_obj)
  names(ngb_list)   = node_vector
  
  ## calculate degree for all
  cat('# Getting the degree per node\n')
  degree_nam_vect        = degree(igraph_object, node_vector)
  names(degree_nam_vect) = node_vector
  
  # calculate topological overlap
  cat('# Calculating pairwise wTO\n')
  all_scores = pbapply(all_pairwise, 2, function(x) wTO.two.nodes(igraph_object, x[1], x[2], node_vector, ngb_list), cl=clust_obj)
  stopCluster(clust_obj)
  
  # append
  all_pairwise_list[['wTO']] = all_scores
  
  # return
  return(all_pairwise_list)
}

#' Tidy-up a list of wTO scores, as given by wTO.network
#'
#' Tidy up a list of wTO scores, returning a complete pairwise matrix
#' @param score_list A list containing three vectors: The names of node1, names of node2 and pairwise wTO scores
#' @keywords wTO, weighted topological overlap, pairwise
#' @export
#' @examples
#' 
#' # Create graph
#' random_graph           = sample_pa(n=200, power = 1.2, directed=F)
#' E(random_graph)$weight = runif(n=length(E(random_graph)))
#' 
#' # Calculate wTO
#' wTO_list   = bowtie::wTO.network(node_vector = as.vector(V(random_graph)), igraph_object = random_graph, thread_numb = 4)
#' wTO_matrix = bowtie::from.list.to.df(wTO_list)
#' 
#' # Ordering of nodes
#' hclust_object = hclust(as.dist(1-wTO_matrix), method = 'average')
#' node_order    = hclust_object$labels[hclust_object$order]
#' 
#' # Visualize
#' image(wTO_matrix[node_order,node_order], useRaster = T, col = colorRampPalette(brewer.pal(9,"YlGnBu"))(49))
#' 

from.list.to.df = function(score_list){
  # convert
  topo_df           = data.frame(matrix(unlist(score_list), ncol = 3, byrow = F), stringsAsFactors = F)
  colnames(topo_df) = c('int1', 'int2', 's')
  topo_df[,3]       = as.numeric(topo_df[,3])
  
  # copy
  copy_df           = topo_df[,c(2,1,3)]
  colnames(copy_df) = c('int1', 'int2', 's')
  
  # add self-combinations
  unq_proteins = unique(c(topo_df[,1], topo_df[,2]))
  self_combos  = data.frame(int1=unq_proteins, int2=unq_proteins, s=rep(NA, length(unq_proteins)))
  
  # assemble
  topo_df = rbind(topo_df, copy_df)
  topo_df = rbind(topo_df, self_combos)
  
  # acast
  full_mat = acast(topo_df, formula = int1~int2, value.var = 's')
  
  return(full_mat)
}
