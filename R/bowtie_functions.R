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

wTO_two_nodes = function(igraph_obj, node1_id, node2_id, node_vector, neighbour_list){
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

wTO_network = function(igraph_object, node_vector, thread_numb=1){

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

from_list_to_df = function(score_list){
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

#' Locate submatrices (protein complexes) on the diagonal of a first-order direct interaction matrix
#'
#' Locate high-scoring square submatrices along the diagonal of a first-order direct interaction matrix.
#' @param full_matrix A first-order direct interaction matrix containing scores between 0-1, where 1 indicates an interaction with highest confidence and 0 indicates no interaction.
#' @param min_conf_score The minimum average score of a submatrix to be accepted as a protein complex (default: 0.75).
#' @param max_complex_size The maximum size of a submatrix (default: 500).
#' @param give_up_expand Give up expanding a submatrix for a given seed when exceeding N proteins in each direction (default: 50). This parameter was implemented to reduce computing time.
#' @keywords submatrices
#' @export
#' @examples
#'

find_complexes = function(full_matrix, min_conf_score=0.75, max_complex_size=500, give_up_expand=50){

  # Collector
  prot_complex_ranges = list()

  # Iterate along diagonal
  for (index in  c(1:dim(full_matrix)[1])){

    biggest_hit         = 0
    hit_count           = 0
    hit_record          = c()

    # We continuesly grow a submatrix around a "protein seed".
    # One hit per seed is reported, i.e. the biggest.
    for (expand in c(1:max_complex_size)){

      low_end  = index-expand
      high_end = index+expand

      # end search when there is not a hit within 5x5 of seed
      if (hit_count == 0 & expand == 5){break}

      # end search if it just doesnt get any better
      if (expand >= give_up_expand & ! 1 %in% tail(hit_record, n=5)){break}

      # if we are not outside matrix, extract submatrix and test
      if (low_end < 1 | high_end >= dim(full_matrix)[1]){
        break
      }else{
        sub_mat      = full_matrix[low_end:high_end, low_end:high_end]
        dim(sub_mat) = NULL
        data_vect    = abs(sub_mat) # the absolute value
        mean_value   = mean(data_vect, na.rm = TRUE)

        # test if interesting
        if (mean_value >= min_conf_score){
          hit_count   = hit_count + 1
          biggest_hit = paste(low_end, high_end, sep=',')
          hit_record  = c(hit_record, 1)
        }else{
          hit_record  = c(hit_record, 0)
        }
      }
    }

    # append biggest hit
    if (!biggest_hit == 0){
      start_index                                = as.numeric(strsplit(biggest_hit, split=',')[[1]][1])
      end_index                                  = as.numeric(strsplit(biggest_hit, split=',')[[1]][2])
      prot_complex_ranges[[as.character(index)]] = c(start_index, end_index)
    }
  }
  return(prot_complex_ranges)
}

#' Merge overlapping submatrices (protein complexes) produced by find.complexes
#'
#' Merge submatrices given by find.complexes if they overlap and apply a minimum size threshold.
#' @param prot_complex_ranges A list of start and end ranges that indicate areas of protein complexes.
#' @param merge_if_overlap Complexes are merged if they have this number of overlapping proteins (default: 1). A value of 1 will produce no overlapping protein complexes.
#' @param min_complex_size The minimum size of a protein complex (default: 8).
#' @return Returns a list of vectors that specifiy the start and end range of each submatrix
#' @keywords submatrices
#' @export
#' @examples
#'

merge_complexes = function(prot_complex_ranges, merge_if_overlap=1, min_complex_size=8){

  prot_complex_ranges_merged = list()
  range_count                = 0

  # iterate over ranges
  for (range_name in names(prot_complex_ranges)){

    complex_range = prot_complex_ranges[[range_name]]

    # if at the first complex range
    if (which(range_name == names(prot_complex_ranges)) == 1){
      start = complex_range[1]
      end   = complex_range[2]
    }

    # if at a complex range in the middle of matrix
    if (which(range_name == names(prot_complex_ranges)) > 1){
      complex_start = complex_range[1]
      complex_end   = complex_range[2]

      # check intersect with current start/end
      overlap = intersect(c(start:end), c(complex_start:complex_end))

      if (length(overlap) >= merge_if_overlap){
        # set new end
        end = complex_end
      }else{
        # this is a new complex, safe old first
        range_count                               = range_count+1
        prot_complex_ranges_merged[[range_count]] = c(start, end)
        start = complex_start
        end   = complex_end
      }
    }

    # if at the final complex range
    if (which(range_name == names(prot_complex_ranges)) == length(prot_complex_ranges)){
      complex_start = complex_range[1]
      complex_end   = complex_range[2]

      # check intersect with current start/end
      overlap = intersect(c(start:end), c(complex_start:complex_end))

      if (length(overlap) > merge_if_overlap){
        # set new end
        end = complex_end
      }else{
        # this is a new complex, save old first
        range_count = range_count+1
        prot_complex_ranges_merged[[range_count]] = c(start, end)
        start = complex_start
        end   = complex_end
      }
      # final save
      range_count = range_count+1
      prot_complex_ranges_merged[[range_count]] = c(start, end)
    }
  }

  # apply min size filter
  count    = 0
  prot_complex_ranges_merged_filt = list()
  for (range in prot_complex_ranges_merged){
    if (range[2]-range[1] >= min_complex_size){
      count = count+1
      prot_complex_ranges_merged_filt[[count]] = range
    }
  }

  return(prot_complex_ranges_merged_filt)
}

#' Cuts a vector of indicies into sub-vectors that are consecutive (increase by +1)
#'
#' Cuts a vector of indicies into sub-vectors that are consecutive (increase by +1)
#' @param string_vect A vector of indicies
#' @return Returns a list of vectors containing indicies that are consecutive
#' @keywords submatrices
#' @export
#' @examples
#'

cut_into_chunks = function(string_vect){
  diff_vector  = c(0, string_vect)-c(string_vect,0)
  index_list   = list()

  # merge
  ct = 0
  subset_indicies = c()
  for (i in diff_vector){
    ct = ct+1
    if (ct > 1 & i < -1){
      subset_indicies = c(subset_indicies, c(ct-1, ct))
    }
  }

  subset_indicies = c(1, subset_indicies, length(string_vect))
  subset_pairs    = split(subset_indicies, ceiling(seq_along(subset_indicies)/2))
  index_list      = list()

  for (pair in subset_pairs){index_list[[paste(pair, collapse='_')]] = string_vect[pair[1]:pair[2]]}
  return(index_list)
}

#' Find sub-vectors of a certain minimum size and average value in a vector
#'
#' Using a sliding window approach, this function returns  sub-vectors of a minimum size and average value found in a numeric vector containing values between 0-1
#' @param input_vector A numeric vector containing values between 0 and 1.
#' @param window_length The minimum size of the desired sub-vectors.
#' @param min_value The minimum average value of the sub-vectors.
#' @return Returns a list of indicies that specify the location of found sub-vectors in the input vector
#' @keywords sub-vectors
#' @export
#' @examples
#'

consecutive_str = function(input_vector, window_length, min_value){

  # problem
  if (length(input_vector) < window_length){
    cat('ERROR: The window length exceeds the size of the input vector\n')
    index_list = NULL
  }else{
    # sliding window
    index_range        = 1:(length(input_vector)-10)
    vector_of_indicies = c()

    for (i_start in index_range){
      i_end        = i_start+(window_length-1)
      sub_string   = input_vector[i_start:i_end]
      string_score = sum(sub_string, na.rm = T)/window_length

      if (string_score >= min_value){vector_of_indicies = c(vector_of_indicies, c(i_start:i_end))}
    }

    unq_indicies = sort(unique(vector_of_indicies))
    index_list   = list()

    if (length(unq_indicies) > 5){index_list = cut_into_chunks(unq_indicies)}
  }
  return(index_list)
}

#' Find strings of high-confidence interactions between protein complexes in a first-order direct interaction matrix
#'
#' Find strings of high-confidence interactions between protein complexes in a first-order direct interaction matrix.
#' @param full_matrix A first-order direct interaction matrix containing scores between 0-1, where 1 indicates an interaction with highest confidence and 0 indicates no interaction.
#' @param protein_complex_areas Start and end ranges along the matrix diagonal annotated to be protein complex submatrices.
#' @param int_fan_min_size Minimum number of interactions in a bow-tie interaction fan (default: 10).
#' @param int_fan_avg_conf_score Minimum average confidence score across the interactions in a bow-tie interaction fan (default: 0.9).
#' @return Returns a list of "knot" proteins and their interactions fans indicated by start and end range indicies.
#' @keywords submatrices
#' @export
#' @examples
#'

find_bowties = function(full_matrix, protein_complex_areas, int_fan_min_size=10, int_fan_avg_conf_score=0.9){

  # 1A) Find lines of correct size and score
  cat('# Searching for interaction fans of bow-ties...\n')
  collected_bowties = list()
  for (row_count in 1:dim(full_matrix)[1]){
    row_vect   = full_matrix[row_count,]

    if (sum(row_vect) >= (int_fan_min_size*int_fan_avg_conf_score)){

      # find potential interaction fans
      index_list = consecutive_str(row_vect, int_fan_min_size, int_fan_avg_conf_score)

      # save interaction fans
      ct = 1
      for (index_set in index_list){
        knot = all_proteins[row_count]
        collected_bowties[[knot]][[ct]] = c(min(index_set), max(index_set))
        ct = ct+1
      }
    }
  }

  # 1B) Get protein complex ranges
  complex_areas         = list()
  complex_ranges_vector = unlist(lapply(prot_complex_ranges_merged_filt, function(x) paste(x[1], x[2], sep='_')))
  for (i in prot_complex_ranges_merged_filt){
    start    = i[1]
    end      = i[2]
    proteins = colnames(int_matrix_all)[start:end]

    for (prot in proteins){
      complex_areas[[prot]] = c(start, end)
    }
  }

  # 2) Interaction fans must not be within protein complexes
  cat('# Filter 1: Keep interaction fans that are outside protein complexes...\n')
  collected_bowties_filt1 = list()
  for (knot_protein in names(collected_bowties)){

    if (knot_protein %in% names(complex_areas)){
      own_complex       = complex_areas[[knot_protein]]
      own_complex_range = c(own_complex[1]:own_complex[2])
      count             = 1

      for (index_set in collected_bowties[[knot_protein]]){
        bowtie_fan_range_old = index_set[1]:index_set[2]
        bowtie_fan_range_new = bowtie_fan_range_old[!bowtie_fan_range_old %in% own_complex_range]
        range_chunks         = cut_into_chunks(bowtie_fan_range_new)

        for (bow_tie_range in range_chunks){
          if (length(bow_tie_range) >= int_fan_min_size){
            collected_bowties_filt1[[knot_protein]][[count]] = c(min(bow_tie_range), max(bow_tie_range))
            count = count+1
          }
        }
      }
    }
  }

  # 3) Interaction fans should be located between two protein complexes
  cat('# Filter 2: Keep interaction fans that are located between two protein complexes...\n')
  collected_bowties_filt2 = list()
  for (knot_protein in names(collected_bowties_filt1)){
    ct = 1
    for (index_set in collected_bowties_filt1[[knot_protein]]){
      bowtie_fan_range = index_set[1]:index_set[2]

      # is the knot-protein in a protein complex?
      knot_in_complex = 'yes' # already checked above
      complex_numb1   = which(complex_ranges_vector %in% paste(complex_areas[[knot_protein]], collapse='_'))

      # is the terminal-fan proteins in a protein complex?
      fans_in_complex = 'no'
      complex_numb2   = 0
      c_count         = 0
      for (prot_comp_range in prot_complex_ranges_merged_filt){
        c_count      = c_count+1
        all_indicies = c(prot_comp_range[1]:prot_comp_range[2])
        status       = bowtie_fan_range %in% all_indicies
        no_of_true   = unname(table(status)['TRUE'])

        if(is.na(no_of_true) == FALSE){
          if (no_of_true >= int_fan_min_size){
            fans_in_complex = 'yes'
            complex_numb2   = c_count
          }
        }
      }

      # check difference
      complexes_are_different = 'no'
      if (complex_numb1 != complex_numb2){
        complexes_are_different = 'yes'
      }

      if (knot_in_complex == 'yes' & fans_in_complex == 'yes' & complexes_are_different == 'yes'){
        collected_bowties_filt2[[knot_protein]][[ct]] = c(index_set[1], index_set[2])
        ct = ct+1
      }
    }
  }

  return(collected_bowties_filt2)
}
