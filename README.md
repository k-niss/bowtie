# Bowtie

Construction and visualization of weighted topological overlap matrices based on networks. This, in order to decode the topology of networks, especially very large ones of > 1000 nodes or with a low average shortest path.

![A great test simle](https://i.cbc.ca/1.4986551.1548087152!/fileImage/httpImage/image.jpg_gen/derivatives/16x9_780/smile-emoji.jpg)

## Getting Started

Make sure you have devtools installed and loaded:
```
install("devtools")
library("devtools")
```

Then install the bowtie package directly from github
```
install_github("k-niss/bowtie")
library("bowtie")
```

You may have to install the github repertoire using R in the terminal instead of Rstudio. After that you can just restart Rstudio and load the library in the standard way using library(bowtie).

### Prerequisites

Remember to install and then load the following packages before using the bowtie functions.

```
library(igraph)       # Graph package
library(reshape2)     # Matrix manipulation
library(pbapply)      # Apply functions with multiprocessing and progress bar
library(parallel)     # For multiprocessing
library(RColorBrewer) # Nice colors
```

## Example

How to create a weighted topological matrix of an igraph object

### Tutorial

Create a toy network using the igraph function sample_pa() and visualize it:
```
random_graph           = sample_pa(n=100, power = 1.2, directed=F)
E(random_graph)$weight = runif(n=length(E(random_graph)))
plot(random_graph, 
     vertex.size=2, 
     layout = igraph::layout.gem(random_graph))
```


Calculate pariwise weighted topological overlap (wTO) for all node pairs:
```
wTO_list = wTO.network(node_vector = as.vector(V(random_graph)), 
                       igraph_object = random_graph, 
                       thread_numb = 2)
```


Convert the list format into a symmetric wTO matrix:
```
wTO_matrix = from.list.to.df(wTO_list)
```


Order the columns and rows of the matrix, to make patterns stand out:
```
hclust_object = hclust(as.dist(1-wTO_matrix), method = 'average')
node_order    = hclust_object$labels[hclust_object$order]
```


Visualize the matrix to get an overview of the topology:
```
image(wTO_matrix[node_order,node_order], 
      useRaster = T, 
      col = colorRampPalette(brewer.pal(9,"YlGnBu"))(49))
```

## Authors

* **PhD Student Kristoffer Niss** - *Coding and conceptual work* 
* **Professor Søren Brunak** - *Conceptual work and supervision* 

## Affiliation

Translational Disease Systems Biology Group, Novo Nordisk Foundation center for Protein Research (CPR), University of Copenhagen, 2200 Copenhagen, Denmark. 

## Scientific paper

This R package was published in connection with the article titlted XXX in journal of YYY (2019), doi:.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thank you to Tibor Varga and Grigorii Nos for help with setting up the R package
