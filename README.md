# Bowtie

One Paragraph of project description goes here

## Getting Started

Make sure you have devtools installed and loaded:
```
install("devtools")
library("devtools")
```

Then install the bowtie package directly from github
```
install_github("KrisNiss/bowtie")
library("bowtie")
```

You may have to install the github repertoire using R in the terminal instead of Rstudio. After that you can just restart Rstudio and load the library in the standard way using library(bowtie).

### Prerequisites

Remember to install and then load the following packages before running

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

Explain what these tests test and why

```
## 1) Create random graph with weights
random_graph           = sample_pa(n=100, power = 1.2, directed=F)
E(random_graph)$weight = runif(n=length(E(random_graph)))

plot(random_graph, vertex.size=2, layout = igraph::layout.gem(random_graph))

## 2) Calculate wTO
wTO_list = wTO.network(node_vector = as.vector(V(random_graph)), igraph_object = random_graph, thread_numb = 4)

## 3) Turn list into matrix
wTO_matrix = from.list.to.df(wTO_list)

## 4) Hierarchical ordering of the matrix
hclust_object = hclust(as.dist(1-wTO_matrix), method = 'average')
node_order    = hclust_object$labels[hclust_object$order]

## 5) Visualize matrix
image(wTO_matrix[node_order,node_order], useRaster = T, col = colorRampPalette(brewer.pal(9,"YlGnBu"))(49))
```

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## Scientific paper

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
