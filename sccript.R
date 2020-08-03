load('dataset.RData')


library(igraph)

x <- graph_from_adjacency_matrix(X, mode = 'undirected')
y <- graph_from_adjacency_matrix(Y, weighted = TRUE, mode = 'undirected')


list.edge.attributes(y)


# visualise the static and weigthed network
igraph.options(vertex.label = NA,
               vertex.size = 2, 
               edge.width = 0.5, 
               edge.color = "grey50", 
               edge.curved = 0.1)

plot(x, layout = layout)
plot(y, layout = layout, edge.label = E(y)$weight, edge.label.cex = 0.5)

# Binary degree of X and the weighted degrees of Y. 
d_bin <- degree(x) 
d_weight <- strength(y)

# 5 nodes with the highest binary degree and their weighted degrees.
best_bin <- order(d_bin, decreasing = TRUE)[1:5]
d_weight[best_bin]




