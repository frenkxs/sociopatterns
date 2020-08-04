load('dataset.RData')


library(igraph)

x <- graph_from_adjacency_matrix(X, mode = 'undirected')
y <- graph_from_adjacency_matrix(Y, weighted = TRUE, mode = 'undirected')


list.edge.attributes(y)
E(y)$weight

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


# Create a visual representation of the network X, where the size of the nodes is proportional to the weighted degree, and the nodes are coloured according to time_frame. NOTE: all of the network plots should use layout for the coordinates of the nodes.

igraph.options(vertex.label = NA,
               vertex.size = d_weight * 0.05, 
               vertex.color = time_frame,
               edge.width = 0.5, 
               edge.color = "grey50", 
               edge.curved = 0.1)

plot(x, layout = layout)


# Consider the distribution of the binary degrees. Create a graphical representation of the degree distribution both in natural scale and in log-log scale. Is there any evidence of a powerlaw behaviour in the tail of the distribution?

hist(d_bin, col="blue", 
     xlab="Degree", ylab="Frequency", 
     main="Degree Distribution")

d_dist <- degree.distribution(x)
d <- 0:max(d_bin) 
ind <- (d_dist != 0)
t <- d[ind]
u <- d_dist[ind]
barplot(log(d[ind]), log(d_dist[ind]), col = "grey50",
     xlab = c("Log-Degree"), ylab = c("Log-Proportion"), 
     main="Log-Log Degree Distribution")

