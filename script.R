rm(list=ls())
load('dataset.RData')


library(igraph)

x <- graph_from_adjacency_matrix(X, mode = 'undirected')
y <- graph_from_adjacency_matrix(Y, weighted = TRUE, mode = 'undirected')



list.edge.attributes(y)
E(y)$weight

# number of nodes
n_nodes <- length(V(x))
n_edges <- length(E(x))

# visualise the static and weigthed network
igraph.options(vertex.label = NA,
               vertex.size = 2, 
               edge.width = 0.5, 
               edge.color = "grey50", 
               edge.curved = 0.1)

plot(x, layout = layout)
plot(y, layout = layout, edge.label = E(y)$weight, edge.label.cex = 0.5)

#--------------------------------------------

# Binary degree of X and the weighted degrees of Y. 
d_bin <- degree(x) 
d_weight <- strength(y)

# 5 nodes with the highest binary degree and their weighted degrees.
best_bin <- order(d_bin, decreasing = TRUE)[1:5]
d_weight[best_bin]

#--------------------------------------------

# Visual representation of the network X, where the size of the nodes is proportional to the weighted degree, and the nodes are coloured according to time_frame.

igraph.options(vertex.label = NA,
               vertex.size = d_weight * 0.05, 
               vertex.color = time_frame,
               edge.width = 0.5, 
               edge.color = "grey50", 
               edge.curved = 0.1)

plot(x, layout = layout)

#--------------------------------------------

# degree distribution: natural scale 
hist(d_bin, col="blue", 
     xlab="Degree", ylab="Frequency", 
     main="Degree Distribution")

# log-log scale 
d_dist <- degree.distribution(x)
d <- 0:max(d_bin) 
ind <- (d_dist != 0)
t <- d[ind]
u <- d_dist[ind]
barplot(log(d[ind]), log(d_dist[ind]), col = "grey50",
     xlab = c("Log-Degree"), ylab = c("Log-Proportion"), 
     main="Log-Log Degree Distribution")


#--------------------------------------------

# average degree
mean(d_bin)

# clustering coefficient
transitivity(x)

# average path length
average.path.length(x)

# average degree of the neighbours
knn(x)

average_deg <- function (graph) {
  n <- length(V(graph))
  dnn <- rep(NA, n)
  for (i in 1:n) {
    dnn[i] <- mean(degree(graph, neighbors(graph, i)))
  }
  return(dnn)
}

average_deg(x) == knn(x)$knn

plot(d_bin, knn(x)$knn)

cor(d_bin, knn(x)$knn)

# Positive correlation indicate assortative mixing > nodes tend to connect with other nodes that have similar degrees

#--------------------------------------------

# Page rank centrality
pr <- page_rank(x)$vector # Page-Rank centrality 

# Betweenness centrality 
bc <- centr_betw(x)$res 

top_pr <- order(pr, decreasing = TRUE)[1:10]
top_bc <- order(bc, decreasing = TRUE)[1:10]

central <- rep(1, length(V(x)))
top <- intersect(top_pr, top_bc)
central[top] <- 2
  

igraph.options(vertex.label = NA,
               vertex.size = pr * 1000, 
               vertex.color = central,
               edge.width = 0.5, 
               edge.color = "grey50", 
               edge.curved = 0.1)

plot(x, layout = layout)

#--------------------------------------------

# spectral clustering with Laplacian matrix
L <- graph.laplacian(x) 
eigen_decomposition <- eigen(L) 

# number of clusters
n_class <- 3 

# select the correct number of eigenvectors
dimensions <- n_nodes:(n_nodes - n_class + 1) 
embedding <- eigen_decomposition$vectors[, dimensions]

# kmeans on the set of selected eigenvectors
members <- kmeans(embedding, n_dim, nstart = 100)$cluster 

# convert to clusters
clust_spect <- make_clusters(x, members) 

# plot
plot(clust_spect, x, layout = layout, vertex.color = time_frame)

# inspect the number of missclassified nodes
table(time_frame, members)

#--------------------------------------------

#Stochastic Block Model on X with class labels corresponding to time_frame. 

# create separate adjacency matrix for each cluster
X_1 <- X[time_frame == 1, time_frame == 1]
X_2 <- X[time_frame == 2, time_frame == 2]
X_3 <- X[time_frame == 3, time_frame == 3]


# compute the MLE, ie. proportion of existing edges to all possible edges.
sum(X_1) / (ncol(X_1) ^ 2)
sum(X_2) / (ncol(X_2) ^ 2)
sum(X_3) / (ncol(X_3) ^ 2)

#--------------------------------------------

library(blockmodels)

sbm_x <- BM_bernoulli(membership_type = "SBM_sym", 
                    adj = X,
                    verbosity = 1, # how much should be printed out while running the algorithm
                    plotting = "",
                    explore_min = 2,
                    explore_max = 8) # maximum number of clusters to consider

# perform VEM
sbm_x$estimate()

# extract the best model according to ICL values
best <- which.max(sbm_x$ICL) 

# exctract probabilities for class memberships and assing nodes to classes based on them
clust_prob <- sbm_x$memberships[[best]]$Z 
clust_members <- apply(clust_prob, 1, which.max) 

# plot group sizes
sbm_x$memberships[[best]]$plot()

# plot connection probabilities
sbm_x$plot_parameters(best) 

# plot the adjacency matrix
image(X[order(clust_members),rev(order(clust_members))], xaxt = "n", yaxt = "n") 

### highlight the blocks (very basic coding, it could be improved)
group_counts <- (as.numeric(table(clust_members)))
abline(v = cumsum(group_counts/sum(group_counts)))
abline(h = 1-cumsum(group_counts/sum(group_counts)))

# plot clusters
clust_sbm <- make_clusters(x, clust_members) 
plot(clust_sbm, x, layout = layout, vertex.label = time_frame)

table(time_frame, clust_members)

#--------------------------------------------

