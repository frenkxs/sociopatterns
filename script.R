rm(list=ls())
load('dataset.RData')


library(igraph)

x <- graph_from_adjacency_matrix(X, mode = 'undirected')
y <- graph_from_adjacency_matrix(Y, weighted = TRUE, mode = 'undirected')


list.edge.attributes(x)
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
pr <- page_rank(x)$vector 

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

# Exponential Random Graph model

library(Bergm)
library(intergraph)
library(network)


# subset of the graph with only the nodes interacting the time 3
x_three <- induced_subgraph(x, V(x)[time_frame == 3])
x_three <- asNetwork(x_three)

x_ergm <- x_three ~ edges + kstar(2) + triangle

summary(x_ergm)

x_bergm <- bergm(x_ergm)
summary(x_bergm)
plot(x_bergm)
plot(bgof(x_bergm))


#--------------------------------------------

# Euclidean Latent Position Model
library(latentnet)

# one (no) clusters
x_lpm1 <- ergmm(x_three ~ euclidean(d = 2, G = 1), verbose = TRUE)

# two clusters
x_lpm2 <- ergmm(x_three ~ euclidean(d = 2, G = 2), verbose = TRUE)

# three clusters
x_lpm3 <- ergmm(x_three ~ euclidean(d = 2, G = 3), verbose = TRUE)

# four clusters
x_lpm4 <- ergmm(x_three ~ euclidean(d = 2, G = 4), verbose = TRUE)

summary(x_lpm1)
extract(x_lpm1, include.bic = TRUE)

bic.ergmm(x_lpm1)$overall
bic.ergmm(x_lpm2)$overall 
bic.ergmm(x_lpm3)$overall
bic.ergmm(x_lpm4)$overall

# model with 2 cluster is the best

G <- 2

plot(x_lpm2, main = "Latent positions model with 2 clusters", print.formula = FALSE, 
     labels = TRUE, label.cex = 0.5, label.pad = 0.1)


#--------------------------------------------

# resilience - random

# Number of iteration
N <- 50 

# Number of nodes to be removed
remove <- 200 

# Initialise containers to store the statistics
component_s <- matrix(NA, N, remove + 1) 
clustering_c <- matrix(NA, N, remove + 1) 
modularity_c <- matrix(NA, N, remove + 1)

# add the initial values
component_s[, 1] <- components(x)$csize
clustering_c[, 1] <- transitivity(x)
modularity_c[, 1] <- modularity(x, time_frame)

for (iter in 1:N) {
  # initiate the network
  network <- x
  
  for (i in 1:remove) {
    # choose a node at random to be removed 
    selected_n <- sample(1:length(V(network)), 1, replace = TRUE) 
     
    # remove the selected node
    network <- delete.vertices(network, selected_n)
    
    # calculate largest component
    component_s[iter, i + 1] <- max(components(network)$csize) 
    
    # calculate clustering coefficient
    clustering_c[iter, i + 1] <- transitivity(network)
    
    # calculate modularity
    modularity_c[iter, i + 1] <- modularity(network, time_frame[-selected_n])
  }
}

### Plotting function

plot_net <- function(stats, ci = FALSE) {
  # calculate the median and quantiles 
  dat <- apply(stats, 2, quantile, probs = c(0.05, 0.5, 0.95)) 
  
  # create the background for the plot
  plot(dat[2,], type = "n", xlab = "nodes removed", ylab = "statistic") 
  
  # add the empirical confidence bounds
  if (ci) {
  polygon(c(1:(remove + 1),(remove + 1):1), c(dat[1, ], rev(dat[3, ])), 
          col = "grey", border = NA) 
  }
  # plot the median as a line
  lines(dat[2,], type = "l", lwd = 2) 
}

plot_net(clustering_c)
plot_net(component_s, ci = TRUE)
plot_net(modularity_c, ci = TRUE)

#--------------------------------------------

# resilience - targeted removal

# Number of nodes to be removed
remove <- 400 

# Initialise containers to store the statistics
s_degree <- s_centrality <- s_pagerank <- rep(NA, N, remove + 1) 

# add the initial values
s_degree[1] <- s_pagerank[1] <- s_centrality[1] <- components(x)$csize

# initialise three networks for simulation
x_degree <- x_centrality <- x_pagerank <- x

for (i in 1:remove) {
  # select nodes to remove 
  node_degree <- which.max(degree(x_degree))
  node_centrality <- which.max(centr_eigen(x_centrality)$vector)
  node_pagerank <- which.max(page_rank(x_pagerank)$vector)
  
  # remove the selected nodes
  x_degree <- delete.vertices(x_degree, node_degree) 
  x_centrality <- delete.vertices(x_centrality, node_centrality)
  x_pagerank <- delete.vertices(x_pagerank, node_pagerank)
    
  # calculate largest component
  s_degree[i + 1] <- max(components(x_degree)$csize) 
  s_centrality[i + 1] <- max(components(x_centrality)$csize) 
  s_pagerank[i + 1] <- max(components(x_pagerank)$csize) 
}

# plot results
plot(c(0:400), s_degree, type = "l", xlab = "Number of nodes removed", 
     ylab = "Size of the largest component", main = "By degree", lwd = 2) 
plot(c(0:400), s_centrality, type = "l", xlab = "Number of nodes removed", 
     ylab = "Size of the largest component", main = "By centrality", lwd = 2) 
plot(c(0:400), s_pagerank, type = "l", xlab = "Number of nodes removed", 
     ylab = "Size of the largest component", main = "By page rank", lwd = 2)
