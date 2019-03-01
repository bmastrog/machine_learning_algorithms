#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are 2 arguments if not, return an error with required fields
if (length(args) != 2) {
  stop("Incorrect arguments. Please input (input file name, output file name).", call.=FALSE)
}  

## put arguments into required forms
input = as.data.frame(read.table(args[1], header=TRUE))
output= args[2]

#load required packages
library(randomizr)
#flip the matrix for convenience
input<- t(input)

save_clust <- c()
metric <- c()
for (p in 1:25){
  print(p)
  Cluster <- complete_ra(N = nrow(input), m = NULL, m_each = c(20,20,20,20,20,20,20), prob = NULL, prob_each = NULL,
                         num_arms = 7, condition_names = c(1, 2, 3, 4, 5, 6, 7), check_inputs = TRUE)
  #old cluster is the cluster assignments at beginning of iteration
  #input<- as.data.frame(cbind(Cluster, input))
  #initialize vectors for loop
  #just use 1 and 2 to initialize different values for while loop
  old_clust <- 1
  new_clust <- 2
  iteration <- 0
  centroid <- c()
  #loop through k-means algorithm until old clusters and new clusters are the same value
  while (identical(old_clust, new_clust) == FALSE){
    iteration <- iteration + 1
    #old cluster is the cluster assignments at beginning of iteration
    old_clust <- Cluster
    #calculate centroids for each cluster
    for (i in 1:7){
      centroid[[i]] <- colMeans(input[Cluster==i, , drop=FALSE])
    }
    
    dist <- c() #initialize euclidean distance vector
    tot_dist<- c()
    for (j in 1:nrow(input)){ #iterate through each row
      for (i in 1:7){ #iterate through each centroid
        dist[i]<- sqrt(sum((input[j, ] - unlist(centroid[i]))^2)) #calculate euclidean distance between sample and each centroid
      }
      new_clust[j] <- which.min(dist) #assign the sample to the centroid with smallest euclidean distance
      tot_dist[j]<- min(dist, na.rm= TRUE)
    }
    Cluster <- new_clust #replace old cluster assignments with new cluster assignments
    metric[p] <- mean(tot_dist, na.rm = TRUE) 
    # need to add calculation for total minimization equation of kmeans.  want to store each run to ensure getting smaller
    
  }
  save_clust[[p]]<- Cluster
  #input<- input[,-1]
}

x<- c(as.character(save_clust[[which.min(metric)]]), metric[which.min(metric)])

write(x, file= output)
