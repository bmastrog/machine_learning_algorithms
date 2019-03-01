#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are 2 arguments if not, return an error with required fields
if (length(args) != 2) {
  stop("Incorrect arguments. Please input (input file name, output file name).", call.=FALSE)
}  

## put arguments into required forms
input = as.data.frame(read.table(args[1], header=TRUE))
output= args[2]

#flip the matrix for convenience
input<- t(input)

new_cent<- c()
#initial cluster
Cluster <- rep(1, length= nrow(input))
#initial centroid
new_cent[[1]]<- colMeans(input)

#input<- as.data.frame(cbind(Cluster, input))

metric <- c()
centroid<- c()
new_c <- c()
save_cent <- c()
save_clust <- c()
cluster<- c()
for (k in 2:4){
  print(k)
  metric <- c()
  centroid<- c()
  new_c <- c()
  save_cent <- c()
  save_clust <- c()
  for (t in 1:nrow(input)){ #iterate through each row
    new_c [[t]] <- input[t, ] #the row to be tested as a new centroid
    new_cent[[k]]<- new_c [[t]] #add new centroid to list of existing centroids
    centroid <- new_cent
    #initialize vectors for loop
    #just use 1 and 2 to initialize different values for while loop
    old_clust <- 1
    new_clust <- 2
    iteration <- 0
    #loop through k-means algorithm until old clusters and new clusters are the same value
    while (identical(old_clust, new_clust) == FALSE){
      iteration <- iteration + 1
      #old cluster is the cluster assignments at beginning of iteration
      old_clust <- cluster
      #calculate centroids for each cluster
      if (iteration==1){
        centroid <- new_cent
      }
      else{
        for (i in 1:k) {
          centroid[[i]] <- colMeans(input[cluster==i, ,drop=FALSE])
        }
      }
      
      dist <- c() #initialize euclidean distance vector
      tot_dist<- c()
      for (j in 1:nrow(input)){ #iterate through each row
        for (i in 1:k){ #iterate through each centroid
          #calculate euclidean distance between sample and each centroid
          dist[i]<- sqrt(sum((input[j,] - unlist(centroid[[i]]))^2)) 
        }
        new_clust[j] <- which.min(dist) #assign the sample to the centroid with smallest euclidean distance
        tot_dist[j]<- min(dist) #obtain vector of distances between data point and closest centroid
      }
      cluster <- new_clust #replace old cluster assignments with new cluster assignments
    }
    save_cent[[t]]<- centroid
    save_clust[[t]]<- cluster
    metric[t]<- mean(tot_dist, na.rm = TRUE)
  }
  
  new_cent<- save_cent[[which.min(metric)]]
}

save_clust[[which.min(metric)]]
metric[which.min(metric)]

x<- c(as.character(save_clust[[which.min(metric)]]), metric[which.min(metric)])

write(x, file= output)
