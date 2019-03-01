#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are 2 arguments if not, return an error with required fields
if (length(args) != 3) {
  stop("Incorrect arguments. Please input (input file name, initial seeding file, output file name).", call.=FALSE)
}  

## put arguments into required forms
input = as.data.frame(read.table(args[1], header=TRUE))
seed =  as.data.frame(read.table(args[2], header=FALSE))
output= args[3]

#flip the matrix for convenience
input<- t(input)
#merge data with custer label
input <- cbind(seed, input)
#relabel the column for clarity
colnames(input)[1] <- "Cluster"

#initialize vectors for loop
#just use 1 and 2 to initialize different values for while loop
old_clust <- 1
new_clust <- 2
iteration <- 0
centroid <- c()
metric <- c()
#loop through k-means algorithm until old clusters and new clusters are the same value
while (identical(old_clust, new_clust) == FALSE){
  iteration <- iteration + 1
  print(iteration)
  #old cluster is the cluster assignments at beginning of iteration
  old_clust <- input$Cluster
  #calculate centroids for each cluster
  for (i in 1:7){
     centroid[[i]] <- colMeans(input[input$Cluster==i,-1])
  }
  
  dist <- c() #initialize euclidean distance vector
  tot_dist<- c()
  for (j in 1:nrow(input)){ #iterate through each row
     for (i in 1:7){ #iterate through each centroid
       #calculate euclidean distance between sample and each centroid
        dist[i]<- sqrt(sum((input[j,-1] - unlist(centroid[i]))^2)) 
     }
  new_clust[j] <- which.min(dist) #assign the sample to the centroid with smallest euclidean distance
  tot_dist[j]<- min(dist) #obtain vector of distances between data point and closest centroid
  }
  input$Cluster <- new_clust #replace old cluster assignments with new cluster assignments
  metric[iteration] <- mean(tot_dist) #calculate objective function 

}

x<- c(as.character(input$Cluster), metric[max(iteration)])

write(x, file= output)

    