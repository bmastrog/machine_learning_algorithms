#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are 2 arguments if not, return an error with required fields
if (length(args) != 3) {
  stop("Incorrect arguments. Please input (input file name, random vector, output file name).", call.=FALSE)
}  

## put arguments into required forms
input = as.data.frame(read.table(args[1], header=TRUE))
vector =  as.data.frame(read.table(args[2], header=FALSE))
output= args[3]

#input <- read.table("input_forPCA.txt", header = TRUE)
#vector <- read.table("random_vector.txt", header= FALSE)

#center the data
cent_dat<- scale(input, scale=FALSE, center=TRUE)

#transpose the vector values for convenience
vector <- t(vector)

#set up a function that normalizes parameters to have squares sum to 1
param = function(pars){
  pars^2/sum(pars^2)
}

#create a function that projects samples onto a line and calculates variance
my.var=function(thetas, x.vals) {
  z<- c()
  y<- param(thetas)
  projections <- matrix( , ncol = 7, nrow = 140)
  project <- c()
  for (i in 1:nrow(x.vals)){
    project <- sum(as.numeric(x.vals[i,]) * y)/sum(y*y) 
    projections[i,]<- project * y 
  }
  variance <- mean(colSums(projections)^2)
  return(variance)
}

#initialize vectors
new_dat<- cent_dat
variance <-c()
for (z in 1:3){ #calculate the first 3 principal components 
   best.params=optim(vector,fn = my.var, x.vals=new_dat,control= list(fnscale=-1)) #find optimum value
   variance[z] <- best.params$value #store variance
   new_vec<- best.params$par #obtain new vector
   new_project <- matrix( , ncol = 7, nrow = 140) #set up empty matrix to store new projections
   for (i in 1:nrow(input)){ #calculate new projected data
     new_project[i,] <- (sum(as.numeric(input[i,]) * new_vec)/sum(new_vec*new_vec))*new_vec 
   }
   new_dat<- new_dat - new_project #subtract new points from old points to uncorrelate for calculation of next PC 
} 
  
variance

#write varaince to file 
write(variance, file = output)
