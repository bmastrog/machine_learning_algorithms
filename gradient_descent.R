#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there are 5 arguments if not, return an error with required fields
if (length(args) != 5) {
  stop("Incorrect arguments. Please input (response file, predictor file, alpha, epsilon, outfile).", call.=FALSE)
}  


## put arguments into required forms
response = as.matrix(read.table(args[1], header= FALSE))
predictor = as.matrix(read.table(args[2], header=TRUE))
alpha = as.numeric(args[3])
epsilon = as.numeric(args[4])
outfile= args[5]

#mark as x and y to simplify code
x<- as.matrix(predictor)
y<- as.matrix(response)

#Add ones to x to account for intercept
x <- cbind(1,x)

# initalize theta vector
theta <- as.vector(rep(0, ncol(x)))

# Number of the observations
m <- nrow(x)


#function for linear model
linear= function(x,th){
  (x %*% th)}

#function for cost calculation
cost = function(x,y,th){
  sum((x %*% th - y)^2) / length(y)}

#initial values for gradient descent
converged= FALSE
previousCost= 0.0
#gradient descent 
#will repeat until previouscost-newcost < epsilon
while( converged == FALSE){
  for(j in 2:ncol(x)){
    temp1 <- theta[1] - alpha * (1/m) * sum(((x%*%theta)- y)) 
    temp <- theta[j] - alpha * (1/m) * sum(((x%*%theta)- y)*x[,j])
    theta[1]= temp1
    theta[j] = temp
  cost1 = cost(x,y,theta)
  if(abs(cost1- previousCost) <= epsilon){
    converged= TRUE
  }else{
    previousCost= cost1
    }
  }
}

#calculate predicted y values
z <-linear(x,theta)

#write values to file
write(z, file= outfile, ncolumns = 1)
