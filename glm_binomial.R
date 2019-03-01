#!/usr/bin/env Rscript

#Forward selection algorithm for logistic regression with input single cell RNA seq data


args = commandArgs(trailingOnly=TRUE)

# test if there are 2 arguments if not, return an error with required fields
if (length(args) != 2) {
  stop("Incorrect arguments. Please input (input file name, output file name).", call.=FALSE)
}  

## put arguments into required forms
input = as.data.frame(read.table(args[1], header=TRUE))
output= args[2]

#flip matrix to simplify coding
z<-t(input)

#isolate the gene names
gene_names<- colnames(z)

#remove the numbers from the cell names to get Memory or Naive
sub <- gsub("_.*" , "", rownames(z))
#add these names to the data frame
cell_type<- as.data.frame(sub)
y<- cbind(z, cell_type)

#assign a 0 to Memory cells and a 1 to Naive cells
cell_numb<- c()
for( i in 1:nrow(y)){
  if (y[i,ncol(y)] == "Memory"){
    cell_numb[i]<- 0
  } else {
    cell_numb[i]<-1
  }
}
#add this to the data frame
y<-cbind(y,cell_numb)


#set empty vectors for loop to store output
RSS <-c()
Iteration = c()
Predictor = c()
Log_likelihood = c()
CV= c()

#set the null model
formula <- glm(cell_numb ~ 1, family = binomial(logit), data=y)

#iterate through 50 times, i.e. add 50 genes to the model
for (i in 1:50){ 
  Iteration[i] <- i #store iteration number
  Predictor[i] <- gene_names[which.min(RSS)] #store added gene name
  q<- predict(formula, y, type = "response") #calculate predicted probability
  q<- as.numeric(gsub("/t*", "", q)) #remove text
  Log_likelihood[i]<- log(prod(q)*prod(1-q)) #calculate log-likelihood
  # CV error calculation ------
  MSEperFold = c()
  # stop the CV fold sampling while k folds is less than 5
  while (length(MSEperFold) < 5) {
    # initialize vector of MSEs per fold
    all_rows <- nrow(y)
    AllRows <- 1:all_rows
    
    # select 4/5ths of the data to train GLM on
    train_rows = sample(x = 1:all_rows, size = (0.8 * all_rows), replace = F)
    
    # select 1/5ths of the data to test GLM on
    test_rows = AllRows[-c(train_rows)]
    
    train_rows %in% test_rows
    AllRows %in% train_rows
    AllRows %in% test_rows
    
    # subset the train and test data sets
    train_data = y[train_rows, ]
    test_data = y[test_rows, ]
    
    GLMforCV = multinom( as.formula(paste(formula,'+', Predictor[i])),
                                    data = train_data, trace = F)
    
    
    # how off were the predictions using our best current model
    PredictedProb = predict(object = GLMforCV, test_data, type = 'probs')
    PredictedProbab
    # the output from predict is the likelihoods of your model predicting what you 
    #previously labeled as 1, whether it be memory or naive
    # in my case it is memory cells
    
    # truth labels for test class for this fold
    testclass<- c()
    for (i in 1:length(levels)){
      t<- grep(levels[i], rownames(PredictedProb))
      item <- levels[i]
      for(i in 1:length(t)){
        testclass[t[i]]<- item
      }
    }
    
    # creating matrix for observed prob
    ObsForTestSet = matrix(0, nrow=nrow(PredictedProbab), ncol = length(levels))
    rownames(ObsForTestSet)= rownames(PredictedProb)
    colnames(ObsForTestSet)= levels
    for ( i in 1:ncol(ObsForTestSet)){
      truthlabels <- colnames(ObsForTestSet)[i]
      for (i in 1:length(testclass)){
        if (testclass[i] %in% truthlabels){
          col<- grep(truthlabels, colnames(ObsForTestSet))
          ObsForTestSet[i,col] <- 1
        }
      }
    }
    
    # substract the two to perform ( sum(yhati - y)**2 ) / n = MSE
    MSEperFold = append(MSEperFold, ( ( sum((PredictedProb - ObsForTestSet)**2) 
                                        / length(PredictedProb)) ) )
  }
  
  CV[i] <- ( sum(MSEperFold) / 5 ) 

  for( i in 1:length(gene_names)){ #forward selection
    g <- update(formula, paste("~ . +",gene_names[i])) #add each gene to the formula
    RSS[i]<- sum((y$cell_numb-g$fitted.values)^2) #calculate and store RSS for each gene
  }
  formula<- update(formula, paste("~ . +",gene_names[which.min(RSS)])) #update formula with gene with smallest RSS
  }



#create a data frame with desired output
out <- cbind(Iteration, Predictor, CV, Log_likelihood)

#write output to file
write.table(out, file= output, row.names = FALSE)



