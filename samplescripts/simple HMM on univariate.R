library(mhsmm)
library(readr)

#declaring formatMhsmm
formatMhsmm <- function(data){
  
  nb.sequences = nrow(data)
  nb.observations = length(data)
  
  #transform list to data frame
  data_df <- data.frame(matrix(unlist(data), nrow = nb.sequences, byrow=F))
  
  
  #iterate over these in loops
  rows <- 1:nb.sequences
  observations <- 0:(nb.observations-1)
  
  #build vector with id values
  id = numeric(length = nb.sequences*nb.observations ) 
  
  for(i in rows)
  {
    for (j in observations)
    {
      id[i+j+(i-1)*(nb.observations-1)] = i
    }
  }
  
  #build vector with observation values
  sequences = numeric(length = nb.sequences*nb.observations) 
  
  for(i in rows)
  {
    for (j in observations)
    {
      sequences[i+j+(i-1)*(nb.observations-1)] = data_df[i,j+1]
    }
  }
  
  data.df = data.frame(id, sequences)
  
  #creation of hsmm.data object needed for training
  N <- as.numeric(table(data.df$id))
  train <- list(x = data.df$sequences, N = N)
  class(train) <- "hsmm.data"
  
  return(train)
}

#loading data
dataset <- read_csv("YOUR LOCAL DIRECTORY HERE",header=TRUE)

train <- formatMhsmm(data.frame(dataset$Global_active_power))
#day subset

traindayform <- formatMhsmm(data.frame(train$x[1:5000]))

#testdayform <- formatMhsmm(data.frame(test$Global_active_power[1:1000]))

#end of data
# 4 states HMM    
k=6
#init probabilities
init <- rep(1/k, k)

#transition matrix
P <- matrix(rep(1/k, k*k), nrow = k)

#emission matrix:  here I used a Gaussian distribution, replace muEst and sigmaEst by your initial estimates of mean and variance
b <- list(mu = c(1,4), sigma = c(2,1)) 

#starting model for EM
startmodel <- hmmspec(init = init, trans = P, parms.emis = b, dens.emis = dnorm.hsmm)
startmodel
#EM algorithm fits an HMM to the data
hmm <- hmmfit(traindayform$x, startmodel , mstep = mstep.norm,maxit = 200)

#print resulting HMM parameters
summary(hmm)
plot(hmm$loglik, type="b", ylab="log-likelihood", xlab="Iteration")

#testhmm1 <- formatMhsmm(test)
yhat1 <- predict (hmm,traindayform$x)
#yhat2 <- predict (hmm,testdayform$x)
#modelbased <- predict(startmodel,train,method="smoothed" )
#plot(modelbased)

plot(yhat1)
#addstates(yhat1$s)
#plot(yhat2)
#addstates(yhat2$s)
