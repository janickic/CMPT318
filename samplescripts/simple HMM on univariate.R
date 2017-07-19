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
dataset <- read_csv("C:\Users\Evan Chisholm\Desktop\CMPT318\train.txt",header=TRUE)
gap <- dataset$Global_active_power
gap <- gap[!is.na(gap)]

minInDay <- 24*60
minInWeek <- 7*minInDay
minInMonth <- 30*minInDay #roughly
minInSeason <- 3*minInMonth 

train_start <- 2*minInSeason
train_end <- 3*minInSeason
test_start <- train_end + 1
test_end <- test_start + minInDay


train <- formatMhsmm(data.frame(gap[train_start:train_end]))
test <- formatMhsmm(data.frame(gap[test_start:test_end]))
#day subset

traindayform <- formatMhsmm(data.frame(train$x))

testdayform <- formatMhsmm(data.frame(test$x))

#end of data
# 5 states HMM    
k=5
#init probabilities
init <- rep(1/k, k)

#transition matrix
P <- matrix(rep(1/k, k*k), nrow = k)

#emission matrix:  here I used a Gaussian distribution, replace muEst and sigmaEst by your initial estimates of mean and variance
b <- list(mu = c(0.29, 0.62,1.46,2.56,4.31), sigma = c(2,1, 0.5, 0.25, 0.125)) 

#starting model for EM
startmodel <- hmmspec(init = init, trans = P, parms.emis = b, dens.emis = dnorm.hsmm)

#EM algorithm fits an HMM to the data
hmm <- hmmfit(traindayform$x, startmodel , mstep = mstep.norm,maxit = 200)

#print resulting HMM parameters
summary(hmm)
plot(hmm$loglik, type="b", ylab="log-likelihood", xlab="Iteration")

#testhmm1 <- formatMhsmm(test)
yhat1 <- predict (hmm,traindayform$x)
yhat2 <- predict (hmm,testdayform$x)
#modelbased <- predict(startmodel,train,method="smoothed" )
#plot(modelbased)
addStates(yhat1$s)
plot(yhat1)


#plot(yhat2)
#addStates(yhat2$s)

yhat1$loglik
yhat2$loglik
