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

# load Training Data (only using Global Active Power), 
# delete entries that are not available, 
# and format the data for input to mhsmm library
trainDataset <- read_csv("C:\\Users\\Trevor\\Google Drive\\School\\SFU\\Year 4\\Summer 2017\\CMPT 318\\Project\\CMPT318\\data\\train.txt")
trainGlobalActivePower <- trainDataset$Global_active_power
trainGlobalActivePower <- trainGlobalActivePower[!is.na(trainGlobalActivePower)]
trainFormattedData <- formatMhsmm(data.frame(trainGlobalActivePower))

# do the same thing for Test data
testDataset <- read_csv("C:\\Users\\Trevor\\Google Drive\\School\\SFU\\Year 4\\Summer 2017\\CMPT 318\\Project\\CMPT318\\data\\test1.txt")
testGlobalActivePower <- testDataset$Global_active_power
testGlobalActivePower <- testGlobalActivePower[!is.na(testGlobalActivePower)]
testFormattedData <- formatMhsmm(data.frame(testGlobalActivePower))

# specify initial HMM parameter values
# values were determined from analysis of data set
# These are currently set to whatever values Evan was using in his analysis
numStates <- 5
initialStateProbabilities <- rep(1/numStates, numStates)
transitionMatrix <- matrix(rep(1/numStates, numStates*numStates), nrow = numStates)
emissionMatrix <- list(mu = c(0.85,3.426), sigma = c(2,1)) 

# create a starting model with initial parameters before performing Expectation Maximization algorithm (training)
startmodel <- hmmspec(init = initialStateProbabilities, 
                      trans = transitionMatrix, 
                      parms.emis = emissionMatrix, 
                      dens.emis = dnorm.hsmm)

# perform Expectation Maximization to obtain a trained HMM model
hmm <- hmmfit(trainFormattedData$x, startmodel, mstep = mstep.norm, maxit = 200)

# print resulting HMM parameters
summary(hmm)

# visualization of the Expectation Maximization process.
# The algorithm tries to maximize the log-likelihood, which you can see eventually plateaus after a number of iterations
# the iteration the log-likelihood plateaus is where we could stop the algorithm (i.e. set maxit of hmmfit() to that value)
plot(hmm$loglik, type="b", ylab="log-likelihood", xlab="Iteration")

# divide training data into sequences that are each a day in length.
# sequences could be any length. I'm not really sure how to choose the idea sequence length tbh
# Note: the last sequence may not be a full day in length
# Snippet taken from here: https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
minutesPerDay <- 1440
trainSequences <- split(trainFormattedData$x, ceiling(seq_along(trainFormattedData$x)/minutesPerDay))

# Create a range of normal log-likelihood using training data seqeunces defined by min/max loglikelihood variables
# Initialize min/max to dumb values that will never happen
for (sequence in trainSequences) {
  # the last sequence in trainSequences may be less than a day (because the data points don't divide perfectly into days)
  # This ensures such a sequence won't be tested
  # There are probably better ways to deal with this issue, but this will do for now.
  if (length(sequence) == minutesPerDay) {
    yhat <- predict (hmm, sequence)
    if (!exists("minLogLikelihood") || yhat$loglik < minLogLikelihood) {
      minLogLikelihood = yhat$loglik
    }
    if (!exists("maxLogLikelihood") || yhat$loglik > maxLogLikelihood) {
      maxLogLikelihood = yhat$loglik
    }
  }
}

# Divide test data into sequences that are each a day in length
# IMPORTANT: train and test sequences must be the same length (unless their log-likelihoods are scaled)
testSequences <- split(testFormattedData$x, ceiling(seq_along(testFormattedData$x)/minutesPerDay))

# Identify Anomalies
# Calculate log likelihood of test sequences. Any sequence outside the normal range is an anomaly 
for (sequence in testSequences) {
  # similar to before, this ensures the last sequence (which may not be a day in length) won't be tested
  if (length(sequence) == minutesPerDay) {
    yhat <- predict (hmm, sequence)
    if (yhat$loglik > maxLogLikelihood || yhat$loglik < minLogLikelihood) {
      print("ANOMALY DETECTEDDDDD")
    }
  }
}
