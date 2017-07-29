library(mhsmm)
library(readr)

#--------------Function Declarations-------------------------------
# declaring formatMhsmm from https://stackoverflow.com/questions/21682619/mhsmm-package-in-r-input-format
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

detectPointAnomalies <- function(hmm, testData, threshold, testOnValidation, corrupt) {
  anomalyPointCount = 0
  truePositiveCount = 0
  yhat <- predict(hmm, testData)
  for (i in 1:length(yhat$s)) {
    observationState <- yhat$s[i]
    observationStateMean <- hmm$model$parms.emission$mu[observationState]
    #This is the same as the test data point
    observationDataPoint <- yhat$x[i]
    if (abs(observationDataPoint - observationStateMean) > threshold) {
      anomalyPointCount = anomalyPointCount + 1
      if (testOnValidation) {
        if (corrupt[i]) {
          truePositiveCount = truePositiveCount + 1
        }
      }
    }
  }
  return(list(anomalyPointCount, truePositiveCount))
}

# Snippet taken from here: https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
splitVector <- function(vector, sequenceSize) {
  return(split(vector, ceiling(seq_along(vector)/sequenceSize)))
}





#-------------Data and HMM Initialization---------------------------------
# load Training Data (only using Global Active Power), 
# delete entries that are not available, 
# and format the data for input to mhsmm library
trainDataset <- read_csv("C:\\Users\\Trevor\\Google Drive\\School\\SFU\\Year 4\\Summer 2017\\CMPT 318\\Project\\CMPT318\\data\\train.txt")
# trainDataset <- read_csv("C:\\Users\\Evan Chisholm\\Desktop\\CMPT318\\train.txt")
trainGlobalActivePower <- trainDataset$Global_active_power
trainGlobalActivePower <- trainGlobalActivePower[!is.na(trainGlobalActivePower)]

testOnValidation <- TRUE
#------ Creating Validation Set from Training Set ------------------------
if(testOnValidation){
  #Take the last 10% of the dataset to use as a validation set
  lastTenPercentRange <- (0.9*length(trainGlobalActivePower)) : length(trainGlobalActivePower)
  validationGlobalActivePower <- trainGlobalActivePower[lastTenPercentRange]
  
  #truncate trainGlobalActivePower as we do not want to train on any data that is in the validation set
  trainRange <- 0.9*length(trainGlobalActivePower)
  trainGlobalActivePower <- trainGlobalActivePower[1:trainRange]
  
  #add noise to our validationGlobalActivePower and keep track of where noise was inserted for scoring purposes
  corrupt <- rbinom(length(validationGlobalActivePower), 1, 0.03) #determine an average of 3% of the data to replace with noise
  #we can reference this vector later to determing which values should have been detected as noise
  corrupt <- as.logical(corrupt)
  noise <- rnorm(sum(corrupt), 15, 3) #creates our noise values to insert by the normal distribution with a mean on 15, sd of 3, which should be pretty obvious
  validationGlobalActivePower[corrupt] <- validationGlobalActivePower[corrupt] + noise
  testFormattedData <- formatMhsmm(data.frame(validationGlobalActivePower))
}

trainFormattedData <- formatMhsmm(data.frame(trainGlobalActivePower))


# do the same thing for Test data
if(!testOnValidation){
  testDataset <- read_csv("C:\\Users\\Trevor\\Google Drive\\School\\SFU\\Year 4\\Summer 2017\\CMPT 318\\Project\\CMPT318\\data\\test1.txt")
  # testDataset <- read_csv("C:\\Users\\Evan Chisholm\\Desktop\\CMPT318\\train.txt")
  testGlobalActivePower <- testDataset$Global_active_power
  testGlobalActivePower <- testGlobalActivePower[!is.na(testGlobalActivePower)]
  testFormattedData <- formatMhsmm(data.frame(testGlobalActivePower))
}
# Specify initial HMM parameter values
# Determined from Evan and Heather's Analysis
numStates <- 5
initialStateProbabilities <- rep(1/numStates, numStates)
transitionMatrix <- matrix(rep(1/numStates, numStates*numStates), nrow = numStates)

#using k-means clustering, estimate initial mu and sigma values for the emission matrix
kgap <- kmeans(trainGlobalActivePower,5,iter.max=10)
muinit <- kgap$centers
varinit = numeric(numStates)
for (i in 1:numStates){
  varinit[i] <- kgap$withinss[i]/(kgap$size[i]-1)
}
emissionMatrix <- list(mu = muinit, sigma = varinit) 

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





#-------------Collective Anomaly Detection Approach Based On Log-Likelihood---------------------
# divide training data into sequences that are each a given length.
# Note: the last sequence may not be a full sequenceSize in length
sequenceSize <- 5
trainSequences <- splitVector(trainFormattedData$x, sequenceSize)

# the last sequence in trainSequences may be less than sequenceSize (because the data points don't divide perfectly by sequenceSize)
# This ensures such a sequence won't be tested
# There are probably better ways to deal with this issue, but this will do for now.
validTrainSequencesCount = if (length(trainSequences[[length(trainSequences)]]) == sequenceSize) length(trainSequences) else (length(trainSequences) - 1)
trainLogLikelihoods <- numeric(validTrainSequencesCount)
# Calculate log-likelihoods of train sequences
for (i in 1:validTrainSequencesCount) {
  sequence <- trainSequences[[i]]
  yhat <- predict (hmm, sequence)
  trainLogLikelihoods[i] <- yhat$loglik
}

# Use the 1.5 x InterQuartile Range Rule for Outliers to create a range of normal log-likelihood
# Description of the rule: http://www.purplemath.com/modules/boxwhisk3.htm
# Note: This is probably not a very good way to do this, because often the range of normal 
# log-likelihood ends up exceeding the minimum or maximum log-likelihood values
# Ideally, we would use a validation data set to analyse how effective our range of normal log-likelihood is,
# and adjust our range until it effectively captures anomalies.
# quartiles <- fivenum(trainLogLikelihoods)
# quartile1Index = 2
# quartile3Index = 4
# interQuartileRange = quartiles[quartile3Index] - quartiles[quartile1Index]
# minLogLikelihood <- quartiles[quartile1Index] - (1.5 * interQuartileRange)
# maxLogLikelihood <- quartiles[quartile3Index] + (1.5 * interQuartileRange)
rangeWidth = 16
trainLogLikelihoodsMedian <- median(trainLogLikelihoods)
minLogLikelihood <- trainLogLikelihoodsMedian - rangeWidth
maxLogLikelihood <- trainLogLikelihoodsMedian + rangeWidth

# Divide test data into sequences that are each a sequenceSize in length
# IMPORTANT: train and test sequences must be the same length (unless their log-likelihoods are scaled)
testSequences <- splitVector(testFormattedData$x, sequenceSize)

if (testOnValidation) {
  corruptSequences <- splitVector(corrupt, sequenceSize)
}

collectiveAnomalyThreshold <- 1
# there must be at least this many point anomalies in a given collective anomaly for it to be considered valid
validCollectiveAnomalyRatio <- .2
truePositiveCollectiveAnomalyCount <- 0
actualCollectiveAnomalyCount <- 0
# Identify Anomalies
# Calculate log-likelihood of test sequences. Any sequence outside the normal range is an anomaly 
anomalyCollectiveCount <- 0
for (i in 1:length(testSequences)) {
  sequence <- testSequences[[i]]
  # similar to before, this ensures the last sequence (which may not be a sequenceSize in length) won't be tested
  if (length(sequence) == sequenceSize) {
    yhat <- predict (hmm, sequence)
    isAnomaly <- yhat$loglik > maxLogLikelihood || yhat$loglik < minLogLikelihood
    if (isAnomaly) {
      anomalyCollectiveCount = anomalyCollectiveCount + 1 
    }
    if (testOnValidation) {
      # TODO: change corrupt to subsequences
      returnValue <- detectPointAnomalies(hmm, sequence, collectiveAnomalyThreshold, testOnValidation, corruptSequences[[i]])
      truePointAnomaliesInSequence <- returnValue[[2]]
      # This differs from isAnomaly, because it checks if a sequence SHOULD have been flagged by the algorithm as an anomaly or not
      isActualAnomaly <- (truePointAnomaliesInSequence / sequenceSize) >= validCollectiveAnomalyRatio
      if (isActualAnomaly) {
        if (isAnomaly) {
          truePositiveCollectiveAnomalyCount = truePositiveCollectiveAnomalyCount + 1  
        }
        actualCollectiveAnomalyCount = actualCollectiveAnomalyCount + 1
      }
    }
  }
}

# Report Results
# similar to before, the last sequence may not be sequenceSize in length, so we may disclude it
validTestSequenceCount <- if (length(testSequences[[length(testSequences)]]) == sequenceSize) length(testSequences) else (length(testSequences) - 1) 
anomalyCollectivePercent <- 100 * (anomalyCollectiveCount / validTestSequenceCount)
cat("Collective Anomaly Percentage: ", anomalyCollectivePercent, "%")
cat("Collective Anomaly Count: ", anomalyCollectiveCount)
if (testOnValidation) {
  precision <- 100 * (truePositiveCollectiveAnomalyCount / anomalyCollectiveCount)
  actualAnomalyCount <- sum(corrupt)
  recall <- 100 * (truePositiveCollectiveAnomalyCount / actualCollectiveAnomalyCount)
  cat("Precision: ", precision, "%\n")
  cat("Recall: ", recall, "%\n")
}





#-----------Point Anomaly Detection------------------------------
# The general idea here is to look at each point individually, and compare it to an expected value.
# Here, the expected value is determined by the most probable state sequence. 
# You look at the most probable state assigned to a data point (in predict()), and find that state's output 
# emission mean (specified in the HMM). This mean represents an expected normal value for the given state
threshold = 2
returnValue <- detectPointAnomalies(hmm, testFormattedData$x, threshold, testOnValidation, corrupt)
anomalyPointCount <- returnValue[[1]]
truePositiveCount <- returnValue[[2]]
# Report Results
anomalyPointPercentage = 100 * (anomalyPointCount/length(testFormattedData$x))
cat("Point Anomaly Percentage: ", anomalyPointPercentage, "%")
cat("Point Anomaly Count: ", anomalyPointCount)
if (testOnValidation) {
  precision <- 100 * (truePositiveCount / anomalyPointCount)
  actualAnomalyCount <- sum(corrupt)
  recall <- 100 * (truePositiveCount / actualAnomalyCount)
  cat("Precision: ", precision, "%\n")
  cat("Recall: ", recall, "%\n")
}

