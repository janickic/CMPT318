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

#declaring detectPointAnomalies for use in our anomaly point detection algorithm
detectPointAnomalies <- function(hmm, testData, threshold, testOnValidation, corrupt, shouldOutputResults) {
  anomalyPointCount = 0
  truePositiveCount = 0
  yhat <- predict(hmm, testData,method="smoothed")
  for (i in 1:length(yhat$s)) {
    observationState <- yhat$s[i]
    observationStateMean <- hmm$model$parms.emission$mu[observationState]
    #This is the same as the test data point
    observationDataPoint <- yhat$x[i]
    
    #We calculate the probability that our observationState matches the true state
    stateConfidence = 0 #p(observationState|observations), confidence that our stated state is actually the state
    for (prevState in 1:hmm$model$J){
      posterior = yhat$p[i-1,prevState] #the posterior distribution of the previous state
      transprob = hmm$model$transition[observationState,prevState] #possibility for state observationState given previous state
      stateConfidence = stateConfidence + posterior*transprob
    }
    
    #If the distance between our observation and expected value is greater than the threshold, it is flagged as an anomaly
    if (abs(observationDataPoint - observationStateMean) > threshold) {
      anomalyPointCount = anomalyPointCount + 1
      if (testOnValidation) {
        if (corrupt[i]) {
          truePositiveCount = truePositiveCount + 1
        }
      }
      if (shouldOutputResults) {
        cat("1,")
      }
    }
    else{
      if (shouldOutputResults) {
        cat("0,")
      }
    }
    if (i==1) {
      stateConfidence=0
    }
    if (shouldOutputResults) {
      cat(round(stateConfidence,digits=2),"\n")
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
#trainDataset <- read_csv("C:\\Users\\Trevor\\Google Drive\\School\\SFU\\Year 4\\Summer 2017\\CMPT 318\\Project\\CMPT318\\data\\train.txt")
# trainDataset <- read_csv("C:\\Users\\Evan Chisholm\\Desktop\\CMPT318\\train.txt")
trainDataset <-read_csv("/home/heather/Code/cmpt-318/train.txt")
trainGlobalActivePower <- trainDataset$Global_active_power
trainGlobalActivePower <- trainGlobalActivePower[!is.na(trainGlobalActivePower)]

trainFormattedData <- formatMhsmm(data.frame(trainGlobalActivePower))

testOnValidation <- FALSE
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

# do the same thing for Test data
if(!testOnValidation){
  #testDataset <- read_csv("C:\\Users\\Trevor\\Google Drive\\School\\SFU\\Year 4\\Summer 2017\\CMPT 318\\Project\\CMPT318\\data\\test1.txt")
  # testDataset <- read_csv("C:\\Users\\Evan Chisholm\\Desktop\\CMPT318\\train.txt")
  testDataset <-read_csv("/home/heather/Code/cmpt-318/test2.txt")
  #testDataset <-read_csv("/home/heather/Code/cmpt-318/test2.txt")
  testGlobalActivePower <- testDataset$Global_active_power
  testGlobalActivePower <- testGlobalActivePower[!is.na(testGlobalActivePower)]
  testFormattedData <- formatMhsmm(data.frame(testGlobalActivePower))
}


# Specify initial HMM parameter values
# Determined from Evan and Heather's Analysis
numStates <- 5
initialStateProbabilities <- rep(1/numStates, numStates)
transitionMatrix <- matrix(rep(1/numStates, numStates*numStates), nrow = numStates)

#using k-means clustering, estimate initial mu and variance values for the emission matrix
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
# When true, if a sequence is identified as a collective anomaly, it will be run through point anomaly detection
# The sequence will only be flagged as an anomaly if it has a certain percentage of point anomalies within it (ex: > 50%)
usePointDetectionOnCollectiveAnomalies <- FALSE
# When true, this will use the range of normal log-likelihoods specified in normalLogLikelihoodRanges
# Rather than use the median of training log-likelihoods +- rangeWidth
useManualNormalRange <- TRUE
normalLogLikelihoodRanges <- list(
  c(-28.12, -4.78),
  c(-112.1, -10.97),
  c(-175, -10),
  c(-200, -9)
)
windowlengths=c(5,10,15,20)
colldetection_outputfilenames=c("/home/heather/Code/cmpt-318/anom_window_5.txt",
                                "/home/heather/Code/cmpt-318/anom_window_10.txt",
                                "/home/heather/Code/cmpt-318/anom_window_15.txt",
                                "/home/heather/Code/cmpt-318/anom_window_20.txt")
for (collrun in 1:length(windowlengths)){
  sequenceSize <- windowlengths[collrun]
  trainSequences <- splitVector(trainFormattedData$x, sequenceSize)
  sink(colldetection_outputfilenames[collrun])
  
  # the last sequence in trainSequences may be less than sequenceSize (because the data points don't divide perfectly by sequenceSize)
  # This ensures such a sequence won't be tested
  # There are probably better ways to deal with this issue, but this will do for now.
  validTrainSequencesCount = if (length(trainSequences[[length(trainSequences)]]) == sequenceSize) length(trainSequences) else (length(trainSequences) - 1)
  if (useManualNormalRange) {
    minLogLikelihood <- normalLogLikelihoodRanges[[collrun]][1]
    maxLogLikelihood <- normalLogLikelihoodRanges[[collrun]][2]
  } else {
    trainLogLikelihoods <- numeric(validTrainSequencesCount)
    # Calculate log-likelihoods of train sequences
    for (i in 1:validTrainSequencesCount) {
      sequence <- trainSequences[[i]]
      yhat <- predict (hmm, sequence)
      trainLogLikelihoods[i] <- yhat$loglik
    }
    # We define a range of normal log-likelihood by centering the range on the median of the train data's log-likelihoods
    # then extend this range in an arbitrary direction (rangeWidth) depending on precision and recall results on the validation data set 
    rangeWidth = 16
    trainLogLikelihoodsMedian <- median(trainLogLikelihoods)
    minLogLikelihood <- trainLogLikelihoodsMedian - rangeWidth
    maxLogLikelihood <- trainLogLikelihoodsMedian + rangeWidth
  }
  
  # Divide test data into sequences that are each a sequenceSize in length
  # IMPORTANT: train and test sequences must be the same length (unless their log-likelihoods are scaled)
  testSequences <- splitVector(testFormattedData$x, sequenceSize)
  
  if (testOnValidation) {
    corruptSequences <- splitVector(corrupt, sequenceSize)
  }
  
  # For collective anomaly detection, we retain the concept of "threshold". Here, it basically means that on average, we should detect 
  # collective anomalies that are made up of a certain percentage (validCollectiveAnomalyRatio) of point anomalies that satisfy the same
  # threshold in the context of point anomaly detection
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
        if (usePointDetectionOnCollectiveAnomalies) {
          returnValue <- detectPointAnomalies(hmm, sequence, collectiveAnomalyThreshold, FALSE, NULL, FALSE)
          pointAnomaliesInSequence <- returnValue[[1]]
          # This differs from isAnomaly, because it checks if a sequence SHOULD have been flagged by the algorithm as an anomaly or not
          hasSufficientPointAnomalies <- (pointAnomaliesInSequence / sequenceSize) >= validCollectiveAnomalyRatio
          if (hasSufficientPointAnomalies) {
            anomalyCollectiveCount = anomalyCollectiveCount + 1
            cat("1,1\n")
          }
        } else{
          anomalyCollectiveCount = anomalyCollectiveCount + 1
          cat("1,1\n")
        }
      }
      else{
        cat ("0,0\n")
      }
      if (testOnValidation) {
        returnValue <- detectPointAnomalies(hmm, sequence, collectiveAnomalyThreshold, testOnValidation, corruptSequences[[i]], FALSE)
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
  sink()
  
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
}





#-----------Point Anomaly Detection------------------------------
# The general idea here is to look at each point individually, and compare it to an expected value.
# Here, the expected value is determined by the most probable state sequence. 
# You look at the most probable state assigned to a data point (in predict()), and find that state's output 
# emission mean (specified in the HMM). This mean represents an expected normal value for the given state
pointdetection_outputfilenames=c("/home/heather/Code/cmpt-318/point_threshold_0point5.txt",
                                 "/home/heather/Code/cmpt-318/point_threshold_0point75.txt",
                                 "/home/heather/Code/cmpt-318/point_threshold_1.txt",
                                 "/home/heather/Code/cmpt-318/point_threshold_1point25.txt",
                                 "/home/heather/Code/cmpt-318/point_threshold_1point5.txt",
                                 "/home/heather/Code/cmpt-318/point_threshold_2.txt")
# pointdetection_outputfilenames=c("C:\\Users\\Trevor\\Google Drive\\School\\SFU\\Year 4\\Summer 2017\\CMPT 318\\Project\\CMPT318\\results\\point_threshold_0point5.txt",
#                                  "C:\\Users\\Trevor\\Google Drive\\School\\SFU\\Year 4\\Summer 2017\\CMPT 318\\Project\\CMPT318\\results\\point_threshold_0point75.txt",
#                                  "C:\\Users\\Trevor\\Google Drive\\School\\SFU\\Year 4\\Summer 2017\\CMPT 318\\Project\\CMPT318\\results\\point_threshold_1.txt",
#                                  "C:\\Users\\Trevor\\Google Drive\\School\\SFU\\Year 4\\Summer 2017\\CMPT 318\\Project\\CMPT318\\results\\point_threshold_1point25.txt",
#                                  "C:\\Users\\Trevor\\Google Drive\\School\\SFU\\Year 4\\Summer 2017\\CMPT 318\\Project\\CMPT318\\results\\point_threshold_1point5.txt",
#                                  "C:\\Users\\Trevor\\Google Drive\\School\\SFU\\Year 4\\Summer 2017\\CMPT 318\\Project\\CMPT318\\results\\point_threshold_2.txt")
thresholds=c(0.5,0.75,1,1.25,1.5,2)

for (threshrun in 1:length(thresholds)){
  threshold=thresholds[threshrun]
  sink(pointdetection_outputfilenames[threshrun])
  returnValue <- detectPointAnomalies(hmm, testFormattedData$x, threshold, testOnValidation, corrupt, TRUE)
  anomalyPointCount <- returnValue[[1]]
  truePositiveCount <- returnValue[[2]]
  sink()
  
  # Report Results
  anomalyPointPercentage = 100 * (anomalyPointCount/length(testFormattedData$x))
  cat("Point Anomaly Percentage: ", anomalyPointPercentage, "%\n")
  cat("Point Anomaly Count: ", anomalyPointCount, "\n")
  if (testOnValidation) {
    precision <- 100 * (truePositiveCount / anomalyPointCount)
    actualAnomalyCount <- sum(corrupt)
    recall <- 100 * (truePositiveCount / actualAnomalyCount)
    cat("Precision: ", precision, "%\n")
    cat("Recall: ", recall, "%\n\n")
  }
}
