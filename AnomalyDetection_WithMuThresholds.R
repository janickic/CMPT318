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
trainDataset <- read_csv("/home/heather/Code/cmpt-318/train.txt")
trainGlobalActivePower <- trainDataset$Global_active_power
trainGlobalActivePower <- trainGlobalActivePower[!is.na(trainGlobalActivePower)]
trainFormattedData <- formatMhsmm(data.frame(trainGlobalActivePower))

# do the same thing for Test data
testDataset <- read_csv("/home/heather/Code/cmpt-318/test1.txt")
testGlobalActivePower <- testDataset$Global_active_power
testGlobalActivePower <- testGlobalActivePower[!is.na(testGlobalActivePower)]
testFormattedData <- formatMhsmm(data.frame(testGlobalActivePower))

# specify initial HMM parameter values
# values were determined from analysis of data set
# These are currently set to whatever values Evan was using in his analysis
numStates <- 5

#using k-means clustering, estimate initial mu and sigma values
kgap <- kmeans(trainGlobalActivePower,5,iter.max=10)
muinit <- kgap$centers
varinit = numeric(numStates)
for (i in 1:numStates){
  varinit[i] <- kgap$withinss[i]/(kgap$size[i]-1)
}

initialStateProbabilities <- rep(1/numStates, numStates)
transitionMatrix <- matrix(rep(1/numStates, numStates*numStates), nrow = numStates)
emissionMatrix <- list(mu = muinit, sigma = varinit) 
#emissionMatrix <- list(mu = c(0.85,3.426), sigma = c(2,1)) 

# create a starting model with initial parameters before performing Expectation Maximization algorithm (training)
startmodel <- hmmspec(init = initialStateProbabilities, 
                      trans = transitionMatrix, 
                      parms.emis = emissionMatrix, 
                      dens.emis = dnorm.hsmm)

# perform Expectation Maximization to obtain a trained HMM model
hmm <- hmmfit(trainFormattedData$x, startmodel, mstep = mstep.norm, maxit = 200)

# print resulting HMM parameters
summary(hmm)

yhatpoint <- predict(hmm,trainFormattedData$x,method="smoothed")

threshold <- 1
inanom <-0

percent<-0
for (t in 2:100000){
  obs<-trainFormattedData$x[t]
  
  probs <- 0
  exp <- 0
  for (currState in 1:5){
    marginalize <- 0 #joint prob for past two states given observations up to now
    for (prevState in 1:5){
      posterior <- yhatpoint$p[t-1,prevState] #the posterior distribution of the previous state
      transprob <- hmm$model$transition[prevState,currState] #possibility for state j given previous state
      marginalize <- marginalize + posterior*transprob
    }
    
    probs <- probs + marginalize*dnorm.hsmm(obs,currState,hmm$model) #likelihood for obs given state
    probs <- min(1,probs)
    exp <- exp + probs*muinit[currState]
  }
  #print(abs(obs-exp))
  if (abs(exp-obs)>threshold){
    if (inanom==0){
      print(t)
      print ("to")
      inanom <-1
    }
  }
  else if (inanom==1) {
    print (t)
    inanom <-0
    print ("-----------")
  }

  if (abs(exp-obs)>threshold){
    percent <- percent+1
  }
}
print(percent)