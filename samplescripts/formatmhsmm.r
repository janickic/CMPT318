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