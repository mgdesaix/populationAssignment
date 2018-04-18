################################################################################################################################
####----------------------------------------------Assignment Likelihood-------------------------------------------------####
################################################################################################################################
## February 23, 2018
# Matt DeSaix



######################## Part 2)  Assignment likelihood ########################  
# import all training data

# using A:B gstudio format


# Two functions required
pop.prob2 <- function(training, holdout){
  p <-  training %>%
    filter( Locus == colnames(holdout), Allele == "A") 
  p <- as.numeric(as.character(p$Frequency))
  if( holdout == "A:B") {
    prob <- (2 * p * (1-p))
    return(prob)
  } else if ( holdout == "A:A") {
    prob <- (p^2)
    return(prob)
  } else if ( holdout == "B:B") {
    prob <- ((1-p)^2)
    return(prob)
  } else {
    prob <- rep(NA, length(unique(training$Stratum)))
    return(prob)
  }
}
multi.pop.prob2 <- function(training_set){
  training <- training_set$Training
  holdout.gen <- training_set$Holdout[,4:ncol(training_set$Holdout)]
  holdout.meta <- training_set$Holdout[1:nrow(holdout.gen), 1:3]
  prob.list <- list()
  for( i in 1: nrow(holdout.gen)){
    ind <- holdout.gen[i,]
    # pbsapply is 'sapply' with a progress bar
    prob.list[[i]] <- pbsapply(1:ncol(holdout.gen), function(x) {pop.prob2(training, ind[x])})
    colnames(prob.list[[i]]) <- colnames(ind)
    rownames(prob.list[[i]]) <- unique(training$Stratum)
  }
  names(prob.list) <- holdout.meta$ID
  return(prob.list)
}



files <- list.files(pattern = "RDS", path = "./training_data")
i <- 0
for(file in files){
  i <- i + 1
  training.data <- readRDS(paste("./training_data/", file, sep = ""))
  training.set.likelihood <- multi.pop.prob2(training_set = training.data)
  
  name <- paste("training.set.likelihood.", i,".RDS", sep = "")
  out.file <- paste("./training_data2/", name, sep = "")
  saveRDS(training.set.likelihood, file = out.file)
}



what_is_this <- readRDS("./training_data/likelihood/training.set.likelihood.10_2.RDS")


files <- list.files(pattern = "RDS", path = "./training_data/likelihood")
i <- 0
col.list <- c()
for( file in files){
  i <- i + 1
  training.data <- readRDS(paste("./training_data/likelihood/", file, sep = ""))
  col.list[i] <- ncol(training.data[[1]])
}











