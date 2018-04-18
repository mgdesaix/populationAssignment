################################################################################################################################
####----------------------------------------------Leave One Out-------------------------------------------------####
################################################################################################################################
## March 16, 2018
# Matt DeSaix

library(dplyr)
library(magrittr)
library(gstudio)

data <- readRDS("df_AB_gstudio_2alleles.RDS")
start <- 1
fst <- readRDS("./hierfstat/loc_fst_ranked_2alleles.RDS")

fst.1000 <- fst[1:1000,]
data.sub <- data[,colnames(data) %in% fst.1000$Locus]
data.loci.all <- cbind(data[,1:3], data.sub)


leave.one.out <- function(data.loci.all, start, fst){
  
  ## A function to calculate genotype likelihood at each locus
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
  
  ## Got this log.sum function from stackexchange
  log.sum <- function(y, precision = log(10^17)){
    if (length(y) == 0) return(1)
    log.plus <- function(a, b) ifelse(abs(b-a) > precision,
                                      max(b,a), max(b,a) + log(1 + exp(-abs(b-a))))
    y <- sort(y)
    x <- y[1]
    for ( z in y[-1]) x <- log.plus(x,z)
    return(x)
  }
  
  ass.mat <- as.data.frame(matrix(nrow = nrow(data.loci.all), ncol  = 5))
  rownames(ass.mat) <- row.names(data.loci.all$ID)
  
 
  for(i in start:nrow(data.loci.all )){
    #system.time({ 
    message( Sys.time())
    message("Processing individual ", i, " of ", nrow(data.loci.all))
    
    df_AB.hold_out <- data.loci.all[i,]
    df_AB.training <- data.loci.all[-i,]
    num_pops <- length(unique(df_AB.training$Pop))
    
    ## list all loci that have both alleles present across all populations
    loci.all.strata.training <- df_AB.training %>%
      frequencies( stratum = "Pop") %>%
      count(Locus) %>%
      # this filters to the frequency count, if both alleles are present then there are two counts for that SNP
      # therefore I want to select only the loci that have a count of twice the number of populations, indicating both alleles are present in each pop
      filter(n == ( num_pops * 2) ) %>%
      select(Locus)
    
    fst.sub <- fst[fst$Locus %in% loci.all.strata.training$Locus,]
    fst.sub <- fst.sub[ order( fst.sub$FST, decreasing = TRUE), ]
    sub.sample <- fst.sub[1:600,1]
    
    # Create Tibble of data with loci of both alleles for every pop
    df_AB_2alleles.training <- df_AB.training %>%
      select(sub.sample) %>%
      cbind(df_AB.training[,1:3], .)
    
    # Allele frequencies
    training.freq <- df_AB_2alleles.training %>%
      frequencies(stratum = "Pop")
    
    # Create hold out data based on frequencies from training
    holdout.gen <- df_AB.hold_out %>%
      select(sub.sample)
    rownames(holdout.gen) <- df_AB.hold_out[3]
    
    
    # determine genotype likelihood
    prob.list <- sapply(1:ncol(holdout.gen), function(x) {pop.prob2(training.freq, holdout.gen[x])})
    rownames(prob.list) <- unique(training.freq$Stratum)
    colnames(prob.list) <- colnames(holdout.gen)
    
    # Atlantic pops
    atl <- c("PopA", "PopE", "PopF", "PopG", "PopH")
    # Mississippi pops
    mrv <- c("PopJ", "PopK", "PopL", "PopM", "PopN", "PopO")
    
    ass.mat[i,5] <- ncol(prob.list)
    
    freq.sum <- apply(prob.list, 1, function(x) log.sum(x) )
    
    # data frame with rows=individuals, columns = decreasing order of population assignment
    pop.mat <- as.data.frame( matrix(nrow = 1, ncol = length(unique(training.freq$Stratum))))
    rownames(pop.mat) <- df_AB.hold_out[3]
    colnames(pop.mat) <- c("MaxLogSum", 2:length(unique(training.freq$Stratum)))
    
    pop.mat[1,] <- names( freq.sum[order(freq.sum, decreasing = TRUE)] )
    
    correct.pop <- substr(rownames(pop.mat), 1, 4) == pop.mat[,1]
    correct.reg <- (substr(rownames(pop.mat), 1, 4) %in% atl) & (pop.mat[,1] %in% atl) | (substr(rownames(pop.mat), 1, 4) %in% mrv) & (pop.mat[,1] %in% mrv)
    
    ass.pop <- names(freq.sum[freq.sum == max(freq.sum)])
    ass.mat[i, 1] <- ass.pop
    ass.mat[i, 2] <- correct.pop
    if (ass.pop %in% atl){
      ass.mat[i,3] <- "Atlantic"
    } else{
      ass.mat[i,3] <- "Mississippi"
    }
    ass.mat[i, 4] <- correct.reg
    
    saveRDS(ass.mat[i,], file = eval(paste("./training_data/leave_one_out/LOO_", data.loci.all[i,3], ".RDS", sep = "")))

   # })
  }
  saveRDS(ass.mat, file = "./training_data/leave_one_out/LOO_all.RDS")
  rownames(ass.mat) <- data.loci.all$ID
  return(ass.mat)
}


hope.this.works <- leave.one.out(data.loci.all = data.loci.all, start = 1, fst = fst.1000)



















