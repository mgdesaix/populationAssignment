################################################################################################################################
####--------------------------------------Assignment with Log Sums, multiple SNP sets---------------------------------------####
################################################################################################################################
# Matt DeSaix
# March 16, 2018




## Functions Needed
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
multi.assign <- function(prob.freq.per.locus, snps, fst){
  # Atlantic pops
  atl <- c("PopA", "PopE", "PopF", "PopG", "PopH")
  # Mississippi pops
  mrv <- c("PopJ", "PopK", "PopL", "PopM", "PopN", "PopO")
  # Assignment list
  assignment.list <- list()
  list.names <- c()
  
  # data frame with individuals, population assigned, T/F, region assigned, T/F
  ass.mat <- as.data.frame( matrix(nrow = length(names(prob.freq.per.locus)), ncol = (4 * length(snps) )))
  rownames(ass.mat) <- names(prob.freq.per.locus)
  
  sum.mat <- as.data.frame(matrix(nrow = length(snps), ncol = 3))
  colnames(sum.mat) <- c("NumberOfSNPs", "PopAssignment", "RegionAssignment")
  
  # Produces the population assignment for each snp set
  for(j in 1:length(snps)){
    # data frame with rows=individuals, columns = decreasing order of population assignment
    pop.mat <- as.data.frame( matrix(nrow = length(names(prob.freq.per.locus)), ncol = 11) )
    rownames(pop.mat) <- names(prob.freq.per.locus)
    colnames(pop.mat) <- c("MaxLogSum", 2:11)
    
    # data frame with rows = individuals, columns = log sum value
    prob.mat <- as.data.frame(matrix(nrow = length(names(prob.freq.per.locus)), ncol = 11))
    rownames(prob.mat) <- names(prob.freq.per.locus)
    colnames(prob.mat) <- c("MaxLogSum", 2:11)
    
    list.names[2*j -1] <- paste("Assigned.Pop", snps[j], sep = "")
    list.names[2*j] <- paste("Log.Sum.Values", snps[j], sep = "")
    
    fst.sub <- fst[fst$Locus %in% colnames( prob.freq.per.locus[[1]]),]
    fst.sub <- fst.sub[ order(fst.sub$FST, decreasing = TRUE), ]
    
    sub.sample <- fst.sub[1:snps[j], 1]
    
    # takes the log.sum of each element in the list of prob.freq.per.locus
    for(i in 1:length(names(prob.freq.per.locus))){
      
      freq.sum <- apply(prob.freq.per.locus[[i]][,sub.sample], 1, function(x) log.sum(x) )
      pop.mat[i,] <- names( freq.sum[order(freq.sum, decreasing = TRUE)] )
      
      prob.mat[i,] <- freq.sum[order(freq.sum, decreasing = TRUE)] 
      
      ass.pop <- names(freq.sum[freq.sum == max(freq.sum)])
      ass.mat[i, 4*j - 3] <- ass.pop
      if (ass.pop %in% atl){
        ass.mat[i,4 * j-1] <- "Atlantic"
      } else{
        ass.mat[i,4 * j-1] <- "Mississippi"
      }
    }
    assignment.list[[2*j - 1]] <- pop.mat
    assignment.list[[2*j]] <- prob.mat
    
    correct.pop <- substr(rownames(pop.mat), 1, 4) == pop.mat[,1]
    correct.reg <- (substr(rownames(pop.mat), 1, 4) %in% atl) & (pop.mat[,1] %in% atl) | (substr(rownames(pop.mat), 1, 4) %in% mrv) & (pop.mat[,1] %in% mrv)
    ass.mat[, (4*j-2)] <- correct.pop
    ass.mat[, (4*j) ] <- correct.reg
    c.pop.val <- length(correct.pop[correct.pop == TRUE]) / length(correct.pop)
    c.reg.val <- length(correct.reg[correct.reg == TRUE]) / length(correct.reg)
    percent.correct.set <- c(snps[j], c.pop.val, c.reg.val)
    
    sum.mat[j,] <- percent.correct.set
  }
  
  assignment.list[[2*j + 1]] <- sum.mat
  assignment.list[[2*j + 2]] <- ass.mat
  
  list.names[2*j + 1] <- "Assignment.Summary"
  list.names[2*j + 2] <- "Assignment.Success"
  
  names(assignment.list) <- list.names
  
  return(assignment.list)
}



# Variables needed

set <- c(200, 400, 600, 800, 1000, 1200, 1500, 2000, 3000)
loc.fst <- readRDS("./hierfstat/loc_fst_ranked.RDS")



files <- list.files("./training_data/likelihood", pattern = ".RDS")


for(file in files){
  data <- readRDS( paste( "./training_data/likelihood/", file, sep = "" ) )
  
  assignment.success <- multi.assign( prob.freq.per.locus = data, snps = set, fst = loc.fst)
  
  file <- gsub(".RDS", "", file)
  name <- paste(file, ".assignment.RDS", sep = "")
  out.file <- paste("./training_data/assignment/", name, sep = "")
  
  saveRDS(assignment.success, file = out.file)
  rm(assignment.success)
  
}









system.time({
  do.it <- multi.assign( prob.freq.per.locus = data, snps = set, fst = loc.fst)
}) 

do.it$Assignment.Summary


