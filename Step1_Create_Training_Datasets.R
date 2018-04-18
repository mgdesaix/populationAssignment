################################################################################################################################
####----------------------------------------------Create Training Data Sets-------------------------------------------------####
################################################################################################################################
## February 23, 2018
# Matt DeSaix

library(gstudio)
library(magrittr)
library(dplyr)
# df_AB_2alleles.rda loads data.loci.all (which I should change the name to match the .rda file...)
# load("df_AB_2alleles.rda")

AB_gstudio <- readRDS("final_snps_AB_LDgstudio.RDS")

create.training.data <- function(data.loci.all){

  # take a random sample of 1/3rd of the data
  ind <- sapply( unique( data.loci.all$Pop), function(x) sample( which(data.loci.all$Pop == x) , length(which(data.loci.all$Pop == x))/3)) %>%
    unlist()
  # 2/3rd as trainings
  df_AB.training <- data.loci.all[-ind,]
  # 1/3rd hold out
  df_AB.hold_out <- data.loci.all[ind,]
  num_pops <- length(unique(df_AB.training$Pop))
  
  ## list all loci that have both alleles present across all populations
  loci.all.strata.training <- df_AB.training %>%
    frequencies( stratum = "Pop") %>%
    count(Locus) %>%
    # this filters to the frequency count, if both alleles are present then there are two counts for that SNP
    # therefore I want to select only the loci that have a count of twice the number of populations, indicating both alleles are present in each pop
    filter(n == ( num_pops * 2) ) %>%
    select(Locus) 
  
  # Create Tibble of data with loci of both alleles for every pop
  df_AB_2alleles.training <- df_AB.training %>%
    select(loci.all.strata.training$Locus) %>%
    cbind(df_AB.training[,1:3], .)
  
  # Allele frequencies
  training.freq <- df_AB_2alleles.training %>%
    frequencies(stratum = "Pop")
  
  # Create hold out data based on frequencies from training
  df_AB_2alleles.hold_out <- df_AB.hold_out %>%
    select(loci.all.strata.training$Locus) %>%
    cbind(df_AB.hold_out[,1:3], .)
  
  output <- list("Training" = training.freq, "Holdout" = df_AB_2alleles.hold_out)
  return(output)
}


# data.sub <- AB_gstudio[c(1:10, 21:30, 81:90), 1:100]

for(i in 6:20){
  message("Processing individual ", i, " of ", 20)
  training <- create.training.data(data.loci.all = AB_gstudio)
  name <- paste("training.set", i, "RDS", sep = ".")
  out.file <- paste("./training_data/", name, sep = "")
  saveRDS(training, file = out.file )
}






