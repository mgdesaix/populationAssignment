################################################################################################################################
####------------------------------------------Leave One Out: Results Display---------------------------------------------####
################################################################################################################################
## March 16, 2018
# Matt DeSaix
library(ggplot2)
library(dplyr)
library(magrittr)
library(ggrepel)


loo.df <- readRDS("./training_data/leave_one_out/LOO_all.RDS")
files <- list.files(pattern = "LOO_Pop", path = "./training_data/leave_one_out")


colnames(loo.df) <- c( "AssignedPopulation", "PopT_F", "AssignedRegion", "RegT_F", "NumberOfTotalSNPs")
rownames(loo.df) <- substr(files, 5, 11)

head(loo.df)


correct_pop <- nrow(loo.df[loo.df$PopT_F == TRUE,])/nrow(loo.df)
correct_reg <- nrow(loo.df[loo.df$RegT_F == TRUE,]) / nrow(loo.df)
correct_pop
correct_reg

pops <- unique(substr(files, 5, 8))
assignment.mat <- matrix(data = 0, nrow = length(pops), ncol = length(pops))
colnames(assignment.mat) <- pops
rownames(assignment.mat) <- pops


actual <- substr(rownames(loo.df), 1, 4)
assigned <- loo.df[,1]

for( i in 1:length(actual)){
  assignment.mat[actual[i], assigned[i]] <- assignment.mat[actual[i], assigned[i]] + 1
}
assignment.mat


assignment.plot <- data.frame("Sampled_pop" = actual, "Inferred_pop" = assigned) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq != 0)

assignment.plot$Inferred_pop <- as.character(assignment.plot$Inferred_pop)
assignment.plot$Sampled_pop <- as.character(assignment.plot$Sampled_pop)

zeros <- rep(0, nrow(assignment.plot))
assignment.plot$Correct_assign <- zeros

for(i in 1:nrow(assignment.plot)){
  if(assignment.plot[i,1] == assignment.plot[i,2]){
    assignment.plot[i,4] <- 1
  }
}

key <- c("VA1", "VA2", "VA3", "NC1", "SC1", "OH1", "LA1", "LA2", "LA3", "AR1", "WI1")
key.df <- data.frame("OldNames" = colnames(assignment.mat), "NewNames" = key)

assignment.plot$Sampled_pop <- key.df$NewNames[ match( assignment.plot$Sampled_pop, key.df$OldNames)]
assignment.plot$Inferred_pop <- key.df$NewNames[ match( assignment.plot$Inferred_pop, key.df$OldNames)]

assignment.plot$Sampled_pop2 <- factor( assignment.plot$Sampled_pop, levels = key)
assignment.plot$Inferred_pop2 <- factor( assignment.plot$Inferred_pop, levels = key)


# Plot with Raw values
ggplot(assignment.plot, aes(x = Inferred_pop2, y = Sampled_pop2)) + 
  geom_point(shape = 1, size = 10, aes(color = factor(Correct_assign))) + 
  geom_text(aes(label = Freq)) +
  theme_bw() +
  ggtitle("Assignment Success") +
  xlab("Inferred Population") +
  ylab("Sampled Population") +
  guides(color = FALSE) + 
  geom_hline( yintercept = 5.5, linetype = "dashed", alpha = 0.5) + 
  geom_vline( xintercept = 5.5, linetype = "dashed", alpha = 0.5)


# Plot with Percent values
assignment.plot$n <- zeros
assignment.plot$Percent_assign <- zeros
for(i in 1:nrow(assignment.plot)){
  pop <- assignment.plot$Sampled_pop[i]
  total <- sum(assignment.plot[assignment.plot$Sampled_pop == pop, 3])
  assignment.plot$Percent_assign[i] <- round(assignment.plot[i,3]/total * 100)
  assignment.plot$n[i] <- total
}


ggplot(assignment.plot, aes(x = Inferred_pop2, y = Sampled_pop2)) + 
  geom_point(shape = 1, size = 15, aes(color = factor(Correct_assign))) + 
  geom_text(size = 5, aes(label = Percent_assign)) +
  theme_bw() +
  ggtitle("Leave-one-out Assignment Success") +
  xlab("Inferred Population") +
  ylab("Sampled Population") +
  guides(color = FALSE) + 
  annotate( "text", x = 2.8, y = 5.2, label = "Atlantic Coast", size = 5) + 
  annotate( "text", x = 8.5, y = 11.2, label = "Mississippi River Valley", size = 5) +
  geom_rect( aes( xmin = 0, xmax = 5.5, ymin = 0, ymax = 5.5), linetype = 2, fill = "gray85", alpha = 0.01) + 
  geom_rect( aes( xmin = 5.5, xmax = 11.6, ymin = 5.5, ymax = 11.6), linetype = 2, fill = "grey85", alpha = 0.01) +
  geom_hline( yintercept = 5.5, linetype = "dashed", alpha = 0.5) + 
  geom_vline( xintercept = 5.5, linetype = "dashed", alpha = 0.5) + 
  theme( axis.text = element_text( size = 12 ),
         axis.title = element_text( size = 16),
         legend.position = "none")
  

ggsave("./quality_plots/thesis/LOO_success2.png", dpi = 300, width = 9, height = 8)


key <- c("VA1", "VA2", "VA3", "NC1", "SC1", "OH1", "LA1", "LA2", "LA3", "AR1", "WI1")
head(assignment.plot)
correct.df <- assignment.plot %>%
  filter(Correct_assign == 1)

correct.df$size_group <- 0

for(i in 1:nrow( correct.df )){
  if( correct.df$n[i] <10){
    correct.df$size_group[i] <- 1
  } else if( correct.df$n[i] >= 10 & correct.df$n[i] < 20){
    correct.df$size_group[i] <- 2
  } else if( correct.df$n[i] >= 20){
    correct.df$size_group[i] <- 3
  }
}
correct.df 

ggplot(correct.df, aes( x = n, y = Percent_assign/100)) + 
  geom_point() +
  geom_label_repel(aes(label = Sampled_pop), size = 4 ) +
  xlab("Sample Size") +
  scale_y_continuous( name = "Assignment Success", limits = c(0,1)) + 
  theme_bw() +
  theme(axis.text = element_text( size = 12),
        axis.title = element_text( size = 16)
        )

ggsave("./quality_plots/thesis/success_by_missing.png", dpi = 300, width = 8, height = 6)

lm1 <- lm(Percent_assign ~ n, data = correct.df)
summary(lm1)
lm2 <- lm(Percent_assign ~ poly(n,2), data = correct.df)
summary(lm2)

## Plot assignment success by sample size of site, with geom_smooth

ggplot(correct.df, aes( x = n, y = Percent_assign/100)) + 
  geom_point() +
  geom_smooth(method = 'lm', formula = y ~ poly(x,2) ) +
  xlab("Sample Size") +
  scale_y_continuous( name = "Assignment Success", limits = c(0,1)) + 
  theme_bw() +
  theme(axis.text = element_text( size = 12),
        axis.title = element_text( size = 16)
  )


summary(correct.df$Percent_assign)



#### transparent background
assignment.plot$color <- rep(1, nrow(assignment.plot))

ggplot(assignment.plot, aes(x = Inferred_pop2, y = Sampled_pop2)) + 
  geom_point(size = 15, aes(color = factor(Correct_assign)), alpha = 0.5) + 
  geom_text(size = 5, aes(label = Percent_assign), color = "white") +
  ggtitle("Leave-one-out Assignment Success") +
  xlab("Inferred Population") +
  ylab("Sampled Population") +
  guides(color = FALSE) + 
  # annotate( "text", x = 2.8, y = 5.2, label = "Atlantic Coast", size = 5, color = "gray85") + 
  # annotate( "text", x = 8.5, y = 11.2, label = "Mississippi River Valley", size = 5, color = "gray85") +
  geom_rect( aes( xmin = 0, xmax = 5.5, ymin = 0, ymax = 5.5), linetype = 2, fill = "gray85", alpha = 0.01) + 
  geom_rect( aes( xmin = 5.5, xmax = 11.6, ymin = 5.5, ymax = 11.6), linetype = 2, fill = "grey85", alpha = 0.01) +
  geom_hline( yintercept = 5.5, linetype = "dashed", alpha = 0.5, color = "grey85") + 
  geom_vline( xintercept = 5.5, linetype = "dashed", alpha = 0.5, color = "gray85") + 
  theme( axis.text = element_text( size = 14, color = "gray85"),
         axis.title = element_text( size = 20, color = "gray85"),
         plot.title = element_text( size = 20, color = "gray85"),
         legend.position = "none",
         panel.background = element_rect(fill = "transparent"),
         plot.background = element_rect(fill = "transparent"),
         panel.grid.major = element_line(color = "gray60", size = 0.5, linetype = "solid"))


ggsave("./quality_plots/LOO_success3_2.png", bg = "transparent", dpi = 300, width = 8, height = 7)



#### Plot assignment success by site size

n <- c(8, 9, 12, 13, 13, 18, 19, 20, 20, 21, 22)
p <- c(.50, .44, .67, .92, 1, .94, .95, .95, .95, 1, 1 )

np.df <- data.frame( "Sample_size" = n, "AssignSuccess" = p)

pq <- ggplot(np.df, aes(x = Sample_size, y = AssignSuccess)) +
  geom_point()
pq



