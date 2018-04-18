library(dplyr)
library(gridExtra)
library(ggplot2)
library(magrittr)

files <- list.files(pattern = ".RDS", path = "./training_data/assignment/")
l <- list()
i <- 0
for(file in files){
  i <- i + 1
  ## the last number is dependent on which is "Assignment.Summary" list
  l[[i]] <- readRDS( paste("./training_data/assignment/", file, sep = "") )[[19]]
}
assignment.df <- Reduce(function(...) merge(...,all = T), l)

## I used df.sub when I made plots for a presentation, otherwise use them all!
df.sub <- assignment.df %>%
  filter( NumberOfSNPs == 200 | NumberOfSNPs == 400 | NumberOfSNPs == 600 | 
            NumberOfSNPs == 1000 | NumberOfSNPs == 2000 | NumberOfSNPs == 3000)
snp.levels <- unique( df.sub$NumberOfSNPs)

df.sub$NumberOfSNPs <- factor( df.sub$NumberOfSNPs, levels = snp.levels)

snp.levels <- unique( assignment.df$NumberOfSNPs)
assignment.df$NumberOfSNPs <- factor(assignment.df$NumberOfSNPs, levels = snp.levels)

assignment.df %>%
  filter( NumberOfSNPs == 600) %>%
  select( RegionAssignment) %>%
  unlist() %>%
  sd()


p1 <- ggplot(assignment.df, aes(x = NumberOfSNPs, 
                                y = PopAssignment, 
                                group = NumberOfSNPs)) + 
  theme_bw() +
  geom_boxplot() +
  ggtitle("Site") +
  xlab("Number Of SNPs") +
  scale_y_continuous( name = "Assignment Success", limits = c(0.4,1)) +
  theme(axis.text = element_text( size = 12),
        axis.title = element_text( size = 16),
        plot.title = element_text( size = 16)
        )

p2 <- ggplot(assignment.df, aes(x = NumberOfSNPs, 
                                y = RegionAssignment, 
                                group = NumberOfSNPs)) + 
  geom_boxplot() + 
  theme_bw() +
  ggtitle("Region") +
  xlab("Number Of SNPs") +
  scale_y_continuous( name = "Assignment Success", limits = c(0.4,1)) +
  theme(axis.text = element_text( size = 12),
        axis.title = element_text( size = 16),
        plot.title = element_text( size = 16)
  )

p3 <- grid.arrange(p1, p2)
plot(p3)

ggsave("./quality_plots/thesis/assignment_success.png", p3, dpi = 300, width = 8, height = 7)

filter(assignment.df, NumberOfSNPs == 600) %>%
  summary()
filter(assignment.df, NumberOfSNPs == 800) %>%
  summary()


##### Transparent


p1 <- ggplot(df.sub, aes(x = NumberOfSNPs, 
                                y = PopAssignment, 
                                group = NumberOfSNPs)) + 
  geom_boxplot( fill = "indianred", color = "indianred3", alpha = 0.60) +
  ggtitle("Sampling Site") +
  ylab("Assignment Success") + 
  xlab("Number Of SNPs") + 
  theme( axis.text = element_text( size = 20, color = "gray85"),
         axis.title = element_text( size = 22, color = "gray85"),
         plot.title = element_text( size = 26, color = "gray85"),
         legend.position = "none",
         panel.background = element_rect(fill = "transparent"),
         plot.background = element_rect(fill = "transparent"),
         panel.grid.major = element_line(color = "gray30", size = 0.5, linetype = "solid"))
p1
ggsave("./quality_plots/assignment_success_tr_site_aos.png", p1, bg = "transparent", dpi = 300, width = 8, height = 4)

p2 <- ggplot(df.sub, aes(x = NumberOfSNPs, 
                                y = RegionAssignment, 
                                group = NumberOfSNPs)) + 
  geom_boxplot( fill = "indianred", color = "indianred3", alpha = 0.60) +
  ggtitle("Region") +
  ylab("Assignment Success") + 
  xlab("Number Of SNPs") + 
  theme( axis.text = element_text( size = 20, color = "gray85"),
         axis.title = element_text( size = 22, color = "gray85"),
         plot.title = element_text( size = 26, color = "gray85"),
         legend.position = "none",
         panel.background = element_rect(fill = "transparent"),
         plot.background = element_rect(fill = "transparent"),
         panel.grid.major = element_line(color = "gray35", size = 0.5, linetype = "solid"))
ggsave("./quality_plots/assignment_success_tr_region_aos.png", p2, bg = "transparent", dpi = 300, width = 8, height = 4)


p3 <- grid.arrange(p1, p2)
plot(p3)
