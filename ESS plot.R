library(ggplot2)
library(tidyverse)
ESS <- data.frame(x=1:162,y=maic.aess.mean)
ggplot(ESS, aes(x = x, y = y)) +
  geom_segment(aes(x = x, xend = x, y = 50, yend = y), color = "grey") +  # Grey lines from 0 to points
  geom_point(size = 1, color = "red",  fill = alpha("grey", 0.3),alpha = 0.7, shape = 21, stroke = 1) +
  geom_hline(yintercept = 50, linetype = "solid", color = "black", size = 0.7) +
  annotate("text", x = 160, y = 60, label = "y = 50", color = "black") +
  xlab("Scenarios in Parameter Combination Order") +  # X-axis label
  ylab("Effective Sample Size") +  # Y-axis label
  scale_y_continuous(limits = c(0, 300), breaks = seq(0, 300, by = 50)) +  # Adjust the y-axis limit and breaks
  scale_x_continuous(breaks = seq(0, 162, by = 9), limits = c(1, 162)) +
  theme_minimal()


 
