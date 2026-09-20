
#Dendrogram Grouped


# Library
library(tidyverse)
library(dendextend)


# Plot

plot(dhc)

# Color in function of the cluster
par(mar=c(1,1,1,7))
dhc %>%
  set("labels_col", value = c("red", "blue" , "Darkgreen", "Purple"), k=4) %>%
  set("branches_k_color", value = c("red", "blue" , "Darkgreen", "Purple"), k = 4) %>%
  plot(horiz=FALSE, axes=FALSE)