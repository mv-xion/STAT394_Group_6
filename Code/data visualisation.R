# install packages
install.packages(c("reshape2", "ggplot2", "ggcorrplot", "GGally", "psych"))
install.packages("devtools")
devtools::install_github("vqv/ggbiplot")

library(reshape2)
library(ggplot2)
library(ggcorrplot)
library(GGally)
library(reshape2)
library(ggbiplot)
library(psych)

#----------------------------------------------------------------
# Data preparation 
# Reading the data
data_spa_waters <- read.csv("C:\\Users\\ixion\\NZ\\PhD\\STAT394_Group_6\\Data\\Hydrochemical Data.csv") # Modify this one
head(data_spa_waters) # Showing data


# Cleaning data
# Removing the first 2 columns (number, place name) 
data_spa_waters_classes = data_spa_waters[-1, -c(1, 2)]
ncolumns_data_spa_waters_classes = ncol(data_spa_waters_classes)

# Consider data as numeric, round to 3 decimals except class column
data_spa_waters_classes[-c(ncolumns_data_spa_waters_classes-1, ncolumns_data_spa_waters_classes)] = lapply(data_spa_waters_classes[-c(ncolumns_data_spa_waters_classes-1, ncolumns_data_spa_waters_classes)], function(x) round(as.numeric(x), 3))
# Set Structure and region as factors
data_spa_waters_classes$ Geological.structure <- as.factor(data_spa_waters_classes$ Geological.structure)
data_spa_waters_classes$ Region <- as.factor(data_spa_waters_classes$ Region)

#Show the data
head(data_spa_waters_classes)
summary(data_spa_waters_classes)

# Remove the last 2 columns to get the raw numerical data
data_spa_waters_raw <- data_spa_waters_classes[,-c(ncolumns_data_spa_waters_classes-1, ncolumns_data_spa_waters_classes)]

#Show the data
head(data_spa_waters_raw)
summary(data_spa_waters_raw)

# Correlation of chemical compounds
cor(data_spa_waters_raw)

# Pairsplot of the dara
ggpairs(data_spa_waters_raw)

#Plot correlogram
ggcorrplot(cor(data_spa_waters_raw), 
           method = "circle",
           type = "lower",
           lab=TRUE)

#Plot correlogram in order
ggcorrplot(cor(data_spa_waters_raw), 
           method = "circle",
           hc.order = TRUE,
           type = "lower",
           lab=TRUE)

#TODO: analysing correlations between compaunds

#--------------------------------------------------------------------------------
# Hamish's part to analyse

# Mahalobis distance for assessing measurements multivariately
(mu.hat <- colMeans(data_spa_waters_raw))
(Sigma.hat <- cov(data_spa_waters_raw))

dM <- mahalanobis(data_spa_waters_raw, center=mu.hat, cov=Sigma.hat)

upper.quantiles <- qchisq(c(.9, .95, .99), df=9)
density.at.quantiles <- dchisq(x=upper.quantiles, df=9)
cut.points <- data.frame(upper.quantiles, density.at.quantiles)

ggplot(data.frame(dM), aes(x=dM)) +
  geom_histogram(aes(y=after_stat(density)), bins=nclass.FD(dM), 
                 fill="white", col="black") +
  geom_rug() +
  stat_function(fun=dchisq, args = list(df=9), 
                col="red", size=2, alpha=.7, xlim=c(0,25)) +
  geom_segment(data=cut.points, 
               aes(x=upper.quantiles, xend=upper.quantiles, 
                   y=rep(0,3), yend=density.at.quantiles),
               col="blue", size=2) +
  xlab("Mahalanobis distances and cut points") +
  ylab("Histogram and density")

# Assessing measurements
data_spa_waters_raw$dM <- dM
data_spa_waters_raw$surprise <- cut(data_spa_waters_raw$dM,
                       breaks= c(0, upper.quantiles, Inf),
                       labels=c("Typical", "Somewhat", "Surprising", "Very"))
table(data_spa_waters_raw$surprise)# Pairs plot

#Pairs plot
ggpairs(data_spa_waters_raw, columns=1:9, 
        ggplot2::aes(col=surprise, alpha=.5),
        upper = list(continuous = "density", combo = "box_no_facet")) +
  ggplot2::scale_color_manual(values=c("lightgray", "green", "blue", "red")) +
  ggplot2::theme(axis.text.x = element_text(angle=90, hjust=1))

#------------------------------------------------------------------------------------

# Analyzing the data with classes

# Melting the data
data_spa_waters_classes.melt <- melt(data=data_spa_waters_classes[-ncolumns_data_spa_waters_classes], 
                  id.vars = "Geological.structure", 
                  variable.name = "Chemical_Compound")

# Box plots with notches
ggplot(data=data_spa_waters_classes.melt, aes(x=Chemical_Compound, y=value)) +
  geom_boxplot(aes(col=Chemical_Compound), notch = TRUE) +
  facet_grid(~ Geological.structure) +
  ylab("Measurements") +
  theme(axis.text.x = element_text(size=7, angle=90, hjust=1),
        legend.position = "none")

# Notice, the data doesn't have enough samples for the notches to appear nicely

# Box plots log scale
ggplot(data=data_spa_waters_classes.melt, aes(x=Chemical_Compound, y=value)) +
  geom_boxplot(aes(col=Chemical_Compound), notch = FALSE) +
  scale_y_log10() +
  facet_grid(. ~ Geological.structure) +
  ylab("Measurements") + 
  theme(axis.text.x = element_text(size=7, angle=90, hjust=1),
        legend.position = "none")

# TODO: Analyse box plots

#PCA
# Set data_spa_waters_raw back to the original:
data_spa_waters_raw <- data_spa_waters_classes[,-c(ncolumns_data_spa_waters_classes-1, ncolumns_data_spa_waters_classes)]

# Calculate PCA
PCA.data_spa_waters_raw <- prcomp(data_spa_waters_raw, center=TRUE, scale=TRUE)
summary(PCA.data_spa_waters_raw)

# Plot proportion of variances of the PCAs (importance of each principal component)
plot(PCA.data_spa_waters_raw, type="l")
# TODO: Say some words about it

options(digits = 4) # Set to show 4 digits
PCA.data_spa_waters_raw$center  # Means of the data
PCA.data_spa_waters_raw$sdev  # Standard deviations of the data
# eigenvectors
PCA.data_spa_waters_raw$rotation # PC contributions of each variable

#TODO: Some words about it

# Showing the biplot PC 1-2
# Pairs plot with unit circle
ggbiplot(PCA.data_spa_waters_raw, circle = TRUE)

# Paris plot with groups of Geological structure
ggbiplot(PCA.data_spa_waters_raw, obs.scale = 1, var.scale = 1,
         groups = interaction(data_spa_waters_classes$Geological.structure), ellipse = TRUE)

# Paris plot with groups of Region
ggbiplot(PCA.data_spa_waters_raw, obs.scale = 1, var.scale = 1,
         groups = interaction(data_spa_waters_classes$Region), ellipse = TRUE)

summary(PCA.data_spa_waters_raw$x)

# Showing the biplot PC 1-3
ggbiplot(PCA.data_spa_waters_raw, obs.scale = 1, var.scale = 1, circle = TRUE, choices = c(1, 3))

# TODO: Analyse the pair plots

print(summary(PCA.data_spa_waters_raw)$importance[, 1:4])

#------------------------------------------------------------------------
# Tried 2 types of factor analysis
# Factor analysis with factanal
n = 4
pvalues <- rep(0, n)

for(f in 1:n)
  pvalues[f] <- factanal(data_spa_waters_raw, factors = f)$PVAL
ggplot(data.frame(pvalues), aes(x=1:n, y=pvalues)) +
  geom_point(size=3) +
  xlab("Number of factors")

factanal(data_spa_waters_raw, factors=4)

fa.diagram(factanal(data_spa_waters_raw, factors=4)$loadings, digits=3)

fa_result <- factanal(data_spa_waters_raw, factors = 4, rotation = "varimax", scores = "regression")
fa_result$scores  # NULL

ggplot(fa_result$score, aes(x = Factor1, y = Factor2)) +
  geom_point()



# Factor Analysis based on Principal Components (Varimax rotation)

# Step 1: Varimax-rotated factor analysis, 4 factors
fa_result <- principal(data_spa_waters_raw, nfactors = 4, rotate = "varimax")

# Check variance accounted for by each rotated factor
fa_result$Vaccounted

fa_result$loadings 




# 1. Build scores data frame from your FA/PCA result
scores_df <- as.data.frame(PCA.data_spa_waters_raw$x)   # or PCA.data_cube$x

# 2. Attach BOTH categorical variables
scores_df$Region <- data_spa_waters_classes$Region
scores_df$Geological.structure <- data_spa_waters_classes$Geological.structure
scores_df$Label <- 1:nrow(scores_df)

# 3. Dynamic color/shape mapping for Region
regions <- unique(scores_df$Region)
my_colors <- c("Pannonian Basin" = "magenta",
               "Serbo-Macedonian massif" = "blue",
               "Carpatho-Balkans" = "black",
               "Vardar Zone" = "red")

my_shapes <- c("Pannonian Basin" = 8,
               "Serbo-Macedonian massif" = 16,
               "Carpatho-Balkans" = 0,
               "Vardar Zone" = 2)

# 4. Plot: points colored/shaped by Region, ellipses drawn around Geological.structure

library(ggnewscale)

ggplot(scores_df, aes(x = PC1, y = PC2 + PC3 + PC4)) +
  stat_ellipse(aes(group = Geological.structure, color = Geological.structure), 
               type = "norm", level = 0.95, 
               linetype = 2, linewidth = 0.8) +
  scale_color_manual(values = c("Volanogenic massif" = "blue", "Cluster II" = "red", "Hydrogeological basin" = "magenta")) +  # adjust names/colors to your actual Geological.structure categories
  new_scale_color() +  # reset color scale so Region gets its own
  geom_point(aes(color = Region, shape = Region), size = 3, stroke = 1.2) +
  geom_text(aes(label = Label, color = Region), vjust = -0.8, hjust = 0.5, 
            size = 3.5, fontface = "bold", show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(x = "PC 1", y = "PC 2,3,4") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.85),
        legend.background = element_rect(color = "black", linewidth = 0.3))


#-------------------------------------------------------------------------------------






#Scores plot:

# Build the scores data frame
scores_df <- as.data.frame(PCA.data_spa_waters_raw$x)
scores_df$Group <- data_spa_waters_classes$Geological.structure
scores_df$Label <- 1:nrow(scores_df)

# Variance explained for axis labels
var_explained <- (PCA.data_spa_waters_raw$sdev^2 / sum(PCA.data_spa_waters_raw$sdev^2)) * 100

# Dynamic color/shape mapping (adjust palette length if you have more groups)
groups <- unique(scores_df$Group)
my_colors <- setNames(c("magenta", "blue", "black", "red", "darkgreen", "orange")[1:length(groups)], groups)
my_shapes <- setNames(c(8, 16, 0, 2, 17, 15)[1:length(groups)], groups)

ggplot(scores_df, aes(x = PC1, y = PC2, color = Group)) +
  stat_ellipse(type = "norm", level = 0.95, linetype = 1, linewidth = 0.8) +
  geom_point(aes(shape = Group), size = 3, stroke = 1.2) +
  geom_text(aes(label = Label), vjust = -0.8, hjust = 0.5, size = 3.5, 
            fontface = "bold", show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(x = paste0("PC1 : ", round(var_explained[1], 2), "%"),
       y = paste0("PC2 : ", round(var_explained[2], 2), "%")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.2, 0.85),
        legend.background = element_rect(color = "black", linewidth = 0.3))


#Scores plot:

# Build the scores data frame
scores_df <- as.data.frame(PCA.data_spa_waters_raw$x)
scores_df$Group <- data_spa_waters_classes$Geological.structure
scores_df$Label <- 1:nrow(scores_df)

# Variance explained for axis labels
var_explained <- (PCA.data_spa_waters_raw$sdev^2 / sum(PCA.data_spa_waters_raw$sdev^2)) * 100

# Dynamic color/shape mapping (adjust palette length if you have more groups)
groups <- unique(scores_df$Group)
my_colors <- setNames(c("magenta", "blue", "black", "red", "darkgreen", "orange")[1:length(groups)], groups)
my_shapes <- setNames(c(8, 16, 0, 2, 17, 15)[1:length(groups)], groups)

ggplot(scores_df, aes(x = PC1, y = PC2, color = Group)) +
  stat_ellipse(type = "norm", level = 0.95, linetype = 1, linewidth = 0.8) +
  geom_point(aes(shape = Group), size = 3, stroke = 1.2) +
  geom_text(aes(label = Label), vjust = -0.8, hjust = 0.5, size = 3.5, 
            fontface = "bold", show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(x = paste0("Factor 1 (PC1) : ", round(var_explained[1], 2), "%"),
       y = paste0("Factor 2 (PC2) : ", round(var_explained[2], 2), "%")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.2, 0.85),
        legend.background = element_rect(color = "black", linewidth = 0.3))


#Clustering:

library(ggplot2)

# 1. Build the scores data frame
scores_df <- as.data.frame(PCA.data_spa_waters_raw$x)
scores_df$Group <- data_spa_waters_classes$Geological.structure
scores_df$Label <- 1:nrow(scores_df)

# 2. Variance explained for axis labels
var_explained <- (PCA.data_spa_waters_raw$sdev^2 / sum(PCA.data_spa_waters_raw$sdev^2)) * 100

# 3. Run clustering on the PC scores (choose number of clusters, e.g. 3)
set.seed(123)
k <- 3  # adjust based on your data
km <- kmeans(scores_df[, c("PC1", "PC2")], centers = k)
scores_df$Cluster <- factor(km$cluster)

# 4. Dynamic color/shape mapping for Group
groups <- unique(scores_df$Group)
my_colors <- setNames(c("magenta", "blue", "black", "red", "darkgreen", "orange")[1:length(groups)], groups)
my_shapes <- setNames(c(8, 16, 0, 2, 17, 15)[1:length(groups)], groups)

ggplot(scores_df, aes(x = PC1, y = PC2)) +
  stat_ellipse(aes(group = Cluster), type = "norm", level = 0.95, 
               color = "grey30", linetype = 2, linewidth = 0.8) +
  geom_point(aes(color = Group, shape = Group), size = 3, stroke = 1.2) +
  geom_text(aes(label = Label, color = Group), vjust = -0.8, hjust = 0.5, 
            size = 3.5, fontface = "bold", show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(x = paste0("Factor 1 (PC1) : ", round(var_explained[1], 2), "%"),
       y = paste0("Factor 2 (PC2) : ", round(var_explained[2], 2), "%")) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.2, 0.85),
        legend.background = element_rect(color = "black", linewidth = 0.3))

