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


head(data_Chemical_classes)

summary(data_Chemical_classes)

