# Reading the data
data_Chemical <- read.csv("C:\\Users\\ixion\\NZ\\PhD\\STAT394_Group_6\\Data\\Hydrochemical Data.csv")
head(data_Chemical) # Showing data

# Cleaning data
# removing the 1st 2 (number, place name) and the last column as well as the 1st row (unit)
data_Chemical_clean <- data_Chemical[-1, -c(1, 2, ncol(data_Chemical))]
head(data_Chemical_clean)
# Change data types to be all numbers
data_cube <- as.data.frame(lapply(data_Chemical_clean, as.integer))
head(data_cube)

summary(data_cube)

# Correlation of chemical measurements
cor(data_cube)
library(ggplot2)
library(ggcorrplot)
ggcorrplot(cor(data_cube), 
           method = "circle",
           hc.order = TRUE,
           type = "lower")

# Mahalobis distance

(mu.hat <- colMeans(data_cube))
(Sigma.hat <- cov(data_cube))

dM <- mahalanobis(data_cube, center=mu.hat, cov=Sigma.hat)

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

data_cube$dM <- dM
data_cube$surprise <- cut(data_cube$dM,
                       breaks= c(0, upper.quantiles, Inf),
                       labels=c("Typical", "Somewhat", "Surprising", "Very"))
table(data_cube$surprise)

library(GGally)
ggpairs(data_cube, columns=1:9, 
        ggplot2::aes(col=surprise, alpha=.5),
        upper = list(continuous = "density", combo = "box_no_facet")) +
  ggplot2::scale_color_manual(values=c("lightgray", "green", "blue", "red")) +
  ggplot2::theme(axis.text.x = element_text(angle=90, hjust=1))
