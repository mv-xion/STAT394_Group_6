# packages
install.packages("rstatix")
library(rstatix)
library(MASS)
library(car)

#----------------------------------------------------------------
# Data preparation 
# Reading the data
data_spa_waters <- read.csv("Desktop/Hydrochemical Data.csv") # Modify this one
head(data_spa_waters) # Showing data


# Cleaning data
# Removing the first 2 columns (number, place name) 
data_spa_waters_classes = data_spa_waters[-1, -c(1, 2)]
ncolumns_data_spa_waters_classes = ncol(data_spa_waters_classes)

# Consider data as numeric, round to 3 decimals except class column
data_spa_waters_classes[-c(ncolumns_data_spa_waters_classes-1, ncolumns_data_spa_waters_classes)] = lapply(data_spa_waters_classes[-c(ncolumns_data_spa_waters_classes-1, ncolumns_data_spa_waters_classes)], function(x) round(as.numeric(x), 3))
# Set Structure and region as factors
data_spa_waters_classes$ Geological.structure <- as.factor(data_spa_waters_classes$ Geological.structure)


#Show the data
head(data_spa_waters_classes)
summary(data_spa_waters_classes)

# Assessing measurements
data_spa_waters_raw$dM <- dM
data_spa_waters_raw$surprise <- cut(data_spa_waters_raw$dM,
                                    breaks= c(0, upper.quantiles, Inf),
                                    labels=c("Typical", "Somewhat", "Surprising", "Very"))
table(data_spa_waters_raw$surprise)# Pairs plot


# pairs transformed
library(GGally)
ggpairs(data_cube, columns=1:9, 
        ggplot2::aes(col = data_cube$region, alpha=.5),
        upper = list(continuous = "density", combo = "box_no_facet")) +
  ggplot2::scale_color_manual(values=c("black", "green", "blue", "red")) +
  ggplot2::theme(axis.text.x = element_text(angle=90, hjust=1))

#data transformation prep
shapiro_test(data_cube, vars = )
box.data <- powerTransform(data_cube)

summary(box.data)

# data transformation
# Natural Log Transformations (Rounded Pwr = 0.00)
log_vars <- c("T", "EC", "TS", "Ca2.", "Cl.", "SO2.4", "HCO.3", "SiO2")
data_cube[log_vars] <- lapply(data_cube[log_vars], log)

# Cube Root Transformations (Rounded Pwr = 0.33)
cube_vars <- c("Na.", "K.")
data_cube[cube_vars] <- lapply(data_cube[cube_vars], function(x) x ^ (1/3))

# Custom Power Transformation (Rounded Pwr = 0.16)
data_cube$Mg2. <- data_cube$Mg2. ^ 0.16


data_cube$region <- factor(rep(c(
  "Hydro", "Karstic", "Volcano", "Meta"), times = c(5, 5, 14, 6)))

data_cube

num <- sapply(data_cube, is.numeric)

# Mahalobis distance
(mu.hat <- colMeans(data_cube[num]))
(Sigma.hat <- cov(data_cube[num]))

dM <- mahalanobis(data_cube[num], center=mu.hat, cov=Sigma.hat)

upper.quantiles <- qchisq(c(.9, .95, .99), df=9)
density.at.quantiles <- dchisq(x=upper.quantiles, df=9)
cut.points <- data.frame(upper.quantiles, density.at.quantiles)

ggplot(data.frame(dM), aes(x=dM)) +
  geom_histogram(aes(y=after_stat(density)), bins=10, 
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
data_cube$dM <- dM
data_cube$surprise <- cut(data_cube$dM,
                          breaks= c(0, upper.quantiles, Inf),
                          labels=c("Typical", "Somewhat", "Surprising", "Very"))

# Investigating surprising observations.
table(data_cube$surprise)
outliers <- which(data_cube$dM > 14.68366)
outliers

