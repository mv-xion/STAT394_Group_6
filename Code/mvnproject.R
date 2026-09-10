hydrodata<-read.csv("/Users/jackhardy/Documents/GitHub/STAT394_Group_6/Data/Hydrochemical Data.csv")

summary(hydrodata)

hydroclean<-hydrodata[-1, 3:14]
as.numeric(hydroclean$T)->hydroclean$T
as.numeric(hydroclean$EC)->hydroclean$EC
as.numeric(hydroclean$TS)->hydroclean$TS
as.numeric(hydroclean$Ca2.)->hydroclean$Ca2.

head(hydroclean)
summary(hydroclean)


#Histograms for each
hist(hydroclean$T, breaks = 25)
hist(hydroclean$pH, breaks = 25)
hist(hydroclean$EC, breaks = 25)
#and so on (not very normal)


#QQ plots for each variables
qqnorm(hydroclean$T)
qqline(hydroclean$T)
qqnorm(hydroclean$pH, main = "pH")
qqline(hydroclean$pH)
#T has problem at tails (visible in histograms) and pH is not that great (okay)


#Shapiro-Wilk Test H0: Normal
shapiro.test(hydroclean$T)
shapiro.test(hydroclean$pH)
shapiro.test(hydroclean$EC)
#and so on (these two have failed)
#Note in the paper they applied these tests and transformed some variables using 
#the Box-Cox transformation


#Pair scatter
plot(hydroclean,
     col = "black",
     main = "Pairwise Scatterplots for assessing Bivariate Normality")
#For bivariate normal want scatter plots that are elliptical, smooth and 
#continuous distribution, no obvious clusters or weird shapes, and no major 
#outliers
#These ones aren't this


#Run Mardia's test
#Formal test for data MVN checking whether multiple variables considered 
#together follow MVN distribution
#We want independent observations, don't want extreme outliers, and also want
#reasonable sample size (ours is small, be careful)
#Testing H0: The data follows a MVN distribution
#Looks at two features: skewness whether overall distribution is asymmetric
#and kurtosis whether overall distribution has heavy/light tails
library(MASS) 
library(lavaan)
library(psych) #I don't know which of these it need
psych::mardia(hydroclean)
# 1. MV Skewness with p-value 4.3*10^-20 so reject NULL of MVN based on skewness
# 2. MV Kurtosis with p-value 2.6*10^-6 so reject NULL of MVN based on kurtosis
#So overall Mardia's test gives strong evidence against MVN


#Run Henze-Zirkler test
#Formally tests whether variables together follows a multivariate normal 
#distribution
#NULL the data is MVN
#Data should be numerical, observation should be independent, designed for 
#multiple variables jointly, can be for small data sets but ours is pretty small
library(mvnTest)
HZ.test(hydroclean)
#We get p-value as 0 so reject NULL and so evidence variables jointly don't
#follow MVN distribution 


#Takeaways
# 1. Individual variables (at least the first few) don't follow normal (my
#note from the paper) - Histograms / QQ-plots / Shapiro-Wilk Test
# 2. Pairs scatters are not very elliptical, etc.
# 3. Mardia's fail (however careful with our small data set)
# 4. Henze-Zirkler fail (we should get more confidence seeing as we reject here 
#too)
