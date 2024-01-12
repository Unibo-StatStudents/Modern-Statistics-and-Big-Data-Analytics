# 3 --------------------------------------------------------------------------------------
# Use a simulation to empirically estimate the relative efficiencies (ratios of the variances of 
# the estimators' sample distributions) of the following three estimators of the population mean: 
# (i) arithmetic mean
# (ii) median
# (iii) Huber's M-estimator with eta = 1:5 and based on scale estimated by the MAD 
# (as automatically done by the huberM-function) 
# For the following three distributions: 
# Use 1000 data sets each, with n = 20.
n = 20
set.seed(2894)

# (a) N(0; 1)-distribution
distr1f <- function(n = n, mean = 0, sd = 1){
  distrN <- list()
  
  for (i in 1:1000) {
    distrN[[i]] <- rnorm(n, mean, sd)
  }
  return(distrN)
}
distr1 <- distr1f(n, 0, 1)
# (i) arithmetic mean
meanf <- function(data = listofdata){
  m <- vector()
  
  for (i in 1:length(data)) {
    m[i] <- mean(data[[i]])
  }
  return(m)
}
mean1 <- meanf(dist1)
varmean1 <- var(mean1)
# (ii) median
medianf <- function(data = listofdata){
  med <- vector()
  
  for (i in 1:length(data)) {
    med[i] <- median(data[[i]])
  }
  return(med)
}
median1 <- medianf(dist1)
varmedian1 <- var(median1)
# (iii) Huber's M-estimator with eta = 1:5 and based on scale estimated by the MAD 
library(MASS)
library(robustbase)
huberf <- function(data = listofdata, k = eta){
  hub <- vector()
  
  for (i in 1:length(data)) {
    hub[i] <- huberM(data[[i]], k)$mu
  }
  return(hub)
}
huber1 <- huberf(dist1, k = 1.5)
varhub1 <- var(huber1)

# (b) t2-distribution
distr2f <- function(n = n, df = 2){
  distrt1 <- list()
  
  for (i in 1:1000) {
    distrt1[[i]] <- rt(n, 2)
  }
  return(distrt1)
}
distr2 <- distr2f(n, 2)
# (i) arithmetic mean
mean2 <- meanf(distr2)
varmean2 <- var(mean2)
# (ii) median
median2 <- medianf(distr2)
varmedian2 <- var(median2)
# (iii) Huber's M-estimator with eta = 1:5 and based on scale estimated by the MAD 
huber2 <- huberf(distr2, k = 1.5)
varhub2 <- var(huber2)

# (c) t4-distribution
distr3 <- distr2f(n, 4)
# (i) arithmetic mean
mean3 <- meanf(distr3)
varmean3 <- var(mean3)
# (ii) median
median3 <- medianf(distr3)
varmedian3 <- var(median3)
# (iii) Huber's M-estimator with eta = 1:5 and based on scale estimated by the MAD 
huber3 <- huberf(distr3, k = 1.5)
varhub3 <- var(huber3)

cbind(varmean1/varmean2, varmean1/varmean3, varmean2/varmean3) # 1 - 3 - 2
cbind(varmedian1/varmedian2, varmedian1/varmedian3, varmedian2/varmedian3) # 1 - 3 - 2
cbind(varhub1/varhub2, varhub1/varhub3, varhub2/varhub3) # 1 - 3 - 2

# 4 --------------------------------------------------------------------------------------
# unicef.txt contains information on 121 countries in 1995/96 on the following variables:
# y Child.Mortality: child mortality rate (deaths before the age of 5 per 1000 births)
# X1 Literacy.Fem: literacy among females compared to males in percent
#                  (ie, 100 means that female literacy equals male literacy)
# X2 Literacy.Ad: literacy among adults in percent
# X3 Drinking.Water: percentage of population with access to safe drinking water
# X4 Polio.Vacc: percentage of one-year olds vaccinated against polio
# X5 Tetanus.Vacc.Preg: percentage of pregnant women vaccinated against tetanus
# X6 Urban.Pop: percentage of population living in urban areas
# X7 Foreign.Aid: received foreign aid as percentage of GDP
# It is of interest to explain child mortality from the variables X1-X7 using linear regression.
# The data set has been analysed, using standard least squares linear regression ((lm in R)
# and the robust regression MM-estimator

# Overview of the dataset
library(robustbase)
unicef <- read.table("unicef.txt",header=TRUE)
pairs(unicef,pch=rownames(unicef))
str(unicef)

# Least Squares regression
uniceflm <- lm(Child.Mortality~Literacy.Fem+Literacy.Ad+Drinking.Water+
                 Polio.Vacc+Tetanus.Vacc.Preg+Urban.Pop+Foreign.Aid, data=unicef)
plot(uniceflm,ask=FALSE)

# Robust regression
unicefmm <- lmrob(Child.Mortality~Literacy.Fem+Literacy.Ad+Drinking.Water+
                    Polio.Vacc+Tetanus.Vacc.Preg+Urban.Pop+Foreign.Aid,
                  data=unicef)
plot(unicefmm,which=c(1,2),ask=FALSE)

# Robustness weights
plot(1:121,unicefmm$rweights,xlab="Observation number",ylab="Weight",
     main="MM-estimator, robustness weights",type="n")
text(1:121,unicefmm$rweights,rownames(unicef),cex=0.7)
plot(1:121,unicefmm$init.S$rweights,xlab="Observation number",ylab="Weight",
     main="S-estimator, robustness weights",type="n")
text(1:121,unicefmm$init.S$rweights,rownames(unicef),cex=0.7)

# Residuals vs. fitted
plot(unicefmm$fitted,unicefmm$residuals,xlab="Fitted values",ylab="Residuals",
     main="MM-estimator, residuals vs. fitted")
plot(unicefmm$init.S$fitted,unicefmm$init.S$residuals,xlab="Fitted values",
     ylab="Residuals",main="S-estimator, residuals vs. fitted")

# Finding outliers using covariance matrix of all variables
cunicef <- cov(unicef) # Covariance matrix
mcdunicef <- covMcd(unicef) # MCD with alpha=0.5
mcdu75 <- covMcd(unicef,alpha=0.75) # MCD with alpha=0.75
# Squared robust Mahalanobis distances
plot(1:121,sqrt(mcdunicef$mah),type="n",xlab="Observation",
     ylab="Squared robust Mahalanobis distance",main="MCD with alpha=0.5")
text(1:121,sqrt(mcdunicef$mah),rownames(unicef),cex=0.7)
abline(sqrt(qchisq(0.99,8)),0,col=2)
plot(1:121,sqrt(mcdu75$mah),type="n",xlab="Observation",
     ylab="Squared robust Mahalanobis distance",main="MCD with alpha=0.75")
text(1:121,sqrt(mcdu75$mah),rownames(unicef),cex=0.7)
abline(sqrt(qchisq(0.99,8)),0,col=2)
# Together with standard Mahalanobis distances
plot(sqrt(mahalanobis(unicef,colMeans(unicef),cunicef)),sqrt(mcdunicef$mah),
     type="n",xlab="Squared standard Mahalanobis distance",
     ylab="Squared robust Mahalanobis distance",main="MCD with alpha=0.5")
text(sqrt(mahalanobis(unicef,colMeans(unicef),cunicef)),sqrt(mcdunicef$mah),
     rownames(unicef),cex=0.7)
abline(sqrt(qchisq(0.99,8)),0,col=2)
abline(v=sqrt(qchisq(0.99,8)),col=2)
plot(sqrt(mahalanobis(unicef,colMeans(unicef),cunicef)),sqrt(mcdu75$mah),
     type="n",xlab="Squared standard Mahalanobis distance",
     ylab="Squared robust Mahalanobis distance",main="MCD with alpha=0.75")
text(sqrt(mahalanobis(unicef,colMeans(unicef),cunicef)),sqrt(mcdu75$mah),
     rownames(unicef),cex=0.7)
abline(sqrt(qchisq(0.99,8)),0,col=2)
abline(v=sqrt(qchisq(0.99,8)),col=2)


# There is an obvious outlier in the data set. Create a new data set identical to the original
# one, but with this outlier removed.
# SaoTP is the evident outlier
unicef120 <- unicef[-(which(rownames(unicef)=="SaoTP")),]

# Run both lm and lmrob on the new data set, and compare the results to the results of
# the same regression method on the original data
# Least Squares regression
unicef120lm <- lm(Child.Mortality~Literacy.Fem+Literacy.Ad+Drinking.Water+
                 Polio.Vacc+Tetanus.Vacc.Preg+Urban.Pop+Foreign.Aid, data=unicef120)

plot(unicef120lm,ask=FALSE)

# Robust regression
unicef120mm <- lmrob(Child.Mortality~Literacy.Fem+Literacy.Ad+Drinking.Water+
                    Polio.Vacc+Tetanus.Vacc.Preg+Urban.Pop+Foreign.Aid,
                  data=unicef120)
plot(unicef120mm,which=c(1,2),ask=FALSE)

# (i) Comment on how t-test results of the variables have changed qualitatively 
# (significances and positive or negative direction of influence) 
ls <- summary(uniceflm) #LS
rob <- summary(unicefmm) #Robust
ls120 <- summary(unicef120lm) #LS120
rob120 <- summary(unicef120mm) #Robust120
rbind(ls$coefficients[,3],ls120$coefficients[,3])
rbind(rob$coefficients[,3],rob120$coefficients[,3])

# (ii) Compare the vectors of estimators of the regression coefficients. 
rbind(ls$coefficients[,1],ls120$coefficients[,1])
rbind(rob$coefficients[,1],rob120$coefficients[,1])

# (iii) Comment on what the results mean regarding the sensitivity of the two regressions to
# outliers.
