# 1 -----------------------------------------------------------------------------------
data <- read.csv("C:/Users/User/Desktop/SCUOLA/UNI/2. Magistrale/II ANNO/primo semestre/modern statistics and big data analytics/slide/datasets/stars5000.txt", sep="")
library(teigen)
library(mixsmsn)
library(cluster)
library(smacof)
library(mclust)
library("robustbase")

plot(data)
str(data)

data300 <- data[sample(5000,300),]
summary(data300)
set.seed(776655)
pr300 <- princomp(data300)
plot(pr300$scores, main = "PCA plot")
mds300 <- cmdscale(dist(data300))
# Assuming mds2000 is your MDS plot object
plot(mds300, main = "MDS plot")
# it seems that a t distribution could fit well, or maybe a t-skewed
hist(data$casn)
lines(density(data$casn), col = "blue", lwd = 2)
points(data$casn, pch = 20)
hist(data$cacont)
lines(density(data$cacont), col = "blue", lwd = 2)
points(data$cacont, pch = 20)
hist(data$kl1)
lines(density(data$kl1), col = "blue", lwd = 2)
points(data$kl1, pch = 20)
hist(data$kl2)
lines(density(data$kl2), col = "blue", lwd = 2)
points(data$kl2, pch = 20)
hist(data$xh1)
lines(density(data$xh1), col = "blue", lwd = 2)
points(data$xh1, pch = 20)
hist(data$xh2)
lines(density(data$xh2), col = "blue", lwd = 2)
points(data$xh2, pch = 20)
plot(data$casn, pch = 16, col = ifelse(data$casn %in% boxplot.stats(data$casn)$out, "red", "blue"), 
     main = "Scatter Plot with Outliers")
plot(data$cacont, pch = 16, col = ifelse(data$cacont %in% boxplot.stats(data$cacont)$out, "red", "blue"), 
     main = "Scatter Plot with Outliers")
plot(data$kl1, pch = 16, col = ifelse(data$kl1 %in% boxplot.stats(data$kl1)$out, "red", "blue"), 
     main = "Scatter Plot with Outliers")
plot(data$kl2, pch = 16, col = ifelse(data$kl2 %in% boxplot.stats(data$kl2)$out, "red", "blue"), 
     main = "Scatter Plot with Outliers")
plot(data$xh1, pch = 16, col = ifelse(data$xh1 %in% boxplot.stats(data$xh1)$out, "red", "blue"), 
     main = "Scatter Plot with Outliers")
plot(data$xh2, pch = 16, col = ifelse(data$xh2 %in% boxplot.stats(data$xh2)$out, "red", "blue"), 
     main = "Scatter Plot with Outliers")

n <-  Mclust(data300, G=3:5)
summary(n)
summary(n$BIC)
n5 <- Mclust(data300, G = 5)
summary(n5)
n.bic <- n5$bic
plot(pr300$scores, cex = 2, col = n5$classification+1, pch = 20)

ns <- smsn.search(data300,nu=1,g.min=3,g.max=5,family="Skew.normal", criteria = "bic")
summary(ns)
ns.bic <- ns$best.model$bic
plot(mds300, col = ns$best.model$group+1, cex = 2, pch = 20)
plot(pr300$scores, col = ns$best.model$group+1, cex = 2, pch = 20)

t <- teigen(data300, Gs=3:5, scale = F)
summary(t)
t.bic <- t$bic
plot(pr300$scores, col = t$classification+1, cex = 2, pch = 20)
plot(t,xmarg=1,ymarg=5,what="contour")
plot(t,xmarg=1,ymarg=5,what="uncertainty")

ts <- smsn.search(data300,nu=4,g.min=3,g.max=5,family="Skew.t")
summary(ts)
ts.bic <- ts$best.model$bic
plot(pr300$scores, col = ts$best.model$group+1, pch = 20, cex = 2)


# Most convincing clustering and why. 
res1 <- rbind(paste("Gmm:",abs(round(n.bic,2)), paste("Skewed Normal:", round(ns$best.model$bic,2)), 
      paste("T-student:",round(abs(t$bic),2)), paste("T-skewed:", round(ts$best.model$bic,2))))
# Comparison with standardized data
library(pcaPP)
s <- ScaleAdv(data300, center = median, scale = mad)
s <- s$x
stand <- scale(data300)

n.stand <- Mclust(stand, G=3:5)
n.stand.bic <- n.stand$bic
ns.stand <- smsn.search(s,nu=1,g.min=2,g.max=5,family="Skew.normal")
ns.stand.bic <- ns.stand$best.model$bic
t.stand <- teigen(s, Gs=3:5, scale = FALSE)
t.stand.bic <- t.stand$bic
ts.stand <- smsn.search(s, nu = 2, g.min = 3, g.max = 5, family = "Skew.t")

res2 <- rbind(paste("Gmm:",abs(round(n.stand.bic,2)), paste("Skewed Normal:", round(ns.stand$best.model$bic,2)), 
                    paste("T-student:",round(abs(t.stand$bic),2)), paste("T-skewed:", "result not found")))
res2

results <- rbind(res1,res2)
results

# I had to reduce dimentionality and boundaries in the functions because in several cases they gave me problems
# and algorithms didn't converge.
# It seems that the moste performing model, both for standardized and original data is the Gaussian Mixture.
# We also notice how performance heavily increases standardizing data.

# 2 ------------------------------------------------------------------------------------------------
data <- read.csv("C:/Users/User/Desktop/SCUOLA/UNI/2. Magistrale/II ANNO/primo semestre/modern statistics and big data analytics/slide/datasets/stars5000.txt", sep="")

# time evaluation
speedyclust <- function(data, n.reduced, g.max){
  n.red <- list()
  n.predict <- list()
  
  data.reduced <- data[sample(nrow(data),n.reduced),]
  n.red <- Mclust(data.reduced, G = 1:g.max)
  n.predict <- predict.Mclust(n.red, data)
  return(n.red)
}

time1 <- system.time(n.speedy <- speedyclust(data, n.reduced = 1000, g.max = 6))
time2 <- system.time(n <- Mclust(data, G = 1:6))
time<- cbind(time1, time2)[3,]
time.gained <- paste("We gained", round(time[2] - time[1], 2), "seconds") 
time.gained

# results evaluation
n <-  Mclust(data, G=1:6)
n.speedy <- speedyclust(data, n.reduced = 1000, g.max = 6)

plot()
data1000 <- data[sample(5000,1000),]
pr1000 <- princomp(data1000)
pr <- princomp(data)
layout(c(1,2))
plot(pr1000$scores, col = n.speedy$classification, pch = 20, cex = 2, 
     main = "Reduced dataset's classification")
plot(pr$scores, col = n$classification, pch = 20, cex = 2, 
     main = "Full dataset classification")
# full dataset seems to be way better fitted than the reduced one
# let's check with BIC measures
rbind(paste("Full dataset", abs(round(n$bic,2))), 
      paste("Reduced dataset", abs(round(n.speedy$bic,2))))
layout(c(1,1))
plot(n)
plot(n.speedy)
# it seems that they have proportionally similar results, but I think that comparing BIC of models
# that use different datasets could not be theoretically correct.
# So, relying on the plots of PCA, I'd say that even though we gained time, this method doesn't give us better
# results, in this specific case.

# 3 ----------------------------------------------------------------------------------------------
# The number of free parameters in mixture models can be calculated based on the number of parameters associated 
# with the mixture components. For a Gaussian mixture model (GMM), skew-normal mixture, and mixture of multivariate t 
# distributions, the number of parameters depends on the choice of covariance matrices and other constraints. 
# Let's break down each case:
# 
# (a) Gaussian Mixture Model (GMM) with Fully Flexible Covariance Matrices:
# Each component in a GMM has parameters for the mean, covariance matrix, and mixing weight.
# 
# Mean Parameters: 10 variables × 4 components = 40 parameters.
# Covariance Matrix Parameters: 10×(10+1)/2 (for a fully flexible covariance matrix) × 4 components = 220 parameters.
# Mixing Weights Parameters: 4???1 (because they sum to 1) = 3 parameters.
# Total parameters for GMM: 40+220+3=263 parameters.
# 
# (b) Gaussian Mixture Model (GMM) Assuming Spherical Covariance Matrices:
# Each component now has a single variance parameter instead of a full covariance matrix.
# 
# Mean Parameters: 40 parameters.
# Variance Parameters: 10 variables x 4 components = 40 parameters.
# Mixing Weights Parameters: 3 parameters.
# Total parameters for GMM with spherical covariance matrices: 40+40+3=83 parameters.
# 
# (c) Fully Flexible Skew-Normal Mixture:
# In a skew-normal distribution, you have mean, scale, and shape parameters for each component.
# 
# Mean Parameters: 40 parameters.
# Scale Parameters: 10×(10+1)/2 x 4 components = 220 parameters.
# Shape Parameters: 4 components = 4 parameters.
# Mixing Weights Parameters: 3 parameters.
# Total parameters for skew-normal mixture: 40+220+4+3=267 parameters.
# 
# (d) Fully Flexible Mixture of Multivariate t Distributions:
# Similar to the Gaussian case but with an additional parameter for the degrees of freedom in the t-distribution.
# 
# Mean Parameters: 40 parameters.
# Scale Parameters: 10×(10+1)/2 x 4 components = 220 parameters.
# Degrees of Freedom Parameters: 4 components = 4 parameters.
# Mixing Weights Parameters: 3 parameters.
# Total parameters for mixture of t distributions: 40+220+4+3=267 parameters.
# 
# (e) Mixture of Skew-t Distributions with Equal Parameters:
# If skewness, degrees of freedom, and covariance matrices are assumed equal across components:
# 
# Mean Parameters: 40 parameters.
# Scale Parameter: 10×(10+1)/2=55 parameters.
# Mixing Weights Parameters: 3 parameters.
# Total parameters for mixture of skew-t distributions with equal parameters: 40+55+3=98 parameters.

# 4 ---------------------------------------------------------------------------------------------------------------
library(cluster)
library(flexclust)
library(mclust)
library(MASS)
library("smacof")
prodcat <- c("PR","MB","AB","N") # product categories
eventcat <- c("S","M","C","P","N") # event categories
prodp <- eventp <- list()
# Category probabilities within the two mixture components (clusters)
prodp[[1]] <- c(0.4,0.4,0.1,0.1)
prodp[[2]] <- c(0.2,0.1,0.4,0.3)
eventp[[1]] <- c(0.5,0.2,0.1,0.1,0.1)
eventp[[2]] <- c(0.1,0.1,0.1,0.3,0.4)
# The first 400 observations from component 1, then 600 from component 2:
n1 <- 400
n2 <- 600
n <- n1+n2
# Initialisation (provide a data frame of the correct type and size
# without already having the data):
consumers <- data.frame(prod=factor(c(prodcat,rep(NA,n-4)),levels=prodcat),
                        event=factor(c(eventcat,rep(NA,n-5)),levels=eventcat))
# Generation of the data; you may want to set a seed here.
set.seed(2398938)
consumers$prod[1:n1] <- sample(prodcat,n1,replace=TRUE,prob=prodp[[1]])
consumers$event[1:n1] <- sample(eventcat,n1,replace=TRUE,prob=eventp[[1]])
consumers$prod[(n1+1):n] <- sample(prodcat,n2,replace=TRUE,prob=prodp[[2]])
consumers$event[(n1+1):n] <- sample(eventcat,n2,replace=TRUE,prob=eventp[[2]])
# You can run table(consumers) or str(consumers) to see what you got.
table(consumers)

# Estimate a latent class mixture model using the flexmixedruns-function
# including estimating the number of clusters (it is enough to go up to 5).
str(consumers)
dist <- daisy(consumers, metric = "gower")
mds <- cmdscale(dist)

mix.cat <- flexmixedruns(consumers, continuous = 0, discrete = 2, n.cluster=1:5)
mix.cat$optimalk
mix.cat$bicvals
plot(1:5,mix.cat$bicvals,typ="l", xlab="Number of clusters",ylab="BIC") # 2 has minimum BIC
plot(mds,col= mix.cat$flexout[[2]]@cluster, pch=clusym[mix.cat$flexout[[2]]@cluster])
plot(consumers, col = mix.cat$flexout[[2]]@cluster)

# Compare the results with the truth that you know about these data 
# (this can concern the clustering vector, but also the estimated parameters).
p1.1 <- c(mix.cat$flexout[[2]]@components$Comp.1[[1]]@parameters$pp[[1]],0)
p2.1 <- mix.cat$flexout[[2]]@components$Comp.1[[1]]@parameters$pp[[2]]
p1.2 <- c(mix.cat$flexout[[2]]@components$Comp.2[[1]]@parameters$pp[[1]],0)
p2.2 <- mix.cat$flexout[[2]]@components$Comp.2[[1]]@parameters$pp[[2]]
p1 <- data.frame(p1.1,p1.2)
p2 <- data.frame(p2.1,p2.2)
prob2 <- data.frame(eventp[[1]],eventp[[2]])
colnames(prob2) <- c("prob2.1", "prob2.2")
prob1 <- data.frame(rbind(data.frame(prodp[[1]], prodp[[2]]),c(0,0)))
colnames(prob1) <- c("prob1.1", "prob1.2")

result.for.probs <- data.frame(round(cbind(p1,p2),1), cbind(prob1,prob2))
result.for.probs

true_clusters <- c(rep(1, n1), rep(2, n2))
estimated_clusters <- mix.cat$flexout[[2]]@cluster
table(true_clusters, estimated_clusters)

# Estimation seems pretty good.

# Also compute the simple matching distance and cluster the data using any distance
# based clustering method (just one), and compare the clustering with the same number of clusters 
# with the optimal clustering from flexmixedruns. 
#You should find that Average Linkage with Simple Matching is somewhat problematic here. Why is
# this not a good method for these data? 
library(nomclust)
sm <- sm(consumers)
hc2 <- hclust(sm, method = "average")
clust2 <- cutree(hc, 2)
adjustedRandIndex(clust2, estimated_clusters)

# Simple matching distance seems to perform worst. This is probably due to the fact that sm combined with 
# average linkage are not that good at treating outliers, and as we've seen this variables are full of outliers.