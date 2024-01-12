# 1 --------------------------------------------------------
# setting
library(prabclus)
data(tetragonula)
ta <- alleleconvert(strmatrix=tetragonula)
tai <- alleleinit(allelematrix=ta)
tai$distmat
plot(tai$amatrix)

# a
library(MASS)
library(cluster)

mds <- cmdscale(tai$distmat, k = 2)
plot(mds, xlab = "Dimension 1", 
     ylab = "Dimension 2", main = "MDS plot",
     labels = rownames(mds), cex = 0.6)

# Hierarchical clustering
hc <- list()
hc$ward <- hclust(as.dist(tai$distmat), method = "ward.D")
plot(hc$ward, main = "Dendrogram for hierarchical clustering - Ward's method")
abline(h = 2, col = "red", lty = 2) # 10
hc$comp <- hclust(as.dist(tai$distmat), method = "complete")
plot(hc$comp, main = "Dendrogram for hierarchical clustering - Complete linkage")
abline(h = 0.75, col = "red", lty = 2) # 11
hc$sing <- hclust(as.dist(tai$distmat), method = "single")
plot(hc$sing, main = "Dendrogram for hierarchical clustering - Single linkage")
abline(h = 0.48, col = "red", lty = 2) # 12
hc$ave <- hclust(as.dist(tai$distmat), method = "average")
plot(hc$ave, main = "Dendrogram for hierarchical clustering - Average")
abline(h = 0.67, col = "red", lty = 2) # 9

# how are clusters' shapes?
plot(mds, xlab = "Dimension 1", 
     ylab = "Dimension 2", main = "MDS plot",
     labels = rownames(mds), cex = 0.6)
# are there outliers in data?
boxplot(tai$distmat)
# yes there are a lot
# are data spherical? 
pca <- prcomp(tai$distmat)
biplot(pca)
library(cluster)
kmeans_fit <- kmeans(tai$distmat, centers = 3)
sil <- silhouette(kmeans_fit$cluster, dist(tai$distmat))
plot(sil)
# yes they are enough spherical

###############################################################################################
# we conclude that considering the potential presence of noise and outliers in genetic data,  # 
# Ward's or Average linkage might be more suitable as they are more robust to noise           #
# and outliers.                                                                               #
###############################################################################################

# Elbow method
a <- numeric(20)
for (i in 2:21) {
  hc_cut <- cutree(hc$ward, k = i)
  a[i - 1] <- sum((tai$distmat)^2 * (as.numeric(factor(hc_cut)) == 1))
}
plot(2:21, wss, type = "b", xlab = "Number of clusters", ylab = "Within groups sum of squares")
wss <- numeric(20)
for (i in 2:21) {
  hc_cut <- cutree(hc$ave, k = i)
  wss[i - 1] <- sum((tai$distmat)^2 * (as.numeric(factor(hc_cut)) == 1))
}
plot(2:21, wss, type = "b", xlab = "Number of clusters", ylab = "Within groups sum of squares")

##############################################################################################
# both methods give us 3 clusters as the best option, even though looking at the dendrograms #
# i would have choosen an higher nunber                                                      #
##############################################################################################

clus <- cutree(hc$ward, 3)
plot(mds, col = clus, pch = clus, main = "Clustering using HC")

# PAM
# ASW to find the right K
pasw <- NA
pclusk <- list()
psil <- list()
# Look at K between 2 and 30:
for (k in 2:30){
  # PAM clustering:
  pclusk[[k]] <- pam(tai$distmat,k)
  # Computation of silhouettes:
  psil[[k]] <- silhouette(pclusk[[k]],dist=tai$distmat)
  # ASW needs to be extracted:
  pasw[k] <- summary(psil[[k]])$avg.width
}
# Plot the ASW-values against K:
plot(1:30,pasw,type="l",xlab="Number of clusters",ylab="ASW", main = "ASW")
abline(v = c(1:30), lty = 2, col = "grey")
# 7 cluster seems to be optimal
taipam <- pam(tai$distmat,7)
plot(mds, col = taipam$clustering, pch = taipam$clustering, main = "clustering using PAM")

#################################################################################################
# PAM method seems to capture better the divisions between clusters, even with a non irrelevant #
# uncertainity, maybe due to the fact that the two dimensions offered by MDS are not enough to  #
# distiguish them completely                                                                    #
#################################################################################################

# b
set.seed(123)
# K-means
gapnc <- function(data,FUNcluster=kmeans,
                  K.max=10, B = 100, d.power = 2,
                  spaceH0 ="scaledPCA",
                  method ="globalSEmax", SE.factor = 2,...){
  # As in original clusGap function the ... arguments are passed on
  # to the clustering method FUNcluster (kmeans).
  # Run clusGap
  gap1 <- clusGap(data,kmeans,K.max, B, d.power,spaceH0,...)
  # Find optimal number of clusters; note that the method for
  # finding the optimum and the SE.factor q need to be specified here.
  nc <- maxSE(gap1$Tab[,3],gap1$Tab[,4],method, SE.factor)
  # Re-run kmeans with optimal nc.
  kmopt <- kmeans(data,nc,...)
  out <- list()
  out$gapout <- gap1
  out$nc <- nc
  out$kmopt <- kmopt
  out
}

kmeansres <- gapnc(mds, K.max = 20, B = 200)
kmeansres$nc

vec <- vector()
for (i in 1:100) {
  gap <- gapnc(mds, spaceH0 = "scaledPCA", SE.factor = 2, K.max = 20)
  vec[i] <- gap$nc
}
nc <- ceiling(mean(vec)) # 9
kmeanstai <- kmeans(mds, 9, nstart = 200)

plot(mds, col = kmeanstai$cluster, pch = kmeanstai$cluster, 
     main = "Clustering using K-means")

# Mixture model
library(mclust)
mix <- Mclust(mds,G=1:15)
summary(mix)
summary(mix$BIC)
plot(mix)
m5 <- Mclust(mds,G=5)
plot(m5)
pca <- princomp(mds)
plot(pca$scores,col=mix$classification,
     pch=mix$classification)

# K-means and Gaussian mixture models assume that the data can be partitioned into clusters. This can be 
# good when clusters have a simple and relatively well distigueed structure.
# They also can be more efficient for larger datasets, as they do not require the computation of a full
# distance matrix, which is often needed in hierarchical clustering. This can make them more practical 
# and feasible for analyzing large genetic datasets.
# in this case I think that we could be in the middle, we have almost but not completely clear datasets
# observing the data with MDS, and even if we know that we are looking for groups of bees, we don't know 
# how many.
# So, hierarchical clustering might be more appropriate due to its ability to capture different cluster 
# shapes and its flexibility in determining the number of clusters without knowing them before.

# Do you think it would be better, for this task, to produce 
# an MDS solution with p > 2?

# 2 ---------------------------------------------------------------------------------------
wdbc <- read.csv("wdbc.txt",header=FALSE)
# The variables I ask you to use are variables 3 to 12,
# so wbbcc is the data set to be clustered:
wdbcc <- wdbc[,3:12]
# There is also a diagnosis whether cancer is benign or malign as variable 2
# in the data set:
wdbcdiag <- as.factor(as.integer(as.factor(wdbc[,2])))
str(wdbcc)

# Clusterings
# Gaussian mixture model (numbers of clusters up to 10)
library(mclust)
mix <- Mclust(wdbcc,G=1:10)
summary(mix)
summary(mix$BIC)
plot(mix)
molive3 <- Mclust(wdbcc,G=4)
plot(molive3)
sprolive <- princomp(wdbcc)
plot(sprolive$scores,col=mix$classification,
     pch=mix$classification)

# Hierarchical clustering
# Gower dissimilarity:
library(cluster)
gdiss <- daisy(wdbcc,metric="gower")
hc <- list()
hc$ward <- hclust(gdiss, method = "ward.D")
plot(hc$ward, main = "Dendrogram for hierarchical clustering - Ward's method")
abline(h = 2.5, col = "red", lty = 2) # 6
hc$comp <- hclust(gdiss, method = "complete")
plot(hc$comp, main = "Dendrogram for hierarchical clustering - Complete linkage")
abline(h = 0.37, col = "red", lty = 2) # 4
hc$sing <- hclust(gdiss, method = "single")
plot(hc$sing, main = "Dendrogram for hierarchical clustering - Single linkage")
abline(h = 0.48, col = "red", lty = 2) # 12
hc$ave <- hclust(gdiss, method = "average")
plot(hc$ave, main = "Dendrogram for hierarchical clustering - Average")
abline(h = 0.2, col = "red", lty = 2) # 5
# decide the best number of clustering justifying why one method is better than the others
clus <- cutree(hc$ward, 6)

library(MASS)

mds <- cmdscale(gdiss, k = 2)
plot(mds, xlab = "Dimension 1", 
     ylab = "Dimension 2", main = "MDS plot",
     labels = rownames(mds), cex = 0.6)
plot(mds, col = clus, pch = clus, main = "Clustering using HC")


# K-means clustering
set.seed(123)
# K-means
gapnc <- function(data,FUNcluster=kmeans,
                  K.max=10, B = 100, d.power = 2,
                  spaceH0 ="scaledPCA",
                  method ="globalSEmax", SE.factor = 2,...){
  # As in original clusGap function the ... arguments are passed on
  # to the clustering method FUNcluster (kmeans).
  # Run clusGap
  gap1 <- clusGap(data,kmeans,K.max, B, d.power,spaceH0,...)
  # Find optimal number of clusters; note that the method for
  # finding the optimum and the SE.factor q need to be specified here.
  nc <- maxSE(gap1$Tab[,3],gap1$Tab[,4],method, SE.factor)
  # Re-run kmeans with optimal nc.
  kmopt <- kmeans(data,nc,...)
  out <- list()
  out$gapout <- gap1
  out$nc <- nc
  out$kmopt <- kmopt
  out
}

kmeansres <- gapnc(wdbcc, K.max = 20, B = 200)
kmeansres$nc # 1

vec <- vector()
for (i in 1:100) {
  gap <- gapnc(wdbcc, spaceH0 = "scaledPCA", SE.factor = 2, K.max = 20)
  vec[i] <- gap$nc
}
nc <- ceiling(mean(vec)) # 1, so 2
kmeanstai <- kmeans(mds, 2, nstart = 200)

plot(wdbcc, col = kmeanstai$cluster, pch = kmeanstai$cluster, 
     main = "Clustering using K-means")
plot(mds, col = kmeanstai$cluster, pch = kmeanstai$cluster, 
     main = "Clustering using K-means")


# PAM
# ASW to find the right K
pasw <- NA
pclusk <- list()
psil <- list()
# Look at K between 2 and 30:
for (k in 2:30){
  # PAM clustering:
  pclusk[[k]] <- pam(gdiss,k)
  # Computation of silhouettes:
  psil[[k]] <- silhouette(pclusk[[k]],dist=gdiss)
  # ASW needs to be extracted:
  pasw[k] <- summary(psil[[k]])$avg.width
}
# Plot the ASW-values against K:
plot(1:30,pasw,type="l",xlab="Number of clusters",ylab="ASW", main = "ASW")
# 2 cluster seems to be optimal
taipam <- pam(gdiss,2)

plot(wdbcc, col = taipam$clustering, pch = taipam$clustering, main = "clustering using PAM")
plot(mds, col = taipam$clustering, pch = taipam$clustering, main = "clustering using PAM")

########################################################################################################
# Comparing                                                                                            #
# hierarchical clustering might be more appropriate because we don't know how many clusters we have to #
# find initially, it's also true that we have a very large dataset, we are on the limit to not consider#
# hierarchical clustering only because it would be computationally difficult to make the distances     #
# matrix.                                                                                              #

# Comparing to wdbcdiag using ARI
adjustedRandIndex(taipam$clustering, wdbcdiag) # 1
adjustedRandIndex(kmeanstai$cluster, wdbcdiag) # 2
adjustedRandIndex(clus, wdbcdiag) # 4
adjustedRandIndex(molive3$classification, wdbcdiag) # 3
########################################################################################################