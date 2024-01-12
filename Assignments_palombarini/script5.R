# 1 ------------------------------------------------------------------------------------
glass <- read.csv("glass.txt",header=TRUE, sep = " ")
mglass <- as.matrix(glasss)
str(glass)
summary(glass)
library(MASS)
library(cluster)
# standardization useful only for k-means and GMM
glassstd <- data.frame(scale(glass))
mds <- cmdscale(dist(glass, "manhattan"), k = 3)
plot(mds, xlab = "Dimension 1", 
     ylab = "Dimension 2", main = "MDS plot",
     labels = rownames(mds), cex = 0.8)
pairs(glass) # data don't seem spherical enough
pca <- prcomp(glass)
biplot(pca)
# data seem only partially spherical, another reason not to choose k-means

# assess variability and outliers
boxplot(glass)
meanvar <- cbind(t(t(as.vector(diag(var(glass))))),t(t(as.vector(summary(glass)[4,]))))
colnames(meanvar) <- c("Variance", "Mean")
rownames(meanvar) <- colnames(glass)
meanvar
plot(glass$RI, pch = 16, col = ifelse(glass$RI %in% boxplot.stats(glass$RI)$out, "red", "blue"), 
     main = "Scatter Plot with Outliers Highlighted")
plot(glass$Na, pch = 16, col = ifelse(glass$Na %in% boxplot.stats(glass$Na)$out, "red", "blue"), 
     main = "Scatter Plot with Outliers Highlighted")
plot(glass$Mg, pch = 16, col = ifelse(glass$Mg %in% boxplot.stats(glass$Mg)$out, "red", "blue"), 
     main = "Scatter Plot with Outliers Highlighted")
plot(glass$Al, pch = 16, col = ifelse(glass$Al %in% boxplot.stats(glass$Al)$out, "red", "blue"), 
     main = "Scatter Plot with Outliers Highlighted")
plot(glass$Si, pch = 16, col = ifelse(glass$Si %in% boxplot.stats(glass$Si)$out, "red", "blue"), 
     main = "Scatter Plot with Outliers Highlighted")
plot(glass$K, pch = 16, col = ifelse(glass$K %in% boxplot.stats(glass$K)$out, "red", "blue"), 
     main = "Scatter Plot with Outliers Highlighted")
plot(glass$Ca, pch = 16, col = ifelse(glass$Ca %in% boxplot.stats(glass$Ca)$out, "red", "blue"), 
     main = "Scatter Plot with Outliers Highlighted")
plot(glass$Ba, pch = 16, col = ifelse(glass$Ba %in% boxplot.stats(glass$Ba)$out, "red", "blue"), 
     main = "Scatter Plot with Outliers Highlighted")
plot(glass$Fe, pch = 16, col = ifelse(glass$Fe %in% boxplot.stats(glass$Fe)$out, "red", "blue"), 
     main = "Scatter Plot with Outliers Highlighted")
# all enough variable / with a relevant number of outliers --> k-means could not be a good method,
# while PAM could.

# skeweness of data
hist(glass$RI)
lines(density(glass$RI), col = "blue", lwd = 2)
hist(glass$Na)
lines(density(glass$Na), col = "blue", lwd = 2)
hist(glass$Mg)
lines(density(glass$Mg), col = "blue", lwd = 2)
hist(glass$Al)
lines(density(glass$Al), col = "blue", lwd = 2)
hist(glass$Si)
lines(density(glass$Si), col = "blue", lwd = 2)
hist(glass$K)
lines(density(glass$K), col = "blue", lwd = 2)
hist(glass$Ca)
lines(density(glass$Ca), col = "blue", lwd = 2)
hist(glass$Ba)
lines(density(glass$Ba), col = "blue", lwd = 2)
hist(glass$Fe)
lines(density(glass$Fe), col = "blue", lwd = 2)
# If the histogram and the kernel density plot exhibit similar patterns, the notion of the density being 
# relatively equal. So, we have not eqaul densities, and k-means coul be not performing.
# variables also don't seem all gaussian distributed and symmetric, this worn us that GMM could not perform well.
# One could prefer other functions, more specific for t distributions and skewed distributions.

# why to choose PAM and Hierarchical (could use hierarchical to initialize PAM)
# we don't know the original clusterization and data of this type could assume a nested configuration, for these
# reasons we could use hierarchical clustering
# also, PAM could give us good results due to its ability to handle skewed, nonlinear and non spherical data.


# Gaussian mixture
# GMM assumes that the underlying distribution of each cluster is Gaussian. If data distribution are 
# significantly different, GMM might not perform well.
library(mclust)
m <- Mclust(glass,G=1:15)
summary(m)
summary(m$BIC)
plot(m)
m4 <- Mclust(glass,G=4)
pcaglass <- princomp(glass)
plot(pcaglass$scores,col=m4$classification,
     pch=m4$classification)
plot(mds, col=m4$classification, pch=m4$classification, main = "MDS by GMM clustering")

# hierarchical clustering
hc <- hclust(dist(glass, "manhattan"), method = "average")
# Average linkage might be more suitable as it is more robust to noise and outliers.
plot(hc, main = "Dendrogram for hierarchical clustering - Average")
abline(h = 8.5, col = "red", lty = 2) # 8
# elbow method
wss <- numeric(20)
for (i in 2:21) {
  hc_cut <- cutree(hc, k = i)
  wss[i - 1] <- sum((dist(glass, "manhattan"))^2 * (as.numeric(factor(hc_cut)) == 1))
}
plot(2:21, wss, type = "b", xlab = "Number of clusters", ylab = "Within groups sum of squares",
     main = "Elbow method")
abline(v = c(7,8,9,10), col = "grey", lty = 2) # 8

# hierarchical clustering
clus <- cutree(hc, 8)
plot(mds, col = clus, pch = clus, main = "MDS by Hierarchical clustering")

# PAM
# ASW to find the right K
pasw <- NA
pclusk <- list()
psil <- list()
# Look at K between 2 and 30:
for (k in 2:30){
  # PAM clustering:
  pclusk[[k]] <- pam(dist(glass,"manhattan"),k)
  # Computation of silhouettes:
  psil[[k]] <- silhouette(pclusk[[k]],dist=dist(glass,"manhattan"))
  # ASW needs to be extracted:
  pasw[k] <- summary(psil[[k]])$avg.width
}
# Plot the ASW-values against K:
plot(1:30,pasw,type="l",xlab="Number of clusters",ylab="ASW", main = "ASW")
abline(v = c(1:10), lty = 2, col = "grey")
# 3 seems to be the best cliustering, anyway it could not make makes sense in this context.
# i decided to take into account also the second best value, that is for k = 9
pam3 <- pam(dist(glass,"manhattan"),3)
plot(mds, col = pam3$clustering, pch = pam3$clustering, main = "clustering using PAM")
pam9 <- pam(dist(glass,"manhattan"),9)
plot(mds, col = pam9$clustering, pch = pam9$clustering, main = "clustering using PAM")
# in 2 dimentions, 3 clusters appear clearer, not forcely the best option.

# K-means
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

kmeansres <- gapnc(glass, K.max = 20, B = 200)
kmeansres$nc # 13

vec <- vector()
for (i in 1:100) {
  gap <- gapnc(glass, spaceH0 = "scaledPCA", SE.factor = 2, K.max = 20)
  vec[i] <- gap$nc
}
nc <- ceiling(mean(vec)) # 18
kmeanstai <- kmeans(glass, 18, nstart = 200)

plot(glass, col = kmeanstai$cluster, pch = kmeanstai$cluster, 
     main = "Clustering using K-means")
plot(mds, col = kmeanstai$cluster, pch = kmeanstai$cluster, 
     main = "Clustering using K-means")


# Compare the clusterings (how meaningful and useful are them)

PAM_HC <- adjustedRandIndex(pam9$clustering, clus)
PAM_GMM <- adjustedRandIndex(pam9$clustering, m4$classification)
GMM_HC <- adjustedRandIndex(m4$classification, clus)
GMM_Kmeans <- adjustedRandIndex(m4$classification, kmeanstai$cluster) 
Kmeans_HC <- adjustedRandIndex(kmeanstai$cluster, clus) 
Kmeans_PAM <- adjustedRandIndex(kmeanstai$cluster, pam9$clustering)
cbind(PAM_HC, PAM_GMM, Kmeans_PAM, Kmeans_HC, GMM_HC, GMM_Kmeans)

# Select the best one 
# I would choose PAM clustering, basing my choice on the mds visualizations by clusters, it seems to be the most
# clear.

# Interpretation of the clusters (Comment on other aspects of the data set that could be relevant)

# 3d visualizations
library(rgl)
rgl.open()
plot3d(mds, col = kmeanstai$cluster, size = 6,
       xlab = "Dimension 1", ylab = "Dimension 2", zlab = "Dimension 3",
       main = "Clustering using K-means")
legend3d("topright", legend = levels(kmeanstai$cluster), col = 1:length(unique(kmeanstai$cluster)))
plot3d(mds, col = pam3$clustering, size = 6,
       xlab = "Dimension 1", ylab = "Dimension 2", zlab = "Dimension 3",
       main = "Clustering using K-means")
legend3d("topright", legend = levels(pam3$clustering), col = 1:length(unique(pam3$clustering)))
plot3d(mds, col = pam9$clustering, size = 6,
       xlab = "Dimension 1", ylab = "Dimension 2", zlab = "Dimension 3",
       main = "Clustering using K-means")
legend3d("topright", legend = levels(clus), col = 1:length(unique(clus)))
plot3d(mds, col = clus, size = 6,
       xlab = "Dimension 1", ylab = "Dimension 2", zlab = "Dimension 3",
       main = "Clustering using K-means")
legend3d("topright", legend = levels(clus), col = 1:length(unique(clus)))
plot3d(mds, col = m4$classification, size = 6,
       xlab = "Dimension 1", ylab = "Dimension 2", zlab = "Dimension 3",
       main = "Clustering using K-means")
legend3d("topright", legend = levels(m4$classification), col = 1:length(unique(m4$classification)))

# 3 ---------------------------------------------------------------------------------------
# loading
library(fpc)
library(factoextra)
library(MASS)
library(cluster)
library(dbscan)
setwd("C:/Users/User/Desktop/SCUOLA/UNI/2. Magistrale/II ANNO/primo semestre/modern statistics and big data analytics/slide/datasets")
data1 <- iris[,-5]
data2ns <- read.csv("glass.txt",header=TRUE, sep = " ")
data2 <- as.data.frame(scale(data2ns))
library(pdfCluster)
data(oliveoil)
data3ns <- oliveoil[,3:10]
data3 <- as.data.frame(scale(data3ns))
str(data1)
str(data2)
str(data3)

# tuning
eps_values <- seq(0.01, 1, by = 0.01)
minPts_values <- seq(4, 10, by = 1)

best_params1 <- c(eps = eps_values[1], minPts = minPts_values[1])
best_silhouette1 <- -Inf

for (eps_val in eps_values) {
  for (minPts_val in minPts_values) {
    dbscan_result <- dbscan::dbscan(data1, eps = eps_val, minPts = minPts_val)
    if (length(unique(dbscan_result$cluster)) > 1) {
      clusters <- dbscan_result$cluster
      distances <- dist(data1)
      silhouette <- cluster.stats(distances, clusters)$avg.silwidth
      if (silhouette > best_silhouette1) {
        best_silhouette1 <- silhouette
        best_params1 <- c(eps = eps_val, minPts = minPts_val)
      }
    }
  }
}
cat("Best Parameters for data1:", "\n")
cat("Epsilon (eps):", best_params1["eps"], "\n")
cat("MinPts:", best_params1["minPts"], "\n")

best_params2 <- c(eps = eps_values[1], minPts = minPts_values[1])
best_silhouette2 <- -Inf

for (eps_val in eps_values) {
  for (minPts_val in minPts_values) {
    dbscan_result <- dbscan::dbscan(data2ns, eps = eps_val, minPts = minPts_val)
    if (length(unique(dbscan_result$cluster)) > 1) {
      clusters <- dbscan_result$cluster
      distances <- dist(data2ns)
      silhouette <- cluster.stats(distances, clusters)$avg.silwidth
      if (silhouette > best_silhouette2) {
        best_silhouette2 <- silhouette
        best_params2 <- c(eps = eps_val, minPts = minPts_val)
      }
    }
  }
}
cat("Best Parameters for data2:", "\n")
cat("Epsilon (eps):", best_params2["eps"], "\n")
cat("MinPts:", best_params2["minPts"], "\n")

best_silhouette3 <- -Inf
best_params3 <- c(eps = eps_values[1], minPts = minPts_values[1])

for (eps_val in eps_values) {
  for (minPts_val in minPts_values) {
    dbscan_result <- dbscan::dbscan(data3, eps = eps_val, minPts = minPts_val)
    if (length(unique(dbscan_result$cluster)) > 1) {
      clusters <- dbscan_result$cluster
      distances <- dist(data3)
      silhouette <- cluster.stats(distances, clusters)$avg.silwidth
      if (silhouette > best_silhouette3) {
        best_silhouette3 <- silhouette
        best_params3 <- c(eps = eps_val, minPts = minPts_val)
      }
    }
  }
}
cat("Best Parameters for data3:", "\n")
cat("Epsilon (eps):", best_params3["eps"], "\n")
cat("MinPts:", best_params3["minPts"], "\n")

# visualizing
dbscan1iter <- dbscan::dbscan(data1, eps = best_params1[1], minPts = best_params1[2])
dbscan2iter <- dbscan::dbscan(data2, eps = best_params2[1], minPts = best_params2[2])
dbscan3iter <- dbscan::dbscan(data3, eps = best_params3[1], minPts = best_params3[2])
dbscan1manual <- dbscan::dbscan(data1, eps = 0.42, minPts = 5)
dbscan2manual <- dbscan::dbscan(data2, eps = 0.5, 10)
dbscan3manual <- dbscan::dbscan(data3, eps = 0.9, 9)
fviz_cluster(dbscan1iter, data = data1)
fviz_cluster(dbscan1manual, data = data1)
# i prefer manual result
fviz_cluster(dbscan2iter, data = data2)
fviz_cluster(dbscan2manual, data = data2)
# i prefer iterative result
fviz_cluster(dbscan3iter, data = data3)
fviz_cluster(dbscan3manual, data = data3)
# i prefer iterative result


mds1 <- cmdscale(dist(data1), k = 2)
mds2 <- cmdscale(dist(data2ns, "manhattan"), k = 2)
mds3 <- cmdscale(dist(data3ns), k = 2)
plot(mds1$scores, col = dbscan1manual$cluster+1, shape = dbscan1manual$cluster)
plot(mds2$scores, col = dbscan2iter$cluster+1, shape = dbscan2iter$cluster)
plot(mds3, col = dbscan3iter$cluster+1, shape = dbscan3iter$cluster)
hullplot(mds1, dbscan1manual)
hullplot(mds2, dbscan2iter)
hullplot(mds3, dbscan3iter)

# it seems that only for glass dataset we obtain poor results
