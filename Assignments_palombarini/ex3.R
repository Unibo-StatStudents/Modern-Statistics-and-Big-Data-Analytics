# 1 ------------------------------------------------------------
x1 = c("blue", 1, 1, 0, 12)
x2 = c("red", 0, 0, NA, NA)
x3 = c("red", 1, 0, NA, 17)
x4 = c("green", 1, 0, 0, 21)
library(cluster)
n <- 4
p <- 5
gowerindex <- list()

# manual try 1
data <- rbind(x1,x2,x3,x4)
gowfunc <- function(data) {
  n <- nrow(data)
  p <- ncol(data)
  types <- sapply(data, function(x) ifelse(is.numeric(x), "numeric", "categorical"))
  
  gower_dist <- matrix(0, n, n)
  w <- matrix(0,n,n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      sum_sij <- 0
      non_missing_values <- 0
      for (k in 1:p) {
        if (is.na(data[i, k]) || is.na(data[j, k])) {
          sum_sij <- sum_sij + 0 
        } else if (types[k] == "numeric") {
          sum_sij <- sum_sij + abs(data[i, k] - data[j, k])
          non_missing_values <- non_missing_values + 1
          w[i,j] <- abs(data[i, k] - data[j, k])
        } else {
          if (data[i, k] == data[j, k]) {
            sum_sij <- sum_sij + 0
            w[i,j] <- 1
          } else {
            sum_sij <- sum_sij + 1
            w[i,j] <- 1
          }
          non_missing_values <- non_missing_values + 1
        }
      }
      gower_dist[i, j] <- (sum_sij/w[i,j])/non_missing_values
    }
  }
  
  rownames(gower_dist) <- colnames(gower_dist) <- 1:n
  return(gower_dist)
}

gow1 <- gowfunc(data)
gowerindex$manual1$d12 <- gow1[1,2]
gowerindex$manual1$d13 <- gow1[1,3]
gowerindex$manual1$d14 <- gow1[1,4]
gowerindex$manual1$d23 <- gow1[2,3]
gowerindex$manual1$d24 <- gow1[2,4]
gowerindex$manual1$d34 <- gow1[3,4]

# # manual try 2
# gowerindex$manual2$d12 <- sum(c(1,1,1,0,0))/sum(1,1,1,0,0)
# gowerindex$manual2$d13 <- sum(c(1,0,1,0,(17-12))/(17-12))/sum(1,1,1,0,1)
# gowerindex$manual2$d14 <- sum(c(1,0,1,1,(21-12))/(21-12))/sum(1,1,1,1,1)
# gowerindex$manual2$d23 <- sum(0,1,1,0,0)/sum(1,1,1,0,0)
# gowerindex$manual2$d24 <- sum(1,1,1,0,0)/sum(1,1,1,0,0)
# gowerindex$manual2$d34 <- sum(c(1,0,1,0,(21-17))/(21-17))/sum(1,1,1,0,1)

# daisy
data <- as.data.frame(matrix(c(x1, x2, x3, x4), 4,5, byrow = T))
data$V1 <- as.factor(data$V1)
data$V2 <- as.factor(data$V2)
data$V3 <- as.factor(data$V3)
data$V4 <- as.factor(data$V4)
data$V5 <- as.numeric(data$V5)
str(data)
daisy <- daisy(data, type=list(asymm=c(2,3,4)))
gowerindex$daisy$d12 <- daisy[1]
gowerindex$daisy$d13 <- daisy[2]
gowerindex$daisy$d14 <- daisy[3]
gowerindex$daisy$d23 <- daisy[4]
gowerindex$daisy$d24 <- daisy[5]
gowerindex$daisy$d34 <- daisy[6]

# visualization
plot(as.matrix(gowerindex$manual1), xlab = "Index", ylab = "Value", col = "darkgreen", pch = 20, xaxt = "n", 
     cex = 3.5, main = "Gower indexes for all combinations of vectors 1-4, manually and with daisy() function", ylim =c(0,1))
axis(1, at = 1:6, labels = names(gowerindex$manual1), las = 1)
# points(y = gowerindex$manual2, x = 1:6, col = "tan2", pch = 20, cex = 2.8)
points(y = gowerindex$daisy, x = 1:6, col = "purple1", pch = 4, cex = 3)
legend(x = "topright", legend = c("manual1", "daisy"), 
       col = c("darkgreen", "purple1"),pch = c(20,4), 
       pt.cex = 3, cex = 0.7, x.intersp = 1, xjust = 1, inset = c(-0.15, 0.05), 
       title = "legend", title.adj = 0.2, xpd = TRUE)

# 2 ------------------------------------------------------------------------
# b
# corr
x = c(1, 2, NaN)
y = c(2, NaN, 3)
z = c(NaN, 4, 5)
d_Cxy = 0.5*(1 - cor(x, y))
d_Cxz = 0.5*(1 - cor(x, z))
d_Cyz = 0.5*(1 - cor(y, z))
cbind(d_Cxy, d_Cxz, d_Cyz) # tiangular inequality doesn't stand

# gower 
library(cluster)
x <- c(1, 2, NaN)
y <- c(2, NaN, 3)
z <- c(4, 4, NaN)
data <- data.frame(x, y, z)
gower_dist <- daisy(data, metric = "gower")
print(gower_dist) # tiangular inequality doesn't stand


# 4 ------------------------------------------------------------------------
### loading data ###
data <- read.table("Covid2021.txt")
data1 <- data[,5:559]
datam <- as.matrix(data1)
str(data1)

# understanding data ----------------------------------------------------------------
plot(1:555,data1[1,],type="l",ylim=c(0,25),
     ylab="New cases over one week per 1000 inhabitants",
     xlab="Day (1 April 2020-7 October 2021)")
for(i in 1:2) {
  points(1:555,data1[i,],type="l")
}
summary(data[,6])
sum(data[,-c(1:4)])


sum4days <- function(vec = c(1:555)){
  sum <- vector()
  
  for (i in vec) {
    s <- sum(data1[,i])
    sum[i] <- s
  }
  return(sum)
}
sum_for_days <- sum4days()

summary4days <- function(vec = c(1:555)){
  mean <- list()
  
  for (i in vec) {
    m <- summary(data1[,i])
    mean[[i]] <- m
  }
  return(mean)
}
summary_for_days <- summary4days()

sum4country <- function(vec = c(1:179)){
  sum <- vector()
  
  for (i in vec) {
    s <- sum(data1[i,])
    sum[i] <- s
  }
  return(sum)
}
sum_for_country <- sum4country()


summary4country <- function(vec = c(1:179)){
  mean <- list()
  
  for (i in vec) {
    m <- summary(data1[i,])
    mean[[i]] <- m
  }
  return(mean)
}
summary_for_country <- summary4country()


# clustering ------------------------------------------------------------------------------
library(mclust)
library(sn)
library(cluster)
library(pdfCluster)

# Step 1: dissimilarity indexes
# Eeuclidean/Manhattan/Mahalanobis
distances <- list()
distances$deuc <- dist(datam,method="euclidean")
distances$dman <- dist(datam,method="manhattan")

# Correlation: d_C/d_CN
distances$d_C <- (1-cor(t(datam)))/2
distances$d_CR <- 1-abs(cor(t(datam)))
distances$d_C <- as.dist(distances$d_C)
distances$d_CR <- as.dist(distances$d_CR)

# Step 2: Hierarchical clustering
# Single Linkage (or nearest neighbour) clustering:
hierclust <- list()
hierclust$single$euclidean <- hclust(distances$deuc,method="single")
plot(hierclust$single$euclidean, cex = 0.4, main = "dendrogram of single linkage with euclidean") 
abline(h = 26, lty = 2, col = "red") # 16
hierclust$single$manhattan <- hclust(distances$dman,method="single")
plot(hierclust$single$manhattan, cex = 0.4, main = "dendrogram of single linkage with manhattan")
abline(h = 400, lty = 2, col = "red") # 12
hierclust$single$C <- hclust(distances$d_C,method="single")
plot(hierclust$single$C, cex = 0.4, main = "dendrogram of single linkage with C")
abline(h = 0.205, lty = 2, col = "red")  # 12
hierclust$single$CR <- hclust(distances$d_CR,method="single")
plot(hierclust$single$CR, cex = 0.4, main = "dendrogram of single linkage with CR")
abline(h = 0.405, lty = 2, col = "red") # 12

# Complete Linkage (or furthest neighbour) clustering:
hierclust$complete$euclidean <- hclust(distances$deuc,method="complete")
plot(hierclust$single$euclidean, cex = 0.4, main = "dendrogram of complete linkage with euclidean") 
abline(h = 26, lty = 2, col = "red") # 16
hierclust$complete$manhattan <- hclust(distances$dman,method="complete")
plot(hierclust$single$manhattan, cex = 0.4, main = "dendrogram of complete linkage with manhattan") 
abline(h = 400, lty = 2, col = "red") # 12
hierclust$complete$C <- hclust(distances$d_C,method="complete")
plot(hierclust$complete$C, cex = 0.4, main = "dendrogram of complete linkage with C")
abline(h = 0.5, lty = 2, col = "red") # 15
hierclust$complete$CR <- hclust(distances$d_CR,method="complete")
plot(hierclust$complete$CR, cex = 0.4, main = "dendrogram of complete linkage with CR")
abline(h = 0.9, lty = 2, col = "red") # 16

# Average Linkage (or UPGMA) clustering:
hierclust$avarage$euclidean <- hclust(distances$deuc,method="average")
plot(hierclust$avarage$euclidean, cex = 0.4, main = "dendrogram of average linkage with euclidean")
abline(h = 40, lty = 2, col = "red") # 10
hierclust$avarage$manhattan <- hclust(distances$dman,method="average")
plot(hierclust$avarage$manhattan, cex = 0.4, main = "dendrogram of average linkage with manhattan")
abline(h = 585, lty = 2, col = "red") # 11
hierclust$avarage$C <- hclust(distances$d_C,method="average")
plot(hierclust$avarage$C, cex = 0.4, main = "dendrogram of single average with C")
abline(h = 0.38, lty = 2, col = "red") # 10
hierclust$avarage$CR <- hclust(distances$d_CR,method="average")
plot(hierclust$avarage$CR, cex = 0.4, main = "dendrogram of single average with CR")
abline(h = 0.75, lty = 2, col = "red") # 10


# Considering Average linkage we have pretty much all similar results, i.e. k = 10

# choosing number of clusters, visualization and interpretation
library(smacof)
mds <- list()
mds$euclidean<- mds(distances$deuc)
mds$manhattan <- mds(distances$dman)
mds$CR <- mds(distances$d_C)
mds$C <- mds(distances$d_CR)

clust <- list()
clust$average$euclidean <- cutree(hierclust$avarage$euclidean, 10)
plot(mds$euclidean$conf, col = clust$average$euclidean, 
     main = "MDS visualization by hierachical clustering with average link method and euclidean distance")
clust$average$manhattan <- cutree(hierclust$avarage$manhattan, 11)
plot(mds$manhattan$conf, col = clust$average$manhattan,
     main = "MDS visualization by hierachical clustering with average link method and manhattan distance")
clust$average$C <- cutree(hierclust$avarage$C, 10)
plot(mds$C$conf, col = clust$average$C,
     main = "MDS visualization by hierachical clustering with average link method and d_C distance")
clust$average$CR <- cutree(hierclust$avarage$CR, 10)
plot(mds$CR$conf, col = clust$average$CR,
     main = "MDS visualization by hierachical clustering with average link method and D_CR distance")

# other plots
library(ggplot2)
library(dplyr)
mds_df <- data.frame(x = mds$euclidean$conf[, 1], y = mds$euclidean$conf[, 2])
mds_df$Country <- row.names(mds$euclidean$conf)
clust_df <- data.frame(Country = names(clust$average$euclidean), Cluster = clust$average$euclidean)
mds_df <- left_join(mds_df, clust_df, by = "Country")
ggplot(mds_df, aes(x = x, y = y)) +
  geom_point(aes(color = as.factor(Cluster))) +
  geom_text(aes(label = as.factor(Cluster)), size = 3, vjust = 1.5) +
  scale_color_discrete(name = "Cluster") +
  labs(title = "MDS visualization by hierarchical clustering with average link method and Euclidean distance")

mds_df_manhattan <- data.frame(x = mds$manhattan$conf[, 1], y = mds$manhattan$conf[, 2])
mds_df_manhattan$Country <- row.names(mds$manhattan$conf)
clust_df_manhattan <- data.frame(Country = names(clust$average$manhattan), Cluster = clust$average$manhattan)
mds_df_manhattan <- left_join(mds_df_manhattan, clust_df_manhattan, by = "Country")
ggplot(mds_df_manhattan, aes(x = x, y = y)) +
  geom_point(aes(color = as.factor(Cluster))) +
  geom_text(aes(label = as.factor(Cluster)), size = 3, vjust = 1.5) +
  scale_color_discrete(name = "Cluster") +
  labs(title = "MDS visualization by hierarchical clustering with average link method and Manhattan distance")

mds_df_C <- data.frame(x = mds$C$conf[, 1], y = mds$C$conf[, 2])
mds_df_C$Country <- row.names(mds$C$conf)
clust_df_C <- data.frame(Country = names(clust$average$C), Cluster = clust$average$C)
mds_df_C <- left_join(mds_df_C, clust_df_C, by = "Country")
ggplot(mds_df_C, aes(x = x, y = y)) +
  geom_point(aes(color = as.factor(Cluster))) +
  geom_text(aes(label = as.factor(Cluster)), size = 3, vjust = 1.5) +
  scale_color_discrete(name = "Cluster") +
  labs(title = "MDS visualization by hierarchical clustering with average link method and d_C distance")

mds_df_CR <- data.frame(x = mds$CR$conf[, 1], y = mds$CR$conf[, 2])
mds_df_CR$Country <- row.names(mds$CR$conf)
clust_df_CR <- data.frame(Country = names(clust$average$CR), Cluster = clust$average$CR)
mds_df_CR <- left_join(mds_df_CR, clust_df_CR, by = "Country")
ggplot(mds_df_CR, aes(x = x, y = y)) +
  geom_point(aes(color = as.factor(Cluster))) +
  geom_text(aes(label = as.factor(Cluster)), size = 3, vjust = 1.5) +
  scale_color_discrete(name = "Cluster") +
  labs(title = "MDS visualization by hierarchical clustering with average link method and D_CR distance")

# Both before and after the visualization i'd say that correlation indexed give us a better prospective

# comparing these plots with the division by continents of data we notice a not so good estimate of the clusters
# (if we consider the continents as thetrue clusters)
plot(mds$euclidean$conf, col = as.factor(data$continent), 
     main = "MDS visualization by hierachical clustering with average link method and euclidean distance")
plot(mds$manhattan$conf, col = as.factor(data$continent),
     main = "MDS visualization by hierachical clustering with average link method and manhattan distance")
plot(mds$C$conf, col = as.factor(data$continent),
     main = "MDS visualization by hierachical clustering with average link method and d_C distance")
plot(mds$CR$conf, col = as.factor(data$continent),
     main = "MDS visualization by hierachical clustering with average link method and D_CR distance")

# Step 3: similarities to evaluate clusterings
# evaluating if they are a lot different and define, if possible, the best one
# Adjusted Rand index
adjustedRandIndex(clust$average$euclidean,clust$average$manhattan)
adjustedRandIndex(clust$average$euclidean,clust$average$CR)
adjustedRandIndex(clust$average$euclidean,clust$average$C)
adjustedRandIndex(clust$average$manhattan,clust$average$CR)
adjustedRandIndex(clust$average$manhattan,clust$average$C)
adjustedRandIndex(clust$average$CR,clust$average$C)
# all pretty different, exeption made for the two correlation indexes

# ---------------------------------- do not consider --------------------------------------------------------
# # K-means gap statistic to find best k and comparing with k found with hclust
# gapnc <- function(data,FUNcluster=kmeans,
#                   K.max=10, B = 100, d.power = 2,
#                   spaceH0 ="scaledPCA",
#                   method ="globalSEmax", SE.factor = 2,...){ # the three dots allows us to specify other parameters that need a
#   # function specification
#   # As in original clusGap function the ... arguments are passed on
#   # to the clustering method FUNcluster (kmeans).
#   # Run clusGap
#   gap1 <- clusGap(data,kmeans,K.max, B, d.power,spaceH0,...)
#   # Find optimal number of clusters; note that the method for
#   # finding the optimum and the SE.factor q need to be specified here.
#   nc <- maxSE(gap1$Tab[,3],gap1$Tab[,4],method, SE.factor)
#   # Re-run kmeans with optimal nc.
#   kmopt <- kmeans(data,nc,...)
#   out <- list()
#   out$gapout <- gap1
#   out$nc <- nc
#   out$kmopt <- kmopt
#   out
# }
# 
# 
# # find another way to reppresent data in order to apply pca and to fin another distance index
# # let's calculate the average monthly increase
# data2 <- data1[,-c(541:555)] # let's remove 7 days to have exactly 18 month of 30 days
# length(data2)
# 540/30
# 
# data3 <- data.frame(matrix(nrow = nrow(data2), ncol = 18))
# for (i in 1:18) {
#   start <- (i - 1) * 30 + 1
#   end <- i * 30
#   col <- paste0("Month", i)
#   data3[, i] <- rowMeans(data2[, start:end])
# }
# 
# str(data3)
# plot(data3)
# pcdata <- prcomp(data3, scale = F)
# summ <- summary(pcdata)
# screeplot(pcdata, type = "lines")
# 
# gap <- gapnc(datam, K.max = 20, B = 100, SE.factor = 1)
# plot(datam, col = gap$kmopt$cluster)
# 
# # run K.means with the best K between the two
# -----------------------------------------------------------------------------------------------------------------

# Heatmaps
heatmap(as.matrix(data1), Rowv = TRUE, Colv = TRUE, col = grey(seq(1, 0, -0.01)))
v <- dist(t(data1),method="binary")
varclust <- hclust(v,method="average")
heatmap(as.matrix(data1), Rowv = as.dendrogram(hierclust$single$euclidean), 
        Colv = as.dendrogram(varclust), col = grey(seq(1, 0, -0.01)))

# here i tried a different method more suitable for time series ---------------------------------------------
library(dtw) # using dinamic tyme warping
covid2021cl <- read.table("Covid2021.txt")
covid2021cl <- covid2021cl[5:559]
# DTW distance matrix
dtw_matrix <- matrix(0, nrow = nrow(covid2021cl), ncol = nrow(covid2021cl))
for (i in 1:(nrow(covid2021cl) - 1)) {
  for (j in (i + 1):nrow(covid2021cl)) {
    dtw_result <- dtw(covid2021cl[i,], covid2021cl[j,], keep = TRUE)
    dtw_matrix[i, j] <- dtw_result$distance
    dtw_matrix[j, i] <- dtw_result$distance
  }
}

# hierarchical clustering
hc_dtw <- hclust(as.dist(dtw_matrix), method = "complete")
# dendrogram
plot(hc_dtw, main = "Hierarchical Clustering with Dynamic Time Warping", xlab = "Countries", sub = "", ylab = "DTW Distance")
abline(h = 47, lty = 2, col = "red")
clusters <- cutree(hc_dtw, h = 14)
cluster_colors <- rainbow(length(unique(clusters)))
data_colors <- cluster_colors[clusters]
plot(1:555, covid2021cl[1,], type = "l", ylim = c(0, 25), ylab = "New cases over one week per 1000 inhabitants", xlab = "Day (1 April 2020-7 October 2021)", col = data_colors[1], lwd = 2)

for (i in 2:179) {
  lines(1:555, covid2021cl[i,], type = "l", col = data_colors[i], lwd = 2)
}


mds <- cmdscale(dtw_matrix, k = 2)
plot(mds, type = "n", xlab = "MDS Dimension 1", ylab = "MDS Dimension 2", main = "MDS Plot of Time Series Data")
points(mds, labels = rownames(dtw_matrix), col = clusters)



 
