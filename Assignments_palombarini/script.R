# 1 ---------------------------------------------------------------

# DATA LOADING AND CLEANING
oliveoil <- read.table("Oliveoil dataset.txt", header = T)
data <- oliveoil[,3:10]
str(data)
set.seed(665544)

# # PCA
# prolive <- princomp(data)
# summary(prolive)
# plot(prolive,main="Olive oil principal components")

# SCALE
library(fpc)
solive <- scale(data)
# sprolive <- princomp(solive) 
# summary(sprolive)
# plot(sprolive,main="Olive oil principal components")

# K-MEANS
c1k3 <- kmeans(data,centers=3,nstart=100)
plot(data,col=c1k3$cluster,pch=clusym[c1k3$cluster])
points(c1k3$centers,pch="M",cex=2,col=4)
table(c1k3$cluster,oliveoil$macro.area)

sc1k3 <- kmeans(solive,centers=3,nstart=100)
plot(solive,col=sc1k3$cluster,pch=clusym[sc1k3$cluster])
points(sc1k3$centers,pch="M",cex=2,col=4)
table(sc1k3$cluster,oliveoil$macro.area)

c1k3$tot.withinss
sc1k3$tot.withinss

c1k9 <- kmeans(data,centers=9,nstart=100)
plot(data,col=c1k9$cluster,pch=clusym[c1k9$cluster])
points(c1k9$centers,pch="M",cex=2,col=4)
table(c1k9$cluster,oliveoil$region)

sc1k9 <- kmeans(solive,centers=9,nstart=100)
plot(solive,col=sc1k9$cluster,pch=clusym[sc1k9$cluster])
points(sc1k9$centers,pch="M",cex=2,col=4)
table(sc1k9$cluster,oliveoil$region)

c1k9$tot.withinss
sc1k9$tot.withinss


# 2 ---------------------------------------------------------------


# 3 ---------------------------------------------------------------
# data insert, cleaning and visualization
data <- read.table("Boston.txt", header = T)
str(data)
datan <- data[, -c(4,9,10)]
str(datan)
plot(datan)
head(data)

# clustering with K-means
apply(datan, 2, sd)
datas <- data.frame(scale(datan))
pca <- prcomp(datas)
pc1 <- pca$x[, 1]
pc2 <- pca$x[, 2]
plot(pc1, pc2, pch = 20, col = 'blue', main = "PCA: First Two Principal Components")
datapca <- as.data.frame(pca$x[, 1:2])
plot(pca)
percpca <- pca$sdev^2 / sum(pca$sdev^2) *100
sum(percpca[c(1,2)])

# plot(1:length(percpca), percpca, type = "b", main = "Scree Plot")

c1k3 <- kmeans(datapca,centers=3,nstart=200)
table(c1k3$cluster)

plot(datapca$PC1, datapca$PC2, pch = 20, col = c1k3$cluster, 
     main = "PCA with K-means Clustering", 
     xlab = "PC1", ylab = "PC2")
points(c1k3$centers,pch="X",cex=2,col=c(3,1,2))

c1k4 <- kmeans(datapca,centers=4,nstart=200)
table(c1k4$cluster)

plot(datapca$PC1, datapca$PC2, pch = 20, col = c1k4$cluster, 
     main = "PCA with K-means Clustering", 
     xlab = "PC1", ylab = "PC2",)
points(c1k4$centers,pch="X",cex=2,col=c(4,3,2,1))


c1k5 <- kmeans(datapca,centers=5,nstart=200)
table(c1k5$cluster)

plot(datapca$PC1, datapca$PC2, pch = 20, col = c1k5$cluster, 
     main = "PCA with K-means Clustering", 
     xlab = "PC1", ylab = "PC2",)
points(c1k5$centers,pch="X",cex=2,col=c(5,4,3,2,1))

plot(datan$age, datan$black, pch = 20, col = c1k5$cluster, 
     main = "PCA with K-means Clustering", 
     xlab = "PC1", ylab = "PC2",)
points(c1k5$centers,pch="X",cex=2,col=c(5,4,3,2,1))


# 4 ---------------------------------------------------------------
data <- read.table("Boston.txt", header = T)
str(data)
datan <- data[, -c(4,9,10)]
datas <- as.matrix(data.frame(scale(datan)))

library(pracma)
kmpp <- function(X, k) {
  n <- nrow(X)
  C <- numeric(k)
  C[1] <- sample(1:n, 1)
  for (i in 2:k) {
    dm <- distmat(X, X[C, ])
    pr <- apply(dm, 1, min); pr[C] <- 0
    C[i] <- sample(1:n, 1, prob = pr)
  }
  kmeans(X, X[C, ])
}

kmpp <- kmpp(datas, 4)

km1 <- kmeans(data,centers=4,nstart=1)
km100 <- kmeans(data,centers=4,nstart=100)

kmpp$tot.withinss
km1$tot.withinss
km100$tot.withinss
