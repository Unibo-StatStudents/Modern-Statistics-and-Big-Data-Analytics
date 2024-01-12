# 1 -----------------------------------------------------------
# a)
data1 <- read.table("Oliveoil dataset.txt", header = T)
data1 <- data1[,-c(1,2)]
data2 <- read.table("Clusterdata2.txt", header = T)
colnames(data2) <- c("v1","v2")
cor(data1)
cor(data2)
library(cluster)
gapnc <- function(data,FUNcluster=kmeans,
                  K.max=10, B = 100, d.power = 2,
                  spaceH0 ="scaledPCA",
                  method ="globalSEmax", SE.factor = 2,...){ # the three dots allows us to specify other parameters that need a
  # function specification
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

cgnc1 <- gapnc(data1,nstart=100, K.max = 20)
# plot(cgnc1$gapout)
# print(cgnc1$gapout,method="globalSEmax",SE.factor=2)
cgnc1$nc
# plot(data1,col=cgnc1$kmopt$cluster)

cgnc2 <- gapnc(data2,nstart=100, K.max = 20)
# plot(cgnc2$gapout)
# print(cgnc2$gapout,method="globalSEmax",SE.factor=2)
cgnc2$nc
# plot(data2,col=cgnc2$kmopt$cluster)

# b)
library(sn)
str(data2)
x1 <- runif(150, min = min(data2$v1), max = max(data2$v1))
x2 <- runif(150, min = min(data2$v2), max = max(data2$v2))
dataunif <- data.frame(x1, x2)

km <- function(data = runif(150), nk = 10, nstart = 100, scatter = c(1, 2, 3)) {
  result <- list()
  w <- list()
  
  for (i in 1:nk) {
    kmeans <- kmeans(x = data, centers = i, nstart = nstart)
    name <- paste("kmeans", i, sep = "")
    result[[name]] <- kmeans$cluster
    w[i] <- log(kmeans$tot.withinss)
  }
  plot(data, col = result[[scatter[1]]], main = paste("Scatterplot for", scatter[1], "clusters", sep = " "))
  plot(data, col = result[[scatter[2]]], main = paste("Scatterplot for", scatter[2], "clusters", sep = " "))
  plot(data, col = result[[scatter[3]]], main = paste("Scatterplot for", scatter[3], "clusters", sep = " "))
  return(list(result = result, logw = w))
  }
result <- km(dataunif, nk = 10, nstart = 100, scatter = c(3,6,10))

plot(1:10,cgnc1$gapout$Tab[1:10,1],xlab="k",ylab="log S_k",type="l", ylim = c(0,20))
points(1:10,result$logw,xlab="k",ylab="log S_k",type="l",lty=2)
legend(4,5,c("log S_k data2","log S_k dataunif"),lty=1:2)

# 2 -------------------------------------------------------------------------------
library(sn)
datasets <- function(n = 100, nobs = 150) {
  lista <- list()
  
  for (i in 1:n) {
    set.seed(i+86428)
    v1 <- c(rnorm(50, 0, 1), rsn(70, 5, 1, 8), rnorm(30, 6, 1))
    v2 <- c(rnorm(50, 0, 1), rsn(70, 0, 1, 8), 8 + rt(30, 5))
    dataset <- cbind(v1, v2)
    name <- paste("data_", i, sep = "")
    lista[[name]] <- as.data.frame(dataset)
  }
  return(lista)
}

datasetslist <- datasets(n = 100)

library(cluster)
gapnc <- function(data,FUNcluster=kmeans,
                  K.max=10, B = 100, d.power = 2,
                  spaceH0 ="scaledPCA",
                  method ="globalSEmax", SE.factor = 2,...){ # the three dots allows us to specify other parameters that need a
  # function specification
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

plots <- function(data) {
  vec1 <- vector()
  vec2 <- vector()
  vec3 <- vector()
  vec4 <- vector()
  result <- data.frame()
  
  for (i in 1:length(data)) {
    gap <- gapnc(data[[i]], spaceH0 = "original", SE.factor = 1)
    vec1[[i]] <- gap$nc
  }
  for (i in 1:length(data)) {
    gap <- gapnc(data[[i]], spaceH0 = "original", SE.factor = 2)
    vec2[[i]] <- gap$nc
  }
  for (i in 1:length(data)) {
    gap <- gapnc(data[[i]], spaceH0 = "scaledPCA", SE.factor = 1)
    vec3[[i]] <- gap$nc
  }
  for (i in 1:length(data)) {
    gap <- gapnc(data[[i]], spaceH0 = "scaledPCA", SE.factor = 2)
    vec4[[i]] <- gap$nc
  }
  result <- data.frame(vec1,vec2,vec3,vec4)
}

result <- plots(datasetslist)

# library(GGally)
# ggpairs(datagroup[,-3], aes(col = as.factor(datagroup)))
# ggpairs(b222[,-5], aes(col = as.factor(b222$cluster)))

library(ggplot2)
library(reshape2)
result2 <- data.frame(rep(1:100,4),melt(result))
colnames(result2)[1] <- "n"
ggplot(result2, aes(x = n, y = value, col = variable)) +
  geom_line() +
  labs(x = "Observation", y = "Value") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal()
plot1 <- ggplot(data.frame(result$vec1), aes(x = 1:100, y = result.vec1)) +
  geom_line(color = "red") +
  labs(x = "Observation", y = "Value") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal()
plot2 <- ggplot(data.frame(result$vec2), aes(x = 1:100, y = result.vec2)) +
  geom_line(color = "blue") +
  labs(x = "Observation", y = "Value") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal()
plot3 <- ggplot(data.frame(result$vec3), aes(x = 1:100, y = result.vec3)) +
  geom_line(color = "green") +
  labs(x = "Observation", y = "Value") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal()
plot4 <- ggplot(data.frame(result$vec4), aes(x = 1:100, y = result.vec4)) +
  geom_line(color = "yellow") +
  labs(x = "Observation", y = "Value") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal()

library(gridExtra)
grid.arrange(plot1, plot2, plot3, plot4, ncol = 2)

# boxplots
ggplot(result2, aes(x = n, y = value, col = variable)) +
  geom_boxplot()
result24 <- subset(result2, result2$value < 4)
ggplot(result24, aes(x = n, y = value, col = variable)) +
  geom_boxplot()
# ve3 seems the best one (less variable)

# evaluate the less variable, that could be the best option
diag(var(result))
mean(result$vec1)
mean(result$vec2)
mean(result$vec3)
mean(result$vec4)

# vec3 is the best, let's visualize an example
gap <- gapnc(datasetslist$data_3, spaceH0 = "scaledPCA", SE.factor = 1)
a <- data.frame(datasetslist$data_3, gap$kmopt$cluster)
colnames(a)[3] <- "cluster"
centers <- as.data.frame(cbind(gap$kmopt$centers,c(NA,NA,NA)))
colnames(centers)[3] <- "cluster"
a2 <- rbind(a, centers)
plotfinal <- ggplot(a2) + aes(x = v1, y = v2) +
  geom_point(aes(shape = ifelse(is.na(a2$cluster), "Centers", "Points")),
             size = ifelse(is.na(a2$cluster), 7, 5),
             col = ifelse(is.na(a2$cluster), "purple", a2$cluster),
             stroke = ifelse(is.na(a2$cluster), 2, 1)) +
  scale_shape_manual(values = c("Centers" = 4, "Points" = 20)) +
  labs(title = "Clusterized observations") +
  guides(shape = guide_legend(title = "Legend")) + 
  theme_minimal()
plotiniz <- ggplot(a) + 
  aes(x = v1, y = v2) + 
  geom_point(col = "grey", size = 5, shape = 20) +
  theme_minimal() +
  labs(title = "Non-clusterized observations")

grid.arrange(plotiniz, plotfinal, ncol = 2)


# let's do the same with a different dataset model
library(mvtnorm)
library(ggplot2)
library(fpc)
library(cluster)
n <- c(50, 100, 150, 200)
# set.seed(2231344)
# means <- matrix(rbind(sample(20:30, 4), 
#                       sample(40:60, 4), 
#                       sample(65:85, 4), 
#                       sample(90:110, 4)), ncol = 4, byrow = TRUE)
# # sigma <- list(matrix(sample(seq(0.01, 0.3, by = 0.01), 16), ncol = 4, byrow = T), 
# #               matrix(sample(seq(0.01, 0.3, by = 0.01), 16), ncol = 4, byrow = T), 
# #               matrix(sample(seq(0.01, 0.3, by = 0.01), 16), ncol = 4, byrow = T), 
# #               matrix(sample(seq(0.01, 0.3, by = 0.01), 16), ncol = 4, byrow = T))
# # sigma[[1]][,1] <- sigma[[1]][1,]
# # sigma[[1]][-c(1,2),2] <- sigma[[1]][2,-c(1,2)]
# # sigma[[1]][4,3] <- sigma[[1]][3,4]
# norm1 <- rmvnorm(n = n[1], mean = means[1, ])#, sigma = sigma[[1]])
# norm2 <- rmvnorm(n = n[2], mean = means[2, ])#, sigma = sigma[[2]])
# norm3 <- rmvnorm(n = n[3], mean = means[3, ])#, sigma = sigma[[3]])
# norm4 <- rmvnorm(n = n[4], mean = means[4, ])#, sigma = sigma[[4]])
# datalist <- list(norm1, norm2, norm3, norm4)
# datanorm <- do.call(rbind, datalist)
# datanormdf <- as.data.frame(datanorm)
# 
# # visualization of k-means application on one example for this type of dataset
# # first two dimentions
# resultnorm1 <- gapnc(datanormdf[,c(1,2)])
# plot(datanormdf[,c(1,2)], pch = 19, col = rep(1:4, n), main = "Datasets")
# plot(datanormdf[,c(1,2)], col = resultnorm1$kmopt$cluster, main = "Data after K-means clustering")
# table(resultnorm1$kmopt$cluster, rep(1:4, n))
# # dimentions 3 and 4
# resultnorm2 <- gapnc(datanormdf[,c(3,4)])
# plot(datanormdf[,c(3,4)], pch = 19, col = rep(1:4, n), main = "Datasets")
# plot(datanormdf[,c(3,4)], col = resultnorm2$kmopt$cluster, main = "Data after K-means clustering")
# table(resultnorm2$kmopt$cluster, rep(1:4, n))
# # all dimentions
# resultnorm <- gapnc(datanormdf)
# plot(datanormdf, col = rep(1:4, n), main = "Datasets")
# plot(datanormdf, col = resultnorm$kmopt$cluster, main = "Data after K-means clustering")

# generating 100 dataset of this kind
datasets2 <- function(n = 100, n.obs = c(50,100,150,200)) {
  lista <- list()
  
  for (i in 1:n) {
    set.seed(i+231344)
    means <- matrix(rbind(sample(20:30, 4), 
                          sample(40:60, 4), 
                          sample(65:85, 4), 
                          sample(90:110, 4)), ncol = 4, byrow = TRUE)
    norm1 <- rmvnorm(n = n.obs[1], mean = means[1, ])
    norm2 <- rmvnorm(n = n.obs[2], mean = means[2, ])
    norm3 <- rmvnorm(n = n.obs[3], mean = means[3, ])
    norm4 <- rmvnorm(n = n.obs[4], mean = means[4, ])
    datalist <- list(norm1, norm2, norm3, norm4)
    datanorm <- do.call(rbind, datalist)
    datanormdf <- as.data.frame(datanorm)
    name <- paste("data_norm", i, sep = "")
    lista[[name]] <- as.data.frame(datanormdf)
  }
  return(lista)
}
datasetslist2 <- datasets2(n = 100)
# dataad <- data.frame(datasetslist2$data_norm1, rep(1:4,n))
# plot(dataad[-5], col = dataad$rep.1.4..n.)
resultsnorm <- plots(datasetslist2)

# let's visualize the results
library(ggplot2)
library(reshape2)
resultsnorm2 <- data.frame(rep(1:100,4),melt(resultsnorm))
colnames(resultsnorm2)[1] <- "n"
ggplot(resultsnorm2, aes(x = n, y = value, col = variable)) +
  geom_line() +
  labs(x = "Observation", y = "Value") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal()
plot1.2 <- ggplot(data.frame(resultsnorm$vec1), aes(x = 1:100, y = resultsnorm.vec1)) +
  geom_line(color = "red") +
  labs(x = "Observation", y = "Value") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal()
plot2.2 <- ggplot(data.frame(resultsnorm$vec2), aes(x = 1:100, y = resultsnorm.vec2)) +
  geom_line(color = "blue") +
  labs(x = "Observation", y = "Value") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal()
plot3.2 <- ggplot(data.frame(resultsnorm$vec3), aes(x = 1:100, y = resultsnorm.vec3)) +
  geom_line(color = "green") +
  labs(x = "Observation", y = "Value") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal()
plot4.2 <- ggplot(data.frame(resultsnorm$vec4), aes(x = 1:100, y = resultsnorm.vec4)) +
  geom_line(color = "yellow") +
  labs(x = "Observation", y = "Value") +
  scale_color_brewer(palette = "Set1") +
  theme_minimal()

library(gridExtra)
grid.arrange(plot1.2, plot2.2, plot3.2, plot4.2, ncol = 2)

# evaluate the less variable, that could be the best option
diag(var(resultsnorm))
mean(resultsnorm$vec1)
mean(resultsnorm$vec2)
mean(resultsnorm$vec3)
mean(resultsnorm$vec4)

# vec1 is the best, let's visualize an example

# # visualization on 2 dimentions
# # we need a PCA to reduce dimentions
# gap2 <- gapnc(datasetslist2$data_norm1[,c(1,2)], spaceH0 = "scaledPCA", SE.factor = 1)
# b <- data.frame(datasetslist2$data_norm1[,c(1,2)], gap2$kmopt$cluster)
# colnames(b)[3] <- "cluster"
# centers2 <- as.data.frame(cbind(gap2$kmopt$centers,c(NA,NA,NA)))
# colnames(centers2)[3] <- "cluster"
# b2 <- rbind(b, centers2)
# plotfinal2 <- ggplot(b2) + aes(x = V1, y = V2) +
#   geom_point(aes(shape = ifelse(is.na(b2$cluster), "Centers", "Points")),
#              size = ifelse(is.na(b2$cluster), 7, 5),
#              col = ifelse(is.na(b2$cluster), "purple", b2$cluster),
#              stroke = ifelse(is.na(b2$cluster), 2, 1)) +
#   scale_shape_manual(values = c("Centers" = 4, "Points" = 20)) +
#   labs(title = "Clusterized observations") +
#   guides(shape = guide_legend(title = "Legend")) + 
#   theme_minimal()
# datahalf <- data.frame(datasetslist2$data_norm1[,c(1,2)], rep(1:4, n))
# colnames(datahalf)[3] <- "n"
# plotiniz2 <- ggplot(datahalf) + 
#   aes(x = V1, y = V2, col = n) + 
#   geom_point(size = 5, shape = 20) +
#   theme_minimal() +
#   labs(title = "Non-clusterized observations") + 
#   theme(legend.position = "none")
# 
# grid.arrange(plotiniz2, plotfinal2, ncol = 2)
# 
# # confusion matrix
# tab <- table(gap2$kmopt$cluster, rep(1:4, n))

# visualization on full dimentions
datagroup <- data.frame(datasetslist2$data_norm1, rep(1:4, n))
pairs(datagroup[,-5], col = datagroup$rep.1.4..n., main = "Data in original groups")

gap22 <- gapnc(datasetslist2$data_norm1, spaceH0 = "scaledPCA", SE.factor = 1)
b22 <- data.frame(datasetslist2$data_norm1, gap22$kmopt$cluster)
colnames(b22)[5] <- "cluster"
centers22 <- as.data.frame(cbind(gap22$kmopt$centers,c(NA,NA,NA)))
colnames(centers22)[5] <- "cluster"
b222 <- rbind(b22, centers22)
pairs(b222[,-5], col = b222$cluster, main = "Data clusterized")

library(GGally)
ggpairs(datagroup[,-5], aes(col = as.factor(datagroup$rep.1.4..n.)))
ggpairs(b222[,-5], aes(col = as.factor(b222$cluster)))

# accuracy
correctly_classified <- sum(rep(1:4, n) == gap22$kmopt$cluster)
tot_obs <- 500
perc <- (correctly_classified/tot_obs)*100
cat("Percentage Correctly Classified:", perc, "%\n")


# 3 ----------------------------------------------------------------
# step 1: single case #
library(cluster)
gapnc <- function(data,FUNcluster=kmeans,
                  K.max=10, B = 100, d.power = 2,
                  spaceH0 ="scaledPCA",
                  method ="globalSEmax", SE.factor = 2,...){ # the three dots allows us to specify other parameters that need a
  # function specification
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

data2 <- read.table("Clusterdata2.txt", header = T)
gap1 <- gapnc(data2, SE.factor = 1, K.max = 20)
gap2 <- gapnc(data2, SE.factor = 2, K.max = 20)
gap1$nc > gap2$nc
cbind(gap1$nc, gap2$nc)

# step 2: general case using 100 datasets#
# generating datasets
library(mvtnorm)
library(fpc)
library(cluster)
datasets2 <- function(n = 100, n.obs = c(50,100,150,200)) {
  lista <- list()
  
  for (i in 1:n) {
    set.seed(i+231344)
    means <- matrix(rbind(sample(20:30, 2), 
                          sample(40:60, 2), 
                          sample(65:85, 2), 
                          sample(90:110, 2)), ncol = 2, byrow = TRUE)
    norm1 <- rmvnorm(n = n.obs[1], mean = means[1, ])
    norm2 <- rmvnorm(n = n.obs[2], mean = means[2, ])
    datalist <- list(norm1, norm2)
    datanorm <- do.call(rbind, datalist)
    datanormdf <- as.data.frame(datanorm)
    name <- paste("data_norm", i, sep = "")
    lista[[name]] <- as.data.frame(datanormdf)
  }
  return(lista)
}
datasetslist2 <- datasets2(n = 100)

plots2 <- function(data) {
  vec1 <- vector()
  vec2 <- vector()
  result <- data.frame()
  
  for (i in 1:length(data)) {
    gap <- gapnc(data[[i]], spaceH0 = "scaledPCA", SE.factor = 1)
    vec1[[i]] <- gap$nc
  }
  for (i in 1:length(data)) {
    gap <- gapnc(data[[i]], spaceH0 = "scaledPCA", SE.factor = 2)
    vec2[[i]] <- gap$nc
  }
  result <- data.frame(vec1,vec2)
}

nc <- plots2(datasetslist2)
nc$vec1 >= nc$vec2
rbind(nc$vec1, nc$vec2)

# 4 ----------------------------------------------------------------
# step1: finding an index #
# index 1, using the relative variation along the two distributions (i think, i hope)
xuno = c(1, 4, 5, 4, 2, 1, 1, 4)
xdue = c(2, 3, 2, 2, 3, 3, 3, 3)
xtre = c(7, 11, 11, 12, 9, 8, 8, 12)
data <- data.frame(xuno,xdue,xtre)
plot(x = 1:8, y = data$xuno, type = "l", col = "green", ylim = c(1,12))
points(x = 1:8, y = data$xdue, type = "l", col = "blue")
points(x = 1:8, y = data$xtre, type = "l", col = "red")
text(x = 3,  y = 2.8,  labels = "x2",  cex = 1, col = "black")
text(x = 4,  y = 5,  labels = "x1",  cex = 1, col = "black")
text(x = 3,  y = 10,  labels = "x3",  cex = 1, col = "black")

distprop <- function(p = 4, x1 = c(1,2,3,4), x2 = c(5,6,7,8)){
  dist <- vector()
  c <- vector()
  d <- vector()
  
  for (i in 1:(p-1)) {
    a <- x1[i+1]-x1[i]+0.1
    b <- x2[i+1]-x2[i]
    c[[i]] <- a
    d[[i]] <- b
    dist <- round(sum(abs(c-d))/(p-1),3)
  }
  return(dist)
}

d12 <- distprop(p = 8, x1 = xuno, x2 = xdue)
d23 <- distprop(p = 8, x1 = xdue, x2 = xtre)
d13 <- distprop(p = 8, x1 = xuno, x2 = xtre)

res <- cbind(d12,d13,d23)
barplot(res, xlab = "Distance index", ylab = "Distance", 
        col = c("red", "blue", "green"), beside = TRUE, ylim = c(0,2), main = "Index 1")
# d23 > d12 > d13

# index 2, using the correlation
distcor <- function(p = 4, x1 = c(1,2,3,4), x2 = c(5,6,7,8)){
  dist = 1 - abs(cor(x1, x2))
  return(dist)
  }
dr12 <- distcor(p = 8, x1 = xuno, x2 = xdue)
dr13 <- distcor(p = 8, x1 = xuno, x2 = xtre)
dr23 <- distcor(p = 8, x1 = xdue, x2 = xtre)
res2 <- cbind(dr12,dr13,dr23)

barplot(res2, xlab = "Distance index", ylab = "Distance", 
        col = c("red", "blue", "green"), beside = TRUE, ylim = c(0,1), main = "Index 2")

# i prefer index 2 because it stays in interval (0,1)

# step 2: prove that it is a dissimilarity #
# 1) d(x; y) = d(y; x) = 0 and
dr21 <- distcor(p = 8, x1 = xdue, x2 = xuno)
dr12 == dr21 # TRUE
# 2) d(x; x) = 0 for x,y in X.
dr11 <- distcor(p = 8, x1 = xuno, x2 = xuno)
dr11 == 0 # TRUE
# 3) triangular equality
# single case
y1 = c(1,3)
y2 = c(5,7)
y3 = c(1,6)
vector <- matrix(c(y1,y2,y3),3,2, byrow = T)
plot(vector, xlim = c(0,6), ylim = c(0,8))
lines(x = c(y1[1], y2[1]),y = c(y1[2], y2[2]), type = "l", col = "red")
lines(x = c(y1[1], y3[1]),y = c(y1[2], y3[2]), type = "l", col = "blue")
lines(x = c(y2[1], y3[1]),y = c(y2[2], y3[2]), type = "l", col = "blue")
text(x = 0.7, y = 3,labels = "y1",cex = 1,  col = "black")
text(x = 5.3,y = 7,labels = "y2",cex = 1,col = "black")
text(x = 0.7,y = 6,labels = "y3",cex = 1,col = "black")
distcor(p = 2, x1 = y1, x2 = y2) + 
  distcor(p = 2, x1 = y1, x2 = y2) >= 
  distcor(p = 2, x1 = y1, x2 = y2) # TRUE

# general case
triangeq <- function(repe = 1000){
  res <- vector()
  n = 3
  v1 <- c(sample(0:1000, 1),sample(0:1000, 1))
  v2 <- c(sample(0:1000, 1),sample(0:1000, 1))
  v3 <- c(sample(0:1000, 1),sample(0:1000, 1))
  vector <- matrix(c(y1,y2,y3),3,2, byrow = T)
  statement <- distcor(p = 2, x1 = v1, x2 = v2) + 
    distcor(p = 2, x1 = v1, x2 = v3) >= 
    distcor(p = 2, x1 = v2, x2 = v3)
  for (i in 1:repe) {
    res[[i]] <- statement
  }
  return(sum(res)/length(res))
}
triangeq(repe = 1000000) # TRUE