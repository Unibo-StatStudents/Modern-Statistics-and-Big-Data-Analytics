# 1 ---------------------------------------------------------------------------
phonemes1000 <- read.table("phonemes1000.txt",header=TRUE)
phonemes1000 <- phonemes1000[1:200,]
data <- as.matrix(phonemes1000[1:200,1:256])

# Represent the data in terms of a suitable B-spline basis. 
library(fda)
library(mclust)
library(sn)
library(cluster)
plot(1:256,data[1,],type="l",ylim=c(0,25), main = "Dataset plot")
for(i in 2:200)
  points(1:256, data[i,],type="l")
i <- 1
plot(1:256, data[i,],cex=0.5,main = colnames(data)[i])
i <- 30 
plot(1:256, data[i,],cex=0.5,main = colnames(data)[i])
# Constructing B-spline basis
bbasis <- create.bspline.basis(c(1,256),nbasis=80)
fdcovid <- Data2fd(1:256,y=t(as.matrix(data)),basisobj=bbasis)
# Plot basis
plot(bbasis, main = "Basis plot")
# Smooth splines for data with smooth mean function:
plot(fdcovid, main = "Smooth splines for data with smooth mean function")
mcovid <- mean.fd(fdcovid)
lines(mcovid,col=2,lwd=5)
# Show how well two exemplary observations are approximated in this way.
plotfit.fd(t(data),1:256,fdcovid,index=1,cex.pch=0.5, main = "smooth fit of 1st obs.")
plotfit.fd(t(data),1:256,fdcovid,index=30,cex.pch=0.5, main = "smooth fit of 30th obs.")

# Run a functional principal component analysis
pca <- pca.fd(fdcovid, nharm = 5)
pca$varprop # Percentage of variance
g <- as.numeric(as.factor(phonemes1000$g))
levels(as.factor(phonemes1000$g))
pairs(pca$scores,col=g,pch=g, main = "PCA scores")
# PCA scores

# Create functional data object of PCA approximations
pcaapprox <- pca$harmonics
i <- 1
pcacoefi <- pca$harmonics$coefs %*% pca$scores[i,] + mcovid$coefs
pcaapprox$coefs <- pcacoefi
for (i in 2:200){
  pcacoefi <- pca$harmonics$coefs %*% pca$scores[i,]+mcovid$coefs
  pcaapprox$coefs <- cbind(pcaapprox$coefs, pcacoefi)
}


# Run a funFEM-clustering 
library(funFEM)
cluster <- funFEM(fdcovid,K=2:6)
set.seed(1234567)
femmodels <- c("DkBk", "DkB", "DBk", "DB", "AkjBk", "AkjB", 
               "AkB", "AkBk", "AjBk", "AjB", "ABk", "AB")
nmodels <- length(femmodels)
femresults <- list() # Save output for all models in femmodels
bestk <- bestbic <- numeric(0)
K=2:6 # Numbers of clusters K to try out.
fembic <- matrix(NA,nrow=nmodels,ncol=max(K))
for (i in 1:nmodels){
  print(femmodels[i])
  femresults[[i]] <- funFEM(fdcovid,model=femmodels[i],K=K)
  fembic[i,K] <- femresults[[i]]$allCriterions$bic
  bestk[i] <- which(fembic[i,]==max(fembic[i,K],na.rm=TRUE))
  bestbic[i] <- max(fembic[i,K],na.rm=TRUE)
}
besti <- which(bestbic==max(bestbic,na.rm=TRUE))
besti
femmodels[besti]
bestk # K=6 optimal for model 11 "ABk"
femresult11 <- femresults[[11]]
# Plot BIC values for all models and K:
i <- 1
plot(1:max(K),fembic[i,],col=i,pch=i, 
     ylim=c(min(fembic,na.rm=TRUE),max(fembic,na.rm=TRUE)),type="n",
     main = "BIC plot")
for(i in 1:nmodels)
  text(1:max(K),fembic[i,],femmodels[i],col=i)
pairs(pca$scores,col=femresult11$cls,pch=19)

# Plot the cluster mean curves
clmeans <- fdcovid
clmeans$coefs <- t(femresult11$prms$my)
plot(clmeans,lwd=3, xlim = c(0,320), main = "Clusters' mean curves")
legend(280,20,legend=1:8,col=c(1:6,1:2),lty=c(1:5,1:3))

# Run a cluster analysis of your choice on the functional principal component scores.
plot(pca$scores) # let's use a GMM
library(mclust)
gmm <- Mclust(data,G=1:20)
summary(gmm)
# 8 clusters
summary(gmm$BIC)
layout(mat = c(1,2))
layout(mat = c(1,1))

# Comparison 
adjustedRandIndex(gmm$classification, phonemes1000$g)
plot(pca$scores,col=gmm$classification,pch = 19)
table(gmm$classification, phonemes1000$g)

adjustedRandIndex(femresult11$cls, phonemes1000$g)
plot(pca$scores,col=femresult11$cls,pch = 19)
table(femresult11$cls, phonemes1000$g)


library(fpc)
?plotcluster()
# 3 ---------------------------------------------------------------------------
# (a) Explain in your own words what discriminant coordinates (\dc") and asymmetric 
# weighted discriminant coordinates (\awc") are, and how they work.

# (b) For the optimal 10-clusters Gaussian mixture clustering of the olive oil data
# (Example 6.5 on the course slides) show 2-dimensional discriminant coordi-
#   nates, and asymmetric weighted discriminant coordinates for all clusters1.
# Comment on how these plots compare to the principal components plot in
# terms of showing the separation of the clusters.
 
# (c) Considering the phoneme data from question 1 and the funFEM-clustering, com-
#   pare the plot of the rst two dimensions of the Fisher discriminating subspace
# from funFEM with what you get when applying discriminant cordinates using
# plotcluster to the data set of the coecients of the full dimensional B-spline
# basis and the funFEM-clustering.


fisher <- cluster$U
fbcoef <- as.matrix(pcaapprox$coefs)
# Apply discriminant coordinates using plotcluster
library(fpc)
plotcluster(x = fbcoef, clvecd = g, method = "dc", main = "Discriminant Coordinates Plot")
# Plot Fisher discriminating subspace and add points from discriminant coordinates plot
plot(fisher, col = g, pch = 19, main = "Fisher Discriminating Subspace")
points(pca$scores[, 1:2], col = g, pch = 17)

