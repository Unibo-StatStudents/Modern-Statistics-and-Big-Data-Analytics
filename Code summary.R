##############################################
## MODERN STATISTICS AND BIG DATA ANALYTICS ##
##              code summary                ##
##############################################

# Packages ####
library(pdfCluster) # Datasets
library(MASS)       # Datasets
library(fpc)        # all
library(flexclust)  # graphics
library("smacof")   # MDS
library(nomclust)   # Distances
library(sn)         # k-means
library(cluster)    # k-means
library(mclust)     # Mixture
library(teigen)     # Mixture for non Gaussian
library(mixsmsn)    # Mixture for non Gaussian (not in the exam)
library(flexmix)    # Mixture for categorical data
library(fda)        # Functional data
library(funFEM)     # Functional data
library(robustbase) # robust statistics
library(mixtools)   # Showing ellipses defined by cov

# Visualization ####
pairs() #multivariate scatterplot
plot(pca_object) # barplot with PCs
table(clustering,true_cluster) # confusion matrix
biplot(prolive,cex=0.7) # PCA 2 dimentions with axes
heatmap() #visualize categorical matrices
# Data preprocessing ####
object_scaled <- scale(object) # Standardise
# MDS & PCA ####
object_pca <- princomp(object)
mds_object <- mds(dist_objecct,ndim=2)
mds_object <- cmdscale(dist_object)
# K-means ####
kmeans_object <- kmeans(object,5,nstart=100)
gap_object <- clusGap(object, kmeans_object,...) # calculates a goodness of clustering measure, the “gap” statistic
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
# The output of clusGap is in component gapout.
# The optimal number of clusters is in component nc.
# The optimal kmeans output is in component kmopt.
# Methods for number of clusters and diagnosys ####
sk <- numeric(0)                                     #
kclusterings <- list()                               #
for (k in 1:10){                                     #
  kclusterings[[k]] <- kmeans(object,k,nstart=100)   # Elbow plot
  sk[k] <- kclusterings[[k]]$tot.withinss            # 
}                                                    #
plot(1:10,sk,xlab="k",ylab="S_k",type="l")           #
adjustedRandIndex(clustering1, clustering2) # difference between clustering
pasw <- NA                                                    # 
pclusk <- list()                                              #
psil <- list()                                                #
# Look at K between 2 and 30:                                 #
for (k in 2:30){                                              #
  # PAM clustering:                                           #
  pclusk[[k]] <- pam(p05manhattan,k)                          #
  # Computation of silhouettes:                               # Average Silhouette Width
  psil[[k]] <- silhouette(pclusk[[k]],dist=p05manhattan)      #
  # ASW needs to be extracted:                                #
  pasw[k] <- summary(psil[[k]])$avg.width                     #
}                                                             #
# Plot the ASW-values against K:                              #
plot(1:30,pasw,type="l",xlab="Number of clusters",ylab="ASW") #
plot(psil[[number_of_cluster]])                               #
# Distances and Hierarchical methods ####
dist_object <- dist(matrix, "euclidean/manhattan")  # euclidean(squared)/manhattan(cityblock)
mahalm <- matrix(0,ncol=572,nrow=572)                         #
datacov <- cov(data)                                          #
for (i in 1:572)                                              # Mahalanobis distance
  mahalm[i,] <- mahalanobis(data,as.numeric(data[i,]),datacov)#
mahalanobis <- as.dist(mahalm)                                #
jaccard <- dist(cateforical_data,method="binary") # jaccard distance
simple_matching <- dist(categorical_data,method="manhattan")/583 # simple matching
simple_matching <- sm(categorical_data) # simple matching
gower_dissimilarity <- daisy(mixed_type_data,c("euclidean", "manhattan", "gower"),
                             type=list(symm=4)) # gower dissimilarity
# Dissimilarities will be computed between the rows of x. Columns of mode numeric 
# (i.e. all columns when x is a matrix) will be recognized as interval scaled variables, columns 
# of class factor will be recognized as nominal variables, and columns of class ordered will be 
# recognized as ordinal variables. Other variable types should be specified with the type argument. 
# Missing values (NAs) are allowed.
# type = "asymm","symm","factor","ordered","logratio","ordratio","numeric"/"integer"
hclust <- hclust(dist_matrix,method="average")                   #
# method = "ward.D", "ward.D2", "single", "complete", "average", #
#          "mcquitty", "median" or "centroid"                    #
plot(hclust)                                                     #
eight_clusters_grouping <- cutree(hclust,8)                      # hierarchical clustering
plot(mds$conf,col=,pch=eight_clusters_gorouping)                 #
heatmap(data_matrix,Rowv=as.dendrogram(clust1),                  #
        Colv=as.dendrogram(clust2),                              #
        col=grey(seq(1,0,-0.01)))                                #

# PAM ####
pam_clustering <- pam(data_matrix,5)
# By default, if pam is called with a data set that is not a
# dist-object, the Euclidean distance is used
# Mixture Models ####
GMM <- Mclust(data,G=1:15, modelNames="VVV") #
summary(GMM)                                 # GMM
summary(GMM$BIC)                             #
tolive <- teigen(olive, Gs=1:12)                   #
# data automatically scaled                        # 
# (scale = F to deny it)                           # t-Student
plot(tolive,xmarg=1,ymarg=5,what="contour")        #
plot(tolive,xmarg=1,ymarg=5,what="uncertainty")    #
estg <- smsn.search(scale(olive),nu=1,g.min=1,g.max=12,family="Skew.normal") # Skewed distrib. 
# Mixture Models for categorical data ####
mixture_categorical <- flexmixedruns(data,continuous=0,discrete=583,n.cluster=1:10)
plot(1:10,mixture_categorical$bicvals,typ="l",
     xlab="Number of clusters",ylab="BIC")
jaccard_distance <- dist(data,method="binary")
mds <- mds(jaccard_distance)
plot(mds$conf,col=mixture_categorical$flexout[[6]]@cluster,
     pch=clusym[mixture_categorical$flexout[[6]]@cluster])
matrix <- as.matrix(data)
jaccard_transposed <- dist(t(matrix),method="binary")
Hierarchical_clustering <- hclust(jaccard_transposed,method="average")
# heatmap, rows ordered by clusters,
# columns by earlier variable clustering
heatmap(matrix[order(mixture_categorical$flexout[[6]]@cluster),],
        Rowv=NA,Colv=as.dendrogram(Hierarchical_clustering),
        RowSideColors=palette()[mixture_categorical$flexout[[6]]@cluster]
        [order(mixture_categorical$flexout[[6]]@cluster)],
        col=c(0,1),scale="none")
# Functional data ####
plot(1:555,covid21v[1,],type="l",ylim=c(0,25),              #
     ylab="New cases over one week per 1000 inhabitants",   #  
     xlab="Day (1 April 2020-7 October 2021)",              #
     main="Covid weekly new cases for 179 countries")       # Visualization
for(i in 2:179)                                             #
  points(1:555,covid21v[i,],type="l")                       #
bbasis <- create.bspline.basis(c(1,555),nbasis=100)               # 
fdcovid <- Data2fd(1:555,y=t(as.matrix(covid21v)),basisobj=bbasis)#
# Plot basis                                                      #
plot(bbasis)                                                      #
# Smooth splines for data with smooth mean function:              #
plot(fdcovid)                                                     # B-spline basis
mcovid <- mean.fd(fdcovid)                                        #
lines(mcovid,col=2,lwd=5)                                         #
# Show smooth fit of individual countries                         #
plotfit.fd(t(covid21v),1:555,fdcovid10,index=79,cex.pch=0.5)      # 
covidpca <- pca.fd(fdcovid, nharm = 5)                        #
plot(covidpca$harmonics) # PCs phi_k                          #
covidpca$varprop # Percentage of variance                     #
cumsum(covidpca$varprop) # Cumulative percentage of variance  # fda PCA
ncontinent <- as.numeric(as.factor(covid21$continent))        #
levels(as.factor(covid21$continent))                          #
pairs(covidpca$scores,col=ncontinent,pch=ncontinent)          #
# PCA scores                                                  #
covidpcaapprox <- covidpca$harmonics                                       #
i <- 1                                                                     #
pcacoefi <- covidpca$harmonics$coefs %*% covidpca$scores[i,]+mcovid$coefs  #
covidpcaapprox$coefs <- pcacoefi                                           #
for (i in 2:179){                                                          #
  pcacoefi <- covidpca$harmonics$coefs %*% covidpca$scores[i,]+mcovid$coefs# Create functional data object of PCA approximations
  covidpcaapprox$coefs <- cbind(covidpcaapprox$coefs, pcacoefi)            #
}                                                                          #
dimnames(covidpcaapprox$coefs)[[2]] <- covid21[,1]                         #
plotfit.fd(t(covid21v),1:555,covidpcaapprox,index=79,cex.pch=0.5)          #

covidcluster <- funFEM(fdcovid,K=2:10)
femmodels <- c("DkBk", "DkB", "DBk", "DB", "AkjBk", "AkjB", 
               "AkB", "AkBk", "AjBk", "AjB", "ABk", "AB")
nmodels <- length(femmodels)
femresults <- list()
bestk <- bestbic <- numeric(0)
# bestk: vector of best K for each model.
# bestbic: Best BIC value for each model.
K=2:10 # Numbers of clusters K to try out.
fembic <- matrix(NA,nrow=nmodels,ncol=max(K))
# fembic will hold all BIC values for models (rows) and K (columns);
# NA for those that cannot be fitted.
for (i in 1:nmodels){
  print(femmodels[i])
  femresults[[i]] <- funFEM(fdcovid,model=femmodels[i],K=K)
  fembic[i,K] <- femresults[[i]]$allCriterions$bic
  bestk[i] <- which(fembic[i,]==max(fembic[i,K],na.rm=TRUE))
  bestbic[i] <- max(fembic[i,K],na.rm=TRUE)
}
besti <- which(bestbic==max(bestbic,na.rm=TRUE))
besti

femmodels[besti] # Sigma_k spherical and equal, beta_k can vary with k.
bestk # K=8 optimal for model 11 "ABk"
femresult11 <- femresults[[11]]
# Plot BIC values for all models and K:
i <- 1
plot(1:max(K),fembic[i,],col=i,pch=i, 
     ylim=c(min(fembic,na.rm=TRUE),max(fembic,na.rm=TRUE)),type="n")
for(i in 1:nmodels)
  text(1:max(K),fembic[i,],femmodels[i],col=i)

# Plot clusters against latitude and longitude
# This prints out the countries in the clusters.
for(i in 1:femresult11$K){
  print(i)
  print(covid21[femresult11$cls==i,1])
}

# Clusters on principal components
pairs(covidpca$scores,col=femresult11$cls,pch=19)
# Visualisation of discriminative subspace U,
# projection of observations on U-space:
fdproj <- t(fdcovid$coefs) %*% femresult11$U
pairs(fdproj,col=femresult11$cls,pch=19)
plot(fdproj,col=femresult11$cls,pch=19,xlab="DC 1",ylab="DC 2")

# Plot the curves and clusters
plot(1:555,covid21v[1,],type="l",ylim=c(0,25),
     ylab="New cases over one week per 1000 inhabitants",
     xlab="Day (1 April 2020-7 October 2021)",lwd=1.5,col=femresult11$cls[1])
for(i in 2:179)
  points(1:555,covid21v[i,],type="l",col=femresult11$cls[i],lwd=1.5)
# Plot the cluster mean curves
clmeans <- fdcovid
clmeans$coefs <- t(femresult11$prms$my)
plot(clmeans,lwd=3) # col doesn't seem to work here, neither lwd
legend(100,10,legend=1:8,col=c(1:6,1:2),lty=c(1:5,1:3))
# Plot individual clusters and mean curves
par(ask=TRUE)
for (k in 1:femresult11$K){
  plot(1:555,covid21v[1,],type="l",ylim=c(0,25),
       ylab="New cases over one week per 1000 inhabitants",
       xlab="Day (1 April 2020-7 October 2021)",col=as.integer(femresult11$cls[1]==k))
  for(i in 2:179)
    points(1:555,covid21v[i,],type="l",col=as.integer(femresult11$cls[i]==k))
  meank <- colMeans(covid21v[femresult11$cls==k,])
  points(1:555,meank,type="l",lwd=5,col=2)
}
par(ask=FALSE)

# Roubust statistics ####
mad(dortmund$birthdeath, constant = 1.4826) # Median Absolute Deviation
huberM(dortmund$birthdeath,k=1.5) # Huber M-estimator
mcdd <- covMcd(dortmund, alpha=0.75) #
colMeans(dortmund)                   # Minimum covariance determinant estimator
mcdd$center                          # 
plot(1:170,sqrt(mcdd$mah),type="n",xlab="Observation",                         #
     ylab="Robust Mahalanobis distance")                                       #
text(1:170,sqrt(mcdd$mah),rownames(dortmund),cex=0.7)                          #
abline(sqrt(qchisq(0.99,7)),0,col=2)                                           #
# Should look at smaller values to see more precisely what's going on:         #
plot(1:170,sqrt(mcdd$mah),type="n",ylim=c(0,30),xlab="Observation",            #
     ylab="Robust Mahalanobis distance")                                       #
text(1:170,sqrt(mcdd$mah),rownames(dortmund),cex=0.7)                          #
abline(sqrt(qchisq(0.99,7)),0,col=2)                                           #
#
# Compare with Mahalanobis distances based on mean and sample cov-matrix       #
cd <- cov(dortmund)                                                            #
plot(sqrt(mahalanobis(dortmund,colMeans(dortmund),cd)),sqrt(mcdd$mah),type="n",# 
     xlab="Standard Mahalanobis distance",                                     # Outliers
     ylab="Robust Mahalanobis distance")                                       # Identification
text(sqrt(mahalanobis(dortmund,colMeans(dortmund),cd)),sqrt(mcdd$mah),         #
     rownames(dortmund),cex=0.7)                                               #
abline(sqrt(qchisq(0.99,7)),0,col=2)                                           #
abline(v=sqrt(qchisq(0.99,7)),col=2)                                           #
# Look closer                                                                  #
plot(sqrt(mahalanobis(dortmund,colMeans(dortmund),cd)),sqrt(mcdd$mah),type="n",#
     xlim=c(0,10),ylim=c(0,30),xlab="Standard Mahalanobis distance",           #
     ylab="Robust Mahalanobis distance")                                       # 
text(sqrt(mahalanobis(dortmund,colMeans(dortmund),cd)),sqrt(mcdd$mah),         #
     rownames(dortmund),cex=0.7)                                               #
abline(sqrt(qchisq(0.99,7)),0,col=2)                                           #
abline(v=sqrt(qchisq(0.99,7)),col=2)                                           #
