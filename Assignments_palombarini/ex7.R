library(teigen)
library(mixsmsn)
library(cluster)
library(smacof)
library(mclust)
library("robustbase")
library(flexclust)
library(MASS)
library(fpc)
library(poLCA)
library(nomclust)


# 1 ----------------------------------------------------------------------------------------------
data(election)

# Survey data from the 2000 American National Election Study. Two sets of six questions with four 
# responses each, asking respondents' opinions of how well various traits (moral, caring, knowledgable, 
# good leader, dishonest, intelligent) describe presidential candidates Al Gore and George W. Bush. 
# The responses are (1) Extremely well; (2) Quite well; (3) Not too well; (4) Not well at all. Many 
# respondents have varying numbers of missing values on these variables.
# The data set also includes potential covariates VOTE3, the respondent's 2000 vote choice (when asked); 
# AGE, the respondent's age; EDUC, the respondent's level of education; GENDER, the respondent's gender; 
# and PARTY, the respondent's Democratic-Republican partisan identification.
# 
# VOTE3 is coded as (1) Gore; (2) Bush; (3) Other.
# 
# EDUC is coded as (1) 8 grades or less; (2) 9-11 grades, no further schooling; (3) High school diploma 
# or equivalency; (4) More than 12 years of schooling, no higher degree; (5) Junior or community college 
# level degree; (6) BA level degrees, no advanced degree; (7) Advanced degree.
# 
# GENDER is coded as (1) Male; (2) Female.
# 
# PARTY is coded as (1) Strong Democrat; (2) Weak Democrat; (3) Independent-Democrat; 
# (4) Independent-Independent; (5) Independent-Republican; (6) Weak Republican; (7) Strong Republican.

election12 <- election[,1:12]
electionwithna <- election12
for (i in 1:12){
  levels(electionwithna[,i]) <- c(levels(election12[,i]),"NA")
  electionwithna[is.na(election12[,i]),i] <- "NA"
}
data <- electionwithna

mds <- cmdscale(sm(data), k = 2)
plot(mds, main = "MDS") # skewed gaussian, spherical

# a
# Compute a latent class clustering with 3 clusters using poLCA (read its help page and decide about 
# potentially useful parameter settings).

# poLCA(formula, data, nclass = 2, maxiter = 1000, graphs = FALSE, 
#       tol = 1e-10, na.rm = TRUE, probs.start = NULL, nrep = 1, 
#       verbose = TRUE, calc.se = TRUE)
# Arguments
# formula: A formula expression of the form response ~ predictors. The details of model specification 
# are given below.
# data:	A data frame containing variables in formula. Manifest variables must contain only integer values, 
# and must be coded with consecutive values from 1 to the maximum number of outcomes for each variable. 
# All missing values should be entered as NA.
# nclass: The number of latent classes to assume in the model. Setting nclass=1 results in poLCA estimating 
# the loglinear independence model. The default is two.
# maxiter: The maximum number of iterations through which the estimation algorithm will cycle.
# graphs:	 Logical, for whether poLCA should graphically display the parameter estimates at the completion of 
# the estimation algorithm. The default is FALSE.
# tol: A tolerance value for judging when convergence has been reached. When the one-iteration change 
# in the estimated log-likelihood is less than tol, the estimation algorithm stops updating and considers 
# the maximum log-likelihood to have been found.
# na.rm: Logical, for how poLCA handles cases with missing values on the manifest variables. If TRUE, those 
# cases are removed (listwise deleted) before estimating the model. If FALSE, cases with missing values are 
# retained. Cases with missing covariates are always removed. The default is TRUE.
# probs.start: A list of matrices of class-conditional response probabilities to be used as the starting 
# values for the estimation algorithm. Each matrix in the list corresponds to one manifest variable, with one
# row for each latent class, and one column for each outcome. The default is NULL, producing random starting 
# values. Note that if nrep>1, then any user-specified probs.start values are only used in the first of the
# nrep attempts.
# nrep: Number of times to estimate the model, using different values of probs.start. The default is one. 
# Setting nrep>1 automates the search for the global-rather than just a local-maximum of the log-likelihood 
# function. poLCA returns the parameter estimates corresponding to the model with the greatest log-likelihood.
# verbose: Logical, indicating whether poLCA should output to the screen the results of the model. If FALSE, 
# no output is produced. The default is TRUE.
# calc.se: Logical, indicating whether poLCA should calculate the standard errors of the estimated 
# class-conditional response probabilities and mixing proportions. The default is TRUE; can only be set 
# to FALSE if estimating a basic model with no concomitant variables specified in formula.
f <- cbind(MORALG,CARESG,KNOWG,LEADG,DISHONG,
           INTELG,INTELG,MORALB,CARESB,KNOWB,
           LEADB,DISHONB,INTELB)~1
Ma <- poLCA(formula = f, data = data, nclass = 3)
summary(Ma)
plot(mds, col = Ma$predclass, pch = 20, main = "MDS with poLCA")

# b
# Compute a latent class clustering with 3 clusters using flexmixedruns.
Mb <- flexmixedruns(data, continuous = 0, discrete = 12, n.cluster = 3)
Mb$bicvals
plot(mds,col= Mb$flexout[[3]]@cluster , pch=20, main = "MDS with flexmixedruns")


BIC_comparison <- as.data.frame(rbind(cbind("Model", "BIC"),cbind("Ma", Ma$bic), 
                                      cbind("Mb", Mb$optsummary@BIC)))
colnames(BIC_comparison) <- c("Model", "BIC")
BIC_comparison <- BIC_comparison[-1,]
BIC_comparison[,2] <- round(as.numeric(BIC_comparison[,2]), 2)
BIC_comparison


# c
# Compute a distance-based clustering of your choice with 3 clusters based on the simple matching distance 
# (you can also choose here between computing the simple matching distance on electionwithna or on election12, 
# which will compute the distance just taking variables into account that are non-missing on both observations, 
# using daisy as explained earlier in class - actually in general then it's a dissimilarity and not a distance 
# as missing values can spoil the triangle inequality).
Mc <- pam(sm(data), k = 3)
summary(Mc)
plot(mds, col = Mc$clustering, pch = 20, main = "MDS with PAM")

# ASW to find the right K
pasw <- NA
pclusk <- list()
psil <- list()
# Look at K between 2 and 30:
for (k in 2:20){
  # PAM clustering:
  pclusk[[k]] <- pam(sm(data),k)
  # Computation of silhouettes:
  psil[[k]] <- silhouette(pclusk[[k]],dist=sm(data))
  # ASW needs to be extracted:
  pasw[k] <- summary(psil[[k]])$avg.width
}
# Plot the ASW-values against K:
plot(1:20,pasw,type="l",xlab="Number of clusters",ylab="ASW", main = "ASW plot")
abline(v = c(1:10), lty = 2, col = "grey")

# 2, 3 and 6 are all good levels for k

# d
# Compute a latent class clustering using flexmixedruns with estimated number of clusters. (poLCA will 
# not estimate the number of cluster automatically, although it gives out the BIC, so this could in 
# principle be implemented easily.)
Md <- flexmixedruns(data, continuous = 0, discrete = 12, n.cluster = 1:10)
Md$optimalk
Md$bicvals[8]
plot(mds, col = Md$flexout[[8]]@cluster, main = "MDS with flexmixedruns 2", pch = 20)

# e
# The original election data set also has the variables AGE and EDUC (these
# are the variables number 14 and 15). Define a data set election14 that has
# the 12 variables already used above and AGE and EDUC. Assuming that these
# can be treated as continuous variables (which is rather questionable at least
# for EDUC, but you can ignore this here), use flexmixedruns to compute
# a clustering based on a latent class mixture model, including estimation of
# the number of clusters, were the continuous variables are assumed Gaussian
# within clusters independently of the categorical variables. In order to use all
# observations for this, impute the missing values of AGE and EDUC with the
# mean of these variables. You can use the MDS already computed above for
# visualising this, but if you are curious, you can also run a new MDS for these
# data using the Gower coefficient.

election14 <- election[,14:15]
electionwithna2 <- election14
for (i in 1:2) {
  electionwithna2[is.na(election14[,i]),i] <- mean(na.omit(election14[,i]))
}
data.e <- as.data.frame(cbind(electionwithna2,data))
for (i in 3:14) {
  data.e[,i] <- as.numeric(data.e[,i])
}

mds.e <- cmdscale(sm(data.e), k = 2)
plot(mds.e, main = "MDS 2")

Me <- flexmixedruns(data.e, continuous = 2, discrete = 12, n.cluster = 1:10)
Me$optimalk
plot(mds.e, col = Me$flexout[[7]]@cluster, main = "MDS by flexmixedruns 3")


# 2 ----------------------------------------------------------------------------------------------
datam <- data
for (i in 1:ncol(datam)) {
  levels(datam[,i]) <- c(1,2,3,4,5)
}
str(datam)
for (i in 1:ncol(datam)) {
  datam[,i] <- as.integer(datam[,i])
}
str(datam)
Md <- flexmixedruns(datam, continuous = 0, discrete = 12, n.cluster = 1:10)
datam <- as.matrix(datam)
Mdhclust <- hclust(dist(t(datam)), method = "average")
plot(Mdhclust)

# Do you nd the clusters convincing? Why or why not? Is there evidence against local independence?

heatmap(datam[order(Md$flexout[[8]]@cluster),], Rowv=NA, Colv=as.dendrogram(Mdhclust),
        RowSideColors=palette()[Md$flexout[[8]]@cluster][order(Md$flexout[[8]]@cluster)], 
        col=c("lightblue", "lightblue", "orange", "orange", "white"), scale="none",
        main = "Heatmap for latent class clustering with flexmixedruns")
legend_labels <- c("Positive Judgment", "Negative Judgment", "No Answer")
legend_colors <- c("lightblue", "orange", "white")
legend("topright", legend = legend_labels, fill = legend_colors, title = "", bty = "n", cex = 0.6)
mtext("Plot seems to underline how dishonesty isn't considered 
a characteristic that describes one of the two candidates, 
exception made for the first cluster, the black one, that 
thinks Gore is enough dishonest.
We also notice that lightblue cluster (number 5) contains almost 
all the NA observations.
While yellow and blue cluster have, in general, a good opinion 
of both the  candidates, violet one has a bad opinion of both.
To conclude, black and green clusters seem to prefer Gore and, 
while red and grey ones seem to prefer Bush.

I personally find the clustering convincing, because the 
distinctions are  pretty clear. 
Data really seem locally dependent", cex = 0.9, at = 0.77, adj = 0, padj = 1.6)

heatmap(datam[order(Md$flexout[[8]]@cluster),], Rowv=NA, Colv=as.dendrogram(Mdhclust),
        RowSideColors=palette()[Md$flexout[[8]]@cluster][order(Md$flexout[[8]]@cluster)], 
        col=c("blue", "white", "white", "red", "white"), scale="none",
        main = "Heatmap for latent class clustering with 
        flexmixedruns ('Strong opinions')")

heatmap(datam[order(Ma$predclass),], Rowv=NA, Colv=as.dendrogram(Mdhclust),
        RowSideColors=palette()[Ma$predclass][order(Ma$predclass)], 
        col=c("lightblue", "lightblue", "orange", "orange", "white"), scale="none",
        main = "Heatmap for latent class clustering with poLCA")

# 3 --------------------------------------------------------------------------------------------------
# Assume a situation with 10 categorical variables. Five variables are
# binary, three variables have three categories, and two variables have five categories.
# What is the number of free parameters for
# (a) a general categorical model that models all possible probabilities,
# (b) a latent class mixture model with 4 mixture components?
#   Remark: Note that the value computed in (a) plus one is the number of possible
# dierent observations! If this number is low, it also implies a strong restriction for
# the possible numbers of clusters to be tted, which cannot be larger and should
# normally be substantially smaller, regardless of the number of observations.
3+(1*3+2*3+4*3)=3+3+6+12
5+6+8

# 4 ----------------------------------------------------------------------------------------------
# Consider the COVID data set analysed in Chapter 7 of the course slides.
# Consider Italy, Haiti, and the USA (country 79, 69, and 164). Produce residual
# plots, i.e., plots with the time points on the x-axis and the residuals (difference
# between data and fit) on the y-axis for these countries for
# (a) the fit by a B-spline basis with p = 100,
# (b) the fit by a 5-dimensional principal components basis as shown on p. 310 of
# the slides.
# Comment on potential model assumption violations from these plots.

library(fda)
covid21 <- read.table("Covid2021.txt")
covid21v <- as.matrix(covid21[c(69,79,164),5:559])
# Raw data plot:
plot(1:555,covid21v[1,],type="l",ylim=c(0,25),ylab="New cases over one
week per 1000 inhabitants",xlab="Day (1 April 2020-7 October 2021)",
     main="Covid weekly new cases for 3 countries")
for(i in 1:3)
  points(1:555,covid21v[i,],type="l")

bbasis <- create.bspline.basis(c(1,555),nbasis=100)
fdcovid <- Data2fd(1:555,y=t(as.matrix(covid21v)),basisobj=bbasis)
plot(fdcovid)
mcovid <- mean.fd(fdcovid)
lines(mcovid,col="blue4",lwd=3)

fit_values <- eval.fd(1:555, fdcovid)
residbspl <- covid21v - t(fit_values)
for (i in 1:3) {
  plot(1:555, residbispl[i,], type="l", 
       ylab="Residuals", xlab="Time",
       main=paste("Residual Plot for", rownames(covid21v)[i]))
}


# Fit a 5-dimensional principal components basis
pca_basis <- pca.fd(fdcovid, nharm = 5)
pca_scores <- pca_basis$scores
pca_functions <- pca_basis$harmonics
eval_pca_functions <- eval.fd(1:555, pca_functions)
fit_values_pca <- pca_scores %*% t(eval_pca_functions)
residuals_pca <- covid21v - fit_values_pca

for (i in 1:3) {
  plot(1:555, residuals_pca[i,], type = "l",
       ylab = "Residuals", xlab = "Day",
       main = paste("Residual Plot for", rownames(covid21v)[i], "using PCA"))
}

# 5 ------------------------------------------------------------------------------------------------
# Representing all countries in the COVID data set by the first functional
# principal component scores only, run a one-way analysis of variance to test whether
# there is evidence that the scores from different continents have different means. Also
# visualise the scores so that the continents can easily be compared, and interpret plot
# and result (try to figure out what larger or smaller/positive or negative scores on
# the first principal component actually mean). This is a simply method to run a test
# comparing groups of functional data objects (obviously relying on the information
# represented by the first principal component only).

covid215 <- as.matrix(covid21[,5:559])
fdcovid5 <- Data2fd(1:555, y=t(as.matrix(covid215)), basisobj=bbasis)
mcovid5 <- mean.fd(fdcovid5)
covidpca5 <- pca.fd(fdcovid5, nharm = 1)
pca_scores <- covidpca5$scores[, 1]
covid_data <- data.frame(
  Country = rownames(covid215),
  Continent = covid21[,2], 
  PCA_Score = pca_scores
)

str(covid21)
anova_result <- aov(PCA_Score ~ Continent, data = covid_data)
summary(anova_result)
boxplot(PCA_Score ~ Continent, data = covid_data,
        ylab = "First Principal Component Scores",
        xlab = "Continent",
        main = "Scores by Continent", cex.axis = 0.7)
scores <- data.frame(cbind(covid21[,2], as.vector(pca_scores)))
scores$X2 <- as.numeric(scores$X2)
scores.average <- data.frame(as.vector(levels(as.factor(covid21$continent))),
                             c(sum(scores$X2[scores$X1=="Africa"])/sum(scores$X1=="Africa"),
                               sum(scores$X2[scores$X1=="Asia"])/sum(scores$X1=="Asia"),
                               sum(scores$X2[scores$X1=="Australia"])/sum(scores$X1=="Australia"),
                               sum(scores$X2[scores$X1=="Europe"])/sum(scores$X1=="Europe"),
                               sum(scores$X2[scores$X1=="North America"])/sum(scores$X1=="North America"),
                               sum(scores$X2[scores$X1=="South America"])/sum(scores$X1=="South America"),
                               sum(scores$X2[scores$X1=="Central America"])/sum(scores$X1=="Central America")))
colnames(scores.average) <- c("Continents", "Scores Averages")
scores.average
































