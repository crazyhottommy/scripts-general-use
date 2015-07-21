##########################################################
#Dataset 2 - NCI Microarray Data
#Understand PCA and Sparse PCA
#PCA solution via the SVD
###########################################################

require("ISLR")

ncidat = t(NCI60$data)
colnames(ncidat) = NCI60$labs

dim(ncidat)
unique(colnames(ncidat))

## this is a pretty old microarray data set. it contains gene expression profile for different cancer
## types.

head(ncidat)

#PCA - take SVD to get solution
#center genes, but don't scale
## it is important to know whether you need to center or scale your data set.

## we are going to center the same gene from different samples, but not scale them
## that's why we transpose the matrix (scale works on columns  ?scale), and then transpose it back.
## be aware of whether you should center/scale your data or not.

## one has to be aware that in the microarray data, columns are samples(observations n), rows are genes(features p)
## usually for a SVD analysis
## X is a n x p matrix (n rows and p columns)
## Xnxp = Unxn Dnxp Vpxp
## we need to first transpose X for the microarray data for the svd analysis

X = t(scale(t(ncidat),center=TRUE,scale=FALSE))

## use svd decomposition to get the vectors. to get an idea of how svd works, I strongly recommend you
## read this :http://genomicsclass.github.io/book/pages/svd.html
## and this:  http://genomicsclass.github.io/book/pages/mds.html
## one can use the base R function princomp (default center and scale), but svd gives you more controls.

## D is a diagnal matrix, use length(D) to look at the length
## by definition D should be a n x p (64 x 6830) matrix, but in this case, it becomes 64 x 64 matrix.
## the reason is that Diagonals of D: d1 >= d2 >= d3 >= ....d(r) where r =rank(X), rank(X) <= min(n,p)
## the rank of a matrix https://en.wikipedia.org/wiki/Rank_(linear_algebra)
## the other ds are all zeros. (no variations after d(r)), that's why D is dropped to 64 x 64 matrix.

sv = svd(t(X))
U = sv$u
V = sv$v
D = sv$d

## in R calculate the rank of a matrix is by
qr(t(X))$rank

#63
length(D)
min(D)
# it is very close to 0, it has to do with the precision of the decimals in computer


## let's plot the first 4 PCs
# PC scatterplots
cols = as.numeric(as.factor(colnames(ncidat)))
K = 3
pclabs = c("PC1","PC2","PC3","PC4")
par(mfrow=c(1,K))
for(i in 1:K){
  j = i+1
  plot(U[,i],U[,j],type="n",xlab=pclabs[i],ylab=pclabs[j])
  text(U[,i],U[,j],colnames(X),col=cols)
}


## U are un-scaled PCs
## we can also plot Z which is scaled PC
## Z = XV or Z=UD
# or Z<- U %*% diag(D)

Z = t(X)%*%V

plot(Z[,1], Z[,2], type ="n")
text(Z[,1], Z[,2], colnames(X), col=cols)

pc_dat<- data.frame(type = rownames(Z), PC1 = Z[,1], PC2= Z[,2])

library(ggplot2)
library(tidyr)

ggplot(pc_dat,aes(x=PC1, y=PC2, col=type)) + geom_point() + geom_text(aes(label = type), hjust=0, vjust=0)

## use directlabels http://directlabels.r-forge.r-project.org/
## use cowplot https://github.com/wilkelab/cowplot



#PC loadings - visualize data by limiting to top genes in magnitude in the PC loadings
## the matrix V contains the weigths for the features, and we can use V to select important
## features(genes) that contribute to the PC2

aa = grep("grey",colors())
bb = grep("green",colors())
cc = grep("red",colors())
gcol2 = colors()[c(aa[1:30],bb[1:20],rep(cc,2))]

j = 2
ord = order(abs(V[,j]),decreasing=TRUE)
x = as.matrix(X[ord[1:250],])
heatmap(x,col=gcol2)

#Variance Explained
varex = 0
cumvar = 0
denom = sum(D^2)
for(i in 1:64){
  varex[i] = D[i]^2/denom
  cumvar[i] = sum(D[1:i]^2)/denom
}

#screeplot
par(mfrow=c(1,2))
plot(1:64,varex,type="l",lwd=2,xlab="PC",ylab="% Variance Explained")
plot(1:64,cumvar,type="l",lwd=2,xlab="PC",ylab="Cummulative Variance Explained")

#######
#Sparse PCA
require("PMA")

spc = SPC(t(X),sumabsv=10,K=4)

#how many genes selected?
apply(spc$v!=0,2,sum)

#PC scatterplots
cols = as.numeric(as.factor(colnames(ncidat)))
K = 3
pclabs = c("SPC1","SPC2","SPC3","SPC4")
par(mfrow=c(1,K))
for(i in 1:K){
  j = i+1
  plot(spc$u[,i],spc$u[,j],type="n",xlab=pclabs[i],ylab=pclabs[j])
  text(spc$u[,i],spc$u[,j],colnames(X),col=cols)
}

#SPC loadings - visualize data by limiting to gene selected by the sparse PC loadings
aa = grep("grey",colors())
bb = grep("green",colors())
cc = grep("red",colors())
gcol2 = colors()[c(aa[1:30],bb[1:20],rep(cc,2))]

j = 1
ind = which(spc$v[,j]!=0)
x = as.matrix(X[ind,])
heatmap(x,col=gcol2)

#variance explained
spc$prop.var.explained

########################################
#apply K-means
K = 9
km = kmeans(t(ncidat),centers=K)

#how do we visualize K-means results?

#PCA - take SVD to get solution
#center genes, but don't scale
X = t(scale(t(ncidat),center=TRUE,scale=FALSE))
sv = svd(t(X));
U = sv$u
V = sv$v
D = sv$d
Z = t(X)%*%V;

## plot
plot(Z[,1],Z[,2],col=km$cluster,type="n")
text(Z[,1],Z[,2],colnames(ncidat),cex=.75,col=km$cluster)
cens = km$centers
points(cens%*%V[,1],cens%*%V[,2],col=1:K,pch=16,cex=3)

#Re-run and see if solution changes, set.seed() if you want reproducible result!
## k-means initialize each observation i to a cluster assignment k
K = 9
km = kmeans(t(ncidat),centers=K)
plot(Z[,1],Z[,2],col=km$cluster,type="n")
text(Z[,1],Z[,2],colnames(ncidat),cex=.75,col=km$cluster)
cens = km$centers
points(cens%*%V[,1],cens%*%V[,2],col=1:K,pch=16,cex=3)

#try different K
K = 5
km = kmeans(t(ncidat),centers=K)
plot(Z[,1],Z[,2],col=km$cluster,type="n")
text(Z[,1],Z[,2],colnames(ncidat),cex=.75,col=km$cluster)
cens = km$centers
points(cens%*%V[,1],cens%*%V[,2],col=1:K,pch=16,cex=3)

##################
#hierarchical clustering

require("ISLR")

ncidat = t(NCI60$data)
colnames(ncidat) = NCI60$labs

dim(ncidat)
unique(colnames(ncidat))

#complete linakge - Euclidean distance
cols = as.numeric(as.factor(colnames(ncidat)))
Dmat = dist(t(ncidat))
com.hclust = hclust(Dmat,method="complete")
plot(com.hclust,cex=.7,main="Complete Linkage")

#single linakge
dev.new()
sing.hclust = hclust(Dmat,method="single")
plot(sing.hclust,cex=.7,main="Single Linkage")

#average linakge
dev.new()
ave.hclust = hclust(Dmat,method="average")
plot(ave.hclust,cex=.7,main="Average Linkage")

#Ward's linakge
dev.new()
ward.hclust = hclust(Dmat,method="ward.D")
plot(ward.hclust,cex=.7,main="Ward's Linkage")

#complete linkage with different distances. one can also use 1- cor(X) as a distance measure!
## commonly used in the clustering of gene expression.

dev.new()
Dmat = dist(t(ncidat),method="manhattan") #L1 distance
com.hclust = hclust(Dmat,method="complete")
plot(com.hclust,cex=.7,main="Complete Linkage - L1 Dist")


##########
#Biclustering - Cluster Heatmap

require("ISLR")
ncidat = t(NCI60$data)
colnames(ncidat) = NCI60$labs

#filter genes using PCA
X = t(scale(t(ncidat),center=TRUE,scale=FALSE))
sv = svd(t(X));
V = sv$v

#PC loadings - visualize data by limiting to top genes in magnitude in the PC loadings
aa = grep("grey",colors())
bb = grep("green",colors())
cc = grep("red",colors())
gcol2 = colors()[c(aa[1:30],bb[1:20],rep(cc,2))]

j = 2
ord = order(abs(V[,j]),decreasing=TRUE)
x = as.matrix(X[ord[1:250],])

#cluster heatmap - uses Ward's linkage (complete is default)
heatmap(x,col=gcol2,hclustfun=function(x)hclust(x,method="ward.D"))


### use heatmap.2, heatmap.3, gapmap, d3heatmap!!!!




######################################################
