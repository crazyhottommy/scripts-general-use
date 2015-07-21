#author data

load("UnsupL.Rdata")

#understand the data a bit
dim(author)
colnames(author)
unique(rownames(author))
TrueAuth = as.factor(rownames(author))


X<- author
S<- svd(author)
U<- S$u
D<- S$d
V<- S$v
Z<- X %*% V

cols = as.numeric(as.factor(rownames(author)))
plot(Z[,1], Z[,2], type ="n")
text(Z[,1], Z[,2], rownames(X), col=cols)

#### using MDS


d<- dist(author)
d<- dist(author, method = "manhattan")

## it looks like manhattan distance is better than eucledian distance to separate the authors

mds<- cmdscale(d, k=2)
plot(mds[,1], mds[,2], type="n")
text(mds[,1], mds[,2], rownames(mds), col=cols)


#### k-means clustering
set.seed(3)
K = 4
km = kmeans(author,centers=K)

plot(Z[,1],Z[,2],col=km$cluster,type="n")
text(Z[,1],Z[,2],rownames(author),cex=.75,col=km$cluster)
cens = km$centers
points(cens%*%V[,1],cens%*%V[,2],col=1:K,pch=16,cex=3)

##### hierachical clustering
library(rafalib)

Dmat = dist(author) ## using default eucledian distance

com.hclust = hclust(Dmat,method="complete")
plot(com.hclust,cex=.7,main="Complete Linkage", labels= rownames(author))
myplclust(com.hclust, labels=rownames(author), lab.col=as.fumeric(rownames(author)))

Dmat<- dist(author, method= "canberra")
#average linakge
dev.new()
ave.hclust = hclust(Dmat, method="average")
myplclust(ave.hclust, labels=rownames(author), lab.col=as.fumeric(rownames(author)))

#Ward's linakge
dev.new()
ward.hclust = hclust(Dmat, method="ward.D")
myplclust(ward.hclust, labels=rownames(author), lab.col=as.fumeric(rownames(author)))


