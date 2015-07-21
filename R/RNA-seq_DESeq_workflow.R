setwd("/Users/Tammy/Downloads")
library("DESeq")
countsTable<- read.delim("qiu_lab_with_header.txt", header=TRUE)
rownames(countsTable)<- countsTable$Gene
countsTable<- countsTable[,-1]
countsTable<- countsTable[which(rowSums(countsTable)>25),]
head(countsTable)
conds<- factor(c("GIE","GIE", "GATA1","V205M", "V205M","2RA"))

abDesign<- data.frame(row.names=colnames(countsTable),
                      condition = c("GIE","GIE","GATA1","V205M","V205M","2RA"),
                      libType = c("pair-end","pair-end","pair-end","pair-end","pair-end","pair-end"))


cds<- newCountDataSet(countsTable, abDesign$condition)
cds<- estimateSizeFactors(cds)
sizeFactors(cds)
head(counts(cds))
head(counts(cds,normalized=TRUE))
cds<- estimateDispersions(cds)

plotDispEsts(cds)

res <- nbinomTest( cds, 'GIE', 'GATA1' )
head(res)
plotMA(res)
# if I use the command above, the red colored points (padj<0.1) are not at the edge of the cloud
# rather they spread around.

# there is a bug in the plotMA function that messes up with the order of the dataframe.
# see a post here https://www.biostars.org/p/103855/#106636
# my question post https://www.biostars.org/p/106588/#106627
# to fix it now:
plotMA(res[order(res$padj),]) 

# mannually creat the MA plot
# x=subset(res, res$baseMean!=0)
# col=ifelse(x$padj>=0.1, "gray32", "red3")
# plot(x$baseMean, x$log2FoldChange, col=col, log="x", ylim=c(-5,5), pch=19, cex=0.5)
# abline(h=0, col='red')

hist(res$padj, breaks=100, col='skyblue', border='slateblue', main='')

resSig <- res[ res$padj < 0.01 & (res$log2FoldChange >1| res$log2FoldChange < -1), ]
resSig <- na.omit(resSig)
head(resSig)
write.table(resSig, "2RA_VS_V205M_2fold.txt", sep="\t", quote=F)

##### GAGE pathway analysis ###########
require(gage)
data(kegg.gs)
#use all the gene fold change 
deseq.fc<- res$log2FoldChange
names(deseq.fc)<- res$id
sum(is.infinite(deseq.fc))  # there are some infinite numbers, if use DESeq2, no such problem.
deseq.fc[deseq.fc>10]=10
deseq.fc[deseq.fc<-10]=-10
exp.fc<- deseq.fc

#kegg.gsets works with 3000 KEGG speicies
data(korg)
head(korg[,1:3], n=20)


#let's get the annotation files for mouse and convert the gene set to gene symbol format
kg.mouse<- kegg.gsets("mouse")
kegg.gs<- kg.mouse$kg.sets[kg.mouse$sigmet.idx]
lapply(kegg.gs[1:3],head)

# egSymb is only for human data, so eg2sym and sym2eg functions are only for human data.
#data(egSymb)
#kegg.gs.sym<- lapply(kegg.gs, eg2sym)
#lapply(kegg.gs.sym[1:3],head)

# to convert IDs among gene/transcript ID to Entrez GeneID or reverse, use eg2id and id2eg in the pathview package written by the same person.
library(pathview)
data(bods)
bods

gene.symbol.eg<- id2eg(ids=names(exp.fc), category='SYMBOL', org='Mm') # convert the gene symbol to Entrez Gene ID
head(gene.symbol.eg, n=100)
head(gene.symbol.eg[,2], n=10)

names(exp.fc)<- gene.symbol.eg[,2]

fc.kegg.p<- gage(exp.fc, gsets= kegg.gs, ref=NULL, samp=NULL)
sel<- fc.kegg.p$greater[,"q.val"] < 0.1 & !is.na(fc.kegg.p$greater[,"q.val"])
table(sel)

sel.l<- fc.kegg.p$less[,"q.val"] < 0.1 & !is.na(fc.kegg.p$greater[,"q.val"])
table(sel.l)

##### PCA analysis ###################
=======
        res = nbinomTest( cds, 'GATA1', 'GIE' )
head(res)
plotMA(res)

hist(res$padj, breaks=100, col='skyblue', border='slateblue', main='')

resSig = res[ res$padj < 0.01 & (res$log2FoldChange >1| res$log2FoldChange < -1), ]
head(resSig)
write.table(resSig, "2RA_VS_V205M_2fold.txt", sep="\t", quote=F)


cdsBlind<- estimateDispersions(cds, method= "blind")
vsd<- varianceStabilizingTransformation(cdsBlind)
plotPCA(vsd, intgroup=c("condition"))

######## heatmap #####################


library("RColorBrewer")
library("gplots")
library("genefilter")


##library size normalized counts heatmap
new_counts_table<- counts(cds, normalize=T)
new_counts_table<- log2(new_counts_table)
rv<- rowVars(new_counts_table)
idx<- order(-rv)[1:500]
hmcols<- colorRampPalette(c("green","green4","red","red4","yellow"))(256)
heatmap.2(new_counts_table[idx,],col=hmcols,trace="none",Colv=T,density.info="none")

###variance stablized counts heatmap
rv<-rowVars(exprs(vsd))  #vsd is an expression set object after variance stablizing transformation
idx<- order(-rv)[1:500]

hmcols<- colorRampPalette(c("green","red"))(256)

hmcols <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

# another color scheme 
hmcols<- rev(brewer.pal(11,name="RdBu"))


mat<- exprs(vsd)
mat<- mat[idx,]
mat<- mat[,c(1,3,4,6)]
hc_row<- hclust(as.dist(1-cor(t(mat))), method="complete")
hc_column<- hclust(as.dist(1-cor(mat)), method="complete")
heatmap.2(mat, col=hmcols, Rowv=as.dendrogram(hc_row), Colv=as.dendrogram(hc_column),trace="none", scale='row', density.info="none", labRow=NA)



heatmap.2(exprs(vsd)[idx,], col=hmcols, hclust=function(x) hclust(x, method='complete'), distfun=function(x) as.dist(1-cor(t(x))), trace="none",margin=c(10,6), scale='row', density.info="none",labRow=NA)
# ?cor:If x and y are matrices then the covariances (or correlations) between the columns of x and the columns of y are computed.
# ?dist: This function computes and returns the distance matrix computed by using the specified distance measure to compute the distances between the rows of a data matrix.
# so dist computes the distance between rows while cor computes the correlation between columns


### use the genes that are differentially expressed between GATA1 vs GIE
mat<- exprs(vsd)
idx<- rownames(mat) %in% resSig$id
mat<- mat[idx,]

mat<- mat[,c(1,3,4,6)]
hc_row<- hclust(as.dist(1-cor(t(mat))), method="complete")
hc_column<- hclust(as.dist(1-cor(mat)), method="complete")
hmcols<- colorRampPalette(c("green","red"))(256)

heatmap.2(mat, col=hmcols, Rowv=as.dendrogram(hc_row), Colv=as.dendrogram(hc_column),trace="none", scale='row', density.info="none", margin=c(10,6), labRow=NA)


#####sample distance  heatmap,  note the transposed the matrix #########
dists<- dist(t(exprs(vsd)))
mat<- as.matrix(dists)
rownames(mat)=colnames(mat)= with(pData(cdsBlind), paste(condition))
heatmap.2(mat, trace="none", col=rev(hmcols), margin=c(13,13))

#multi-dimension scaling  MDS plot
mds<- cmdscale(dists)
km<- kmeans(t(exprs(vsd)), centers=4)
plot(mds, col=km$cluster, pch=16)


##### make a vacano plot ##################

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pval), pch=20, main="Volcano plot", xlim=c(-6,6)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pval), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pval), pch=20, col="orange"))
with(subset(res, padj <.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pval), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot

library(calibrate)
with(subset(res, padj<.05 & abs(log2FoldChange)> 4), textxy(log2FoldChange, -log10(pval), labs=id, cex=.8))


#### heatmap for log2fold change ######
res1 <- nbinomTest( cds, 'GIE', 'GATA1' )
res2 <- nbinomTest( cds, 'GIE', '2RA' )
res3 <- nbinomTest( cds, 'GIE', 'V205M' )

df<- data.frame(id=res1$id, logFC1=res1$log2FoldChange,padj1=res1$padj,
                logFC2=res2$log2FoldChange, padj2=res2$padj,logFC3=res3$log2FoldChange, padj3=res3$padj)


df<- df[abs(df$logFC1)>1 & abs(df$logFC2)>1,]

df <- na.omit(df)

View(df)

m<- as.matrix(df[,c(2,4,6)])
m<- m[!rowSums(!is.finite(m)),]
colnames(m)<- c("GATA1_VS_G1E","2RA_VS_G1E","V205M_VS_G1E")
head(m)

# try different color schemes http://seqanswers.com/forums/showthread.php?t=12022
library(RColorBrewer)
hmcols<- colorRampPalette(c("green","green4","red","red4"))(256)
hmcols<- colorRampPalette(c("white","red","yellow"))(256)
hmcols<- colorRampPalette(c("white","green","green4","violet","purple"))(100)
hmcols<-colorRampPalette(c("red","green","blue"))(256)

# set scale='row'  get a standardized z-score to remove the mean gene expression 
# https://www.biostars.org/p/15285/
heatmap.2(m, col=hmcols,trace="none",Colv=T,density.info="none", labRow=NA, scale="row", margin=c(17,15))

# the default is scaled with row which means for each gene, it scaled in this three different groups.
# One subtle point in the previous examples is that the heatmap function has automatically scaled the 
#colours for each row (i.e. each gene has been individually normalised across patients). This can be 
# disabled using scale="none", which you might want to do if you have already done your own normalisation
# (or this may not be appropriate for your data) http://www2.warwick.ac.uk/fac/sci/moac/people/students/peter_cock/r/heatmap/

# the defaut distance calculating method in hclust is Euclidian, but it is problemetic. genes are upregulated
# did not cluster together. see a post here https://www.biostars.org/p/91978/ and 
# here http://liorpachter.wordpress.com/2014/01/19/why-do-you-look-at-the-speck-in-your-sisters-quilt-plot-and-pay-no-attention-to-the-plank-in-your-own-heat-map/
# let's try dist but using the method single
heatmap.2(m, col=hmcols, scale="row", hclust=function(x) hclust(x, method='single'), distfun=function(x) as.dist(1-cor(t(x))), trace="none",margin=c(17,15), density.info="none", labRow=NA)

# usually the best way to cluster using person correlation distance https://www.biostars.org/p/14156/
# let's try to compute the pearson distances instead
# do something like this http://stackoverflow.com/questions/6719747/heatmap-of-microarray-data-using-pearson-distance
# both row and colum are clustered by pearson distance
heatmap.2(m, col=hmcols, scale="row", hclust=function(x) hclust(x, method='complete'), distfun=function(x) as.dist(1-cor(t(x))), trace="none",margin=c(17,15), density.info="none", labRow=NA)

hm<- heatmap.2(m, col=hmcols, hclust=function(x) hclust(x, method='complete'), distfun=function(x) as.dist(1-cor(t(x))), trace="none",margin=c(17,15), scale='row', density.info="none", labRow=NA)

# to get the the matrix after clustering

names(hm)
# return the maxtrix returned after clustering as in the heatmap
m.afterclust<- m[rev(hm$rowInd),rev(hm$colInd)]

# to extract subgroups that are clustered together
# rowDendrogram is a list object 
# convert the rowDendrogram to a hclust object
hc<- as.hclust(hm$rowDendrogram)

names(hc)
plot(hc)  # rotate the dendrogram 90 degree, it is the same as in the heatmap
rect.hclust(hc,h=0.5) # based on the height of the tree, you can specify h

ct<- cutree(hc,h=0.5)

# get the members of each subgroup in the order of the cluster(left--->right), the row order will
# it is reversed compared to the heatmap.
table(ct)
ct[hc$order]


# get the matrix after clustering in the order of the heatmap (up--->down)

tableclustn<-  data.frame(m.afterclust, rev(ct[hc$order]))
head(tableclustn)
write.table(tableclustn, file="tableclustn.xls", row.names=T, sep="\t")

# remake the heatmap adding the RowSide bar based on the subgroups

png("Rheatmap2.png", width=400, height=800)
mycolhc<- sample(rainbow(256))
mycolhc<-mycolhc[as.vector(ct)]
rowDend<- as.dendrogram(hc)

heatmap.2(m, col=hmcols, RowSideColors=mycolhc,hclust=function(x) hclust(x, method='complete'), distfun=function(x) as.dist(1-cor(t(x))), trace="none",margin=c(20,15), scale='row', density.info="none", labRow=NA)

dev.off()