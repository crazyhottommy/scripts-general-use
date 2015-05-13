
## The web link for this lesson
## https://courses.edx.org/courses/HarvardX/PH525.4x/1T2015/courseware/92dbe89dc80c4fc084cad3f00b0381f7/494a3c869b7d4d9facbc92c1c8ded6c6/

download.file("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData","bottomly_eset.RData")
load("bottomly_eset.RData")
library("Biobase")
pData(bottomly.eset)

# raw counts matrix
head(exprs(bottomly.eset))

total<- nrow(exprs(bottomly.eset))

ratios<- numeric()
for (i in 1:ncol(exprs(bottomly.eset))) {
        zeros<- sum(exprs(bottomly.eset)[,i] == 0)
        ratio<- zeros/total
        ratios<- c(ratios, ratio)
}

g<- as.factor(pData(bottomly.eset)$experiment.number)

library(ggplot2)
d<- data.frame(ratios=ratios, group=g)
ggplot(d) + geom_boxplot(aes(x=group, y=ratios))

## different batches have differnt proportion of 0s for the counts 

## make a multidimensional scaling plot
##  However, when computing distances we prefer
# to have measures with similar variability across samples. 
## Because this is RNAseq data, the variance depends on the mean:
library(genefilter)
A<- rowMeans(exprs(bottomly.eset))
SD<- rowSds(exprs(bottomly.eset))
plot(A,SD)

## log2 transformation to stablize the variance
y<- log2( exprs( bottomly.eset )+0.5)

library(devtools)
install_github("rafalib","ririzarr") 
library(rafalib)
mypar(2,1)
hist(y[,1],nc=100)
hist(y[y[,1]>0,1],nc=100)
abline(v=3)

## use counts bigger than 8
y<- exprs( bottomly.eset )
ind<- which(apply( y>=8, 1, all))
y<- log2( y[ind,] )

## make a multidimensional scaling plot

s<- svd(y)
U<- s$u  ## orthogonal matrix  UtU = I 
V<- s$v  ## orthogonal matrix 
D<- diag(s$d)  # diagnonal matrix

library(rafalib)
mypar2(1,1)
plot(s$d)

Yhat<- U %*% D %*% t(V)
resid<- Y - Yhat
boxplot(resid, ylim = c(-2,2), range=0)


##### normalizing NGS count data for biased GC content, this is just to demonstrate
### that GC bias can affect counts, for real data analysis, GC bias correction methods
### for example:
### EDAseq bioconductor http://www.ncbi.nlm.nih.gov/pubmed/22177264

# Systematic bias is important to keep in mind, especially in experiments with many batches. Comparisons across batches are complicated by technical biases, and it is essential that conditions of interest (treatment vs untreated) be balanced across batches, or else the experiment will uncover a mix of biological differences and systematic bias, with no hope to disentangle the two. This issue comes up over and over again in genomic studies, since the beginning of microarray data and continuing to the present day, with massive sequencing-based genomic data sets.

# The ReCount project provides summaries of a number of RNA-seq experiments as ExpressionSet objects, so that labs around the world do not have to perform alignment just to "re-count" the number of reads which fall in genes. If you use data from this project make sure to cite the authors! Download the ExpressionSet for the Bottomly experiment, here. Load this object into R:

library(Biobase)
load("bottomly_eset.RData")
bottomly.eset
head(exprs(bottomly.eset))

# Examine the histogram of sample-sample correlations of log counts:
        
hist(cor(log(exprs(bottomly.eset) + 1)))

#There are a number of reasons for such high correlations:
        
#Biological correlation: genes are expressed similarly across organisms
# Systematic bias: technical bias affects the same genes similarly (for example, if the bias is due to the sequence of the gene, then this could be common across experiments)
# The presence of many zeros, which anchor the cloud of log counts at (0,0)
#Remove some of the rows with many zeros, and plot the histogram of correlations again:
        
mat = exprs(bottomly.eset)
rs = rowSums(mat)
hist(cor(log(mat[ rs > 10, ] + 1)))

#Subset and rename the ExpressionSet to remove the genes which have 2 or fewer reads across all samples:
        
e = bottomly.eset[ rowSums(exprs(bottomly.eset)) > 2, ]
dim(e)

# We now have 12971 genes and 21 samples. The genes are annotated with ENSEMBL IDs, 
# stored as rownames (and in fData(e)). The quickest way for us to check the GC content
# of each gene is to use the precomputed TxDb for mouse in Bioconductor. If we were doing 
# this analysis not just for quick check, we should load the ENSEMBL genes from the ENSEMBL 
# GTF file and use these for computing GC.

biocLite("TxDb.Mmusculus.UCSC.mm9.knownGene")
biocLite("BSgenome.Mmusculus.UCSC.mm9")
biocLite("org.Mm.eg.db")

library("TxDb.Mmusculus.UCSC.mm9.knownGene")
library("BSgenome.Mmusculus.UCSC.mm9")
library("org.Mm.eg.db")


# Now we will connect the ENSEMBL IDs to Entrez IDs, and then get the sequence for the exons per gene.
## because in the TxDb.Mmusculus.UCSC.mm9.knownGene database, genes are in ENTREZID.

res = select(org.Mm.eg.db, keys=rownames(e), keytype="ENSEMBL", columns="ENTREZID")

#Note there are 1 to many mappings, which will we ignore for this quick analysis. Now add a column to the ExpressionSet, by getting the first match of the ENTREZ IDs to our ENSEMBL IDs.

fData(e)$ENTREZ = res$ENTREZID[ match(rownames(e), res$ENSEMBL) ]

#There are some NA's for those ENSEMBL IDs missing matches, which we will just remove from analysis:

sum(is.na(fData(e)$ENTREZ))
e = e[ !is.na(fData(e)$ENTREZ) , ]

#Now it's only two steps to get the GC content: get the genes, and get the sequence for the genes. First we get the exons for each gene:

txdb = TxDb.Mmusculus.UCSC.mm9.knownGene
grl = exonsBy(txdb, by="gene")

#We need to subset the ExpressionSet again to only those genes for which we have ranges:

e = e[ fData(e)$ENTREZ %in% names(grl), ]

#Now put the exon GRangesList in the order of the counts ExpressionSet and reduce the exons GRangesList. 
## The reduce() call is necessary so that we don't overcount any overlapping exon sequence.

reduced.exons = reduce(grl[ fData(e)$ENTREZ ])

# Get the sequence for the reduced exons using getSeq. Note that the returned object is a 
# DNAStringSetList. For each gene, we have a DNAStringSet of the reduced exon sequence. 
# By calling DNAStringSet(lapply(x, unlist)) on the DNAStringSetList, 'x', we get a DNAStringSet 
# for each gene with the joined reduced exon sequence. Note that this unlisting does take a few minutes to perform though.

#What is the GC-content (a ratio, answer between 0 and 1) for the first gene in reduced.exons? 
# Hint: letterFrequency

?getSeq
mm9_genome<- BSgenome.Mmusculus.UCSC.mm9

exonSeq<- getSeq(mm9_genome, reduced.exons)
exonSeq_unlist<- DNAStringSet(lapply(exonSeq, unlist))
#G/C content for the first gene

sum(letterFrequency(x=exonSeq_unlist[1], letters=c("C","G"), as.prob=T))

gc<- apply(letterFrequency(x=exonSeq_unlist, letters=c("C","G"), as.prob=T), 1, sum)


# Plot the log counts over the GC content for one sample:
        
plot( log(exprs(e)[,1]+1) ~ gc)

# It might be easier to see a dependence with boxplots:
        
boxplot( log(exprs(e)[,1]+1) ~ cut(gc, 20))

# We can calculate the median log count for each bin of GC content for a single sample:
        
sapply(split(log(exprs(e)[,1]+1), cut(gc, 20)), median)

# We can perform this calculation across all samples as well:
        
gc.depend = sapply(1:ncol(e), function(i) sapply(split(log(exprs(e)[,i]+1), cut(gc, 20)), median))

# Plot the GC dependence and color the lines by batch identifier:
        
plot(gc.depend[,1], type="n", ylim=c(0,6))
batch = factor(e$experiment.number)
for (i in 1:ncol(e)) lines(gc.depend[,i], col=batch[i])
legend("bottom", levels(batch), col=1:3, lty=1)




