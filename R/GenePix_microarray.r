### limma for GenePix microarray data
## read https://www.bioconductor.org/packages/3.3/bioc/vignettes/limma/inst/doc/usersguide.pdf
## http://biocourse.wp.sanbi.ac.za/wp-content/uploads/sites/7/2013/01/day3.pdf


## read here on how to install bioconductor packages https://www.bioconductor.org/install/
source("https://bioconductor.org/biocLite.R")
biocLite("limma")

## load the library 
library(limma)

## set the folder containing the gpr files as working directory 
setwd("/Users/mtang1/Downloads/1H315112001/Gpr")

## according to the service report, it is a single channel experiment with 5ug Cy5-labeled aRNA
targets <- readTargets("targets.txt")

# Set up a filter so that any spot with a flag of −99 or less gets zero weight. 

f <- function(x) as.numeric(x$Flags > -99)

# read in the data

## this is single channel data using only Cy5 channel according to the service report , one has to fool limma by:
## read Gorden Smith's reply here https://support.bioconductor.org/p/22328/
## http://permalink.gmane.org/gmane.science.biology.informatics.conductor/51460

## it is a bit complicated for limma to read in single channel 

## Cy3 green channel 
## Cy5 red channel 

Cy3 <- "F532 Median"
Cy5<- "F635 Median"
RG <- read.maimages(targets$FileName, source="genepix", columns = list(R=Cy5, G=Cy5), wt.fun=f)

RG.b<- backgroundCorrect(RG=RG, method='normexp')


## get rid of the fake channel
rownames(RG.b$R) <- RG.b$genes$Name
RG.b$G <- NULL


## design matrix
design <- model.matrix(~ 0+factor(c(1,1,1,2,2,2)))
colnames(design) <- c("group1", "group2")
contrast.matrix <- makeContrasts(group2-group1, levels=design)

library(vsn)
norm.vsn <- normalizeVSN(RG.b$R)
fit.vsn <- lmFit(norm.vsn, design)
fit2.vsn <- contrasts.fit(fit.vsn, contrast.matrix)
fit3.vsn <- eBayes(fit2.vsn)

limma::volcanoplot(fit3.vsn)

toptable<- topTable(fit3.vsn, number= Inf, sort.by = 'P')

hist(toptable$P.Value , breaks=100, col="skyblue", border="slateblue", main="")
hist(toptable$adj.P.Val, breaks=100, col="skyblue", border="slateblue", main="")

toptable %>% arrange(adj.P.Val) %>% head
# Make a basic volcano plot
with(toptable, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot", xlim=c(-6,6)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(toptable, P.Value<.05 ), points(logFC, -log10(P.Value), pch=20, col="red"))
with(subset(toptable, abs(logFC)>1.5), points(logFC, -log10(P.Value), pch=20, col="orange"))
with(subset(toptable, P.Value <.01 & abs(logFC)>1.5), points(logFC, -log10(P.Value), pch=20, col="green"))
abline(h=2)
abline(v=1.5)
abline(v=-1.5)
# Label points with the textxy function from the calibrate plot

library(calibrate)
with(subset(toptable, adj.P.Val<.5 & abs(logFC)> 2), textxy(logFC, -log10(P.Value), labs=ID, cex=.8))

write.table(toptable, "differential_genes.txt", sep="\t", row.names= F, quote = F)


###### PCA analysis and heatmap

colnames(norm.vsn)<- c("C1", "C2", "C3", "PO1", "PO2", "PO3")
plotPCA_SVD(norm.vsn)


subset_by_PC_loadings<- function(x, n){
        X<- t(scale(t(x),center=TRUE,scale=FALSE))
        sv<- svd(t(X))
        U<- sv$u
        V<- sv$v
        D<- sv$d
        k=1
        ord1<- order(abs(V[,k]),decreasing=TRUE)
        x1<- as.matrix(X[ord1[1:n],])
        x1
}

x1<- subset_by_PC_loadings(norm.vsn, 1000)

library(pheatmap)
pheatmap(x1, show_rownames = F, clustering_method = "ward.D")


library(RColorBrewer)
library(gplots)
hmcols<- colorRampPalette(brewer.pal(9,"GnBu"))(100)

colnames(x1)<- c("C1", "C2", "C3", "PO1", "PO2", "PO3")
heatmap.2(x1, distfun=function(x) as.dist(1-cor(t(x))), 
          hclustfun=function(x)hclust(x,method="ward.D"),trace="none", 
          scale = "row", col=hmcols, labCol=colnames(x1), labRow = F, margins = c(10,6), 
          density.info="none", main = "top 1000 regions")
