# read GEO data sets from NCBI

library(Biobase)
library(GEOquery)
#Download GDS file, put it in the current diretory, and load it:
gds<- getGEO('GDS2758', destdir=".") #hypoxia induction for MCF7 cells

# or open an existing GDS file (even if it is compressed):
#gds<- getGEO( filename="GDS2758.soft.gz")

gds # GDS class object many informations there.

colnames(Table(gds))
Table(gds)[1:10,1:6]

#turn the GDS object into an expression set object(using log2)
eset<- GDS2eSet(gds, do.log2=TRUE)
eset
featureNames(eset)[1:10]
sampleNames(eset)
pData(eset)
phenoData(eset)
#anotation file, packages from Bioconductor
annotation(eset)

source("http://www.bioconductor.org/biocLite.R")
biocLite("hgu133a.db")
library("hgu133a.db")
library(annotate)
ID <- featureNames(eset)
Symbol<- getSYMBOL(ID, "hgu133a.db")
fData(eset)<- data.frame(ID, Symbol)


#differentially expressed genes by limma
library(limma)
design<- model.matrix(~0+factor(c(1,1,1,2,2,2,3,3,3)))
colnames(design)<- c("group1", "group2","group3")
fit<- lmFit(eset, design)
contrast.matrix<- makeContrasts(group2-group1,group3-group2,group3-group1, levels=design)
fit2<- contrasts.fit(fit, contrast.matrix)
fit2<- eBayes(fit2)

toptable<- topTable(fit2, coef=1,number=30000, sort.by="P")
head(toptable)
write.table(toptable, "microarray_hypoxia_MCF7_all.txt" ,quote=F, sep="\t")

hist(toptable$adj.P.Val, breaks=100, col='skyblue', border='slateblue', main='adjust P-value distribution for differentially expressed genes ', xlab="adjust p-value")
##### make a vacano plot ##################

# Make a basic volcano plot
with(toptable, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot", xlim=c(-4,6)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(toptable, adj.P.Val<.05 ), points(logFC, -log10(P.Value), pch=20, col="red"))
with(subset(toptable, abs(logFC)>1), points(logFC, -log10(P.Value), pch=20, col="orange"))
with(subset(toptable, adj.P.Val <.05 & abs(logFC)>1), points(logFC, -log10(P.Value), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot

library(calibrate)
with(subset(toptable, adj.P.Val<.05 & abs(logFC)> 4), textxy(logFC, -log10(P.Value), labs=Symbol, cex=.8))



######## make a heatmap ##################


library("gplots")
names(fit2)
head(fit2$p.value)
selected<- p.adjust(fit2$p.value[,1], method='BH') <0.05  # a logical vector to subset the eset
table(selected)

esetSel<- eset[selected,]
mat<- exprs(esetSel)[,1:6]

hmcols<- colorRampPalette(c("green","green4","red","red4"))(256)

heatmap.2(mat, col=hmcols,trace="none", margin=c(10,6), scale='row', density.info="none")

# use pearson correlation distance instead 
heatmap.2(mat, col=hmcols, hclust=function(x) hclust(x, method='complete'), distfun=function(x) as.dist(1-cor(t(x))), trace="none",margin=c(10,6), scale='row', density.info="none", labRow=NA)
