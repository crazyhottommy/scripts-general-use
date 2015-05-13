library(affy)

basedir<- "/Users/Tammy/rawdata-master/"
file.path(basedir, "celfiles/sampleinfo.txt")

sample<- read.delim(file.path(basedir, "celfiles/sampleinfo.txt"), header=TRUE)

sample[sample$filenames=="1521a99hpp_av06.CEL.gz", "X36311_at"]

basedir<- "/Users/Tammy/rawdata-master/celfiles"
fns<- list.celfiles(basedir)
fns %in% sample[,1]
ab<- ReadAffy(filenames= file.path(basedir, sample[,1]), phenoData=sample)

table(probeNames(ab)=="36311_at")
## there are 16 probes for the same feature 36311_at

## gene expression for all the probes
pm(ab)[1:10,1:6]

## note it is different from exprs(ab) which summarise the gene expression at the feature level

which(probeNames(ab)=="36085_at")

## extract all the probes for this feature 36085_at
par(mfrow=c(1,1))
mat<- pm(ab)[which(probeNames(ab)=="36085_at"),]
matplot(mat, type="l")
matplot(log2(mat[,"1532a99hpp_av04.CEL.gz"]), type="l")
matplot(log2(mat[,"1532b99hpp_av04.CEL.gz"]), type="l")
matplot(mat[,"1532a99hpp_av04.CEL.gz"], mat[,"1532b99hpp_av04.CEL.gz"], type="l")

eset<- rma(ab)
g<-  factor(pData(eset)[,2])

library(genefilter)
rowttests(eset, g)["36085_at",]

### spike-ins
sig<- colnames(pData(ab))[-1]

##get rid of the leading X for each feature
sig

rowttests(eset,g)[sig,]
library(limma)
datadir = "/Users/Tammy/rawdata-master/"
basedir = file.path(datadir, "/agilent")
setwd(basedir)
targets = readTargets("TargetBeta7.txt")
RG = read.maimages(targets$FileName, source="genepix")
MA<- MA.RG(RG, bc.method="none")


###################### NGS data
biocLite("pasillaBamSubset")
library(pasillaBamSubset)
library(Rsamtools)

# the path to the bam file
filename <- untreated1_chr4()

# add parathensis will evaluate the expression
(bf <- BamFile(filename))

seqinfo(bf)
(sl <- seqlengths(bf))

# summary of bam file
quickBamFlagSummary(bf)

?scanBamParam

## genomic ranges for chr4
(gr <- GRanges("chr4",IRanges(1, sl["chr4"])))

countBam(bf, param=ScanBamParam(which = gr))

?ScanBamParam
ScanBamParam(which = gr)

reads<- scanBam(BamFile(filename, yieldSize=5))
?scanBam

class(reads)
names(reads[[1]])

reads[[1]]$pos
reads[[1]]$rname
reads[[1]]$strand
reads[[1]]$seq


library(GenomicAlignments)
(ga <- readGAlignments(bf))

length(ga)
granges(ga[1])


gr <- GRanges("chr4", IRanges(700000, 800000))
(fo <- findOverlaps(ga, gr)) # which reads over this range

countOverlaps(gr, ga) # count overlaps of range with the reads

table(ga %over% gr) # logical vector of read overlaps with the range

## If we had run countOverlaps(ga, gr) it would return an integer vector 
## with the number of overlaps for each read with the range in gr.
countOverlaps(ga, gr)


gr <- GRanges("chr4", IRanges(440000, 470000))
countOverlaps(gr, ga)

reads <- scanBam(bf, param=ScanBamParam(what=c("pos","strand","seq"), which=gr))
## access sequences of the reads by
reads$`chr4:440000-470000`$seq


## GC content of the sequences 
apply(letterFrequency(x=reads$`chr4:440000-470000`$seq, letters=c("C","G"), as.prob=T), 1, sum)


mean (apply(letterFrequency(x=reads$`chr4:440000-470000`$seq, letters=c("C","G"), as.prob=T), 1, sum))


biocLite("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
g<- genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene)


## Let's focus on the genes at a particular location:

g2<- g[g %over% GRanges("chr4",IRanges(200000, 300000))]

g2

g2['FBgn0039890']
ranges(g2['FBgn0039890'])

## how many reads mapped on to this gene?
countOverlaps(g2['FBgn0039890'], ga, ignore.strand=FALSE)

##### count reads from a bam file

library(pasillaBamSubset)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
library(GenomicAlignments)



## genes on chr4
g<- genes(txdb)
g<- g[seqnames(g) == "chr4"]

txdb<- TxDb.Dmelanogaster.UCSC.dm3.ensGene
grl<- exonsBy(txdb, by="gene")

#GRangesList for chr4
grl<- grl[names(g)]

# the 10th gene
grl[10]

# access the GRanges
grl[[10]]

#first exon
grl[[10]][1]

# Test for the same names:
all.equal(names(g), names(grl))

bf<- BamFile(untreated1_chr4())
so1<- summarizeOverlaps(features=grl,
                        reads= bf,
                        ignore.strand=TRUE)
head(assay(so1))
colSums(assay(so1))
rowData(so1)
colData(so1)
colData(so1) <- c("one", "two")
metadata(rowData(so1)) # rowData depricated
metadata(rowRanges(so1))

so2<- summarizeOverlaps(features=g,
                        reads= bf,
                        ignore.strand=TRUE)

head(assay(so2))

plot(log2(assay(so2)[,1]), log2(assay(so1)[,1]))
abline(a=0, b=1)

## some reads mapped to exon-exon junctions, we see more reads fro g compared with grl
mean( assay(so1)[assay(so2)[,1] >0,] / assay(so2)[assay(so2)[,1] >0,] )


head(assay(so1))
## FPM (fragments per million counts for each gene (exons))
head( assay(so1)[,1]/colSums(assay(so1)) * 10^6 )


### FPKM 

## reduce can directly work on GRangelist, 
## try reduce on Grange object as well.

## look at the first gene
reduce(grl[[1]])
width(reduce(grl[[1]]))
width(grl[[1]])

## compare with this 
width(reduce(grl))


## sum of exons for each gene
sapply(width(reduce(grl)), sum)
exons_length_per_gene<-sapply(width(reduce(grl)), sum)

## 
summary(sapply(width(reduce(grl)), sum))

### calculate RPKM
FPKM<- (assay(so1)[,1]/colSums(assay(so1)) * 10^6)/exons_length_per_gene * 1000
FPKM["FBgn0002521"]

#### visualize 
# coverage of the GenomicAlignment object ga
fl1<- untreated1_chr4()
fl2<- untreated3_chr4()

x<- BamFile(fl1)
y<- BamFile(fl2)

xcov<- coverage(x)
xcov
z<- GRanges("chr4", IRanges(456500,466000))
# subset the coverage using z
xcov[z]  # it is a Rlelist
## or
xcov$chr4[ranges(z)]  # it is integer-Rle

xnum<- as.numeric(unlist(xcov[z]))
plot(xnum)


ycov<- coverage(y)
ynum<- as.numeric(ycov$chr4[ranges(z)])
plot(xnum, type="l", col="blue", lwd=2)
lines(ynum, col="red", lwd=2)


## zoom in
plot(xnum, type="l", col="blue", lwd=2, xlim=c(6200,6600))
lines(ynum, col="red", lwd=2)


### check gene lgs using biomart
library(biomaRt)
listMarts()
m<- useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
lf<- listFilters(m)
lf[grep("name", lf$description, ignore.case=TRUE),]

map<- getBM(mart=m,
            attributes=c("ensembl_gene_id", "flybasename_gene"),
            filters="flybasename_gene",
            values="lgs")

map


grl<- exonsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene, by="gene")
gene<- grl[[map$ensembl_gene_id[1]]]
gene

### visualize by Gviz
library(Gviz)
gtrack<- GenomeAxisTrack()
atrack<- AnnotationTrack(gene, name="Gene Model")
plotTracks(list(gtrack, atrack))

xgr<- as(xcov, "GRanges")
ygr<- as(ycov, "GRanges")
dtrack1<- DataTrack(xgr[xgr %over% z], name="sample 1")
dtrack2<- DataTrack(ygr[ygr %over% z], name="sample 2")

plotTracks(list(gtrack, atrack, dtrack1, dtrack2))
plotTracks(list(gtrack, atrack, dtrack1, dtrack2), type="polygon")

## use ggbio  built on ggplot2
library(ggbio)
autoplot(gene)
autoplot(fl1, which=z)
autoplot(fl2, which=z)

