source("http://bioconductor.org/biocLite.R")

### DNA sequences for the genome
biocLite("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
BSgenome.Hsapiens.UCSC.hg19

# We can access chromosome 11 like this:
chr11seq <- BSgenome.Hsapiens.UCSC.hg19[["chr11"]]
# Here, for example, is a segment of 25 bases starting  at base 1 million 
subseq(chr11seq,start=10^6,width=25)

## how many ATG in the chr11?
countPattern("ATG", chr11seq)
countPattern("TGA", chr11seq)
countPattern("TAA", chr11seq)
countPattern("TAG", chr11seq)

chr7seq <- BSgenome.Hsapiens.UCSC.hg19[["chr7"]]
alphabetFrequency(chr7seq, as.prob=T)

```{r}

```

### dbSNP database
biocLite("SNPlocs.Hsapiens.dbSNP.20120608")
library(SNPlocs.Hsapiens.dbSNP.20120608)
s17<- getSNPlocs("ch17")
head(s17)
class(s17)
library(dplyr)

#What is the location on chr17 of SNP rs73971683
s17 %>% filter(RefSNP_id=="73971683")

### gene annotations for the human genome 
biocLite("Homo.sapiens")
library(Homo.sapiens)
class(Homo.sapiens)
# possible keys 
keytypes(Homo.sapiens)

#There are also columns in the database, not all of which are keys. To list all the columns:
columns(Homo.sapiens)

head(keys(Homo.sapiens, keytype="ENTREZID"))

# how many unique ERNTREZID?
length(unique(keys(Homo.sapiens, keytype="ENTREZID")))

# how many ensemble id?
length(unique(keys(Homo.sapiens, keytype="ENSEMBL")))

# use select command to look up

select(Homo.sapiens, key="123", keytype="ENTREZID", columns=c("SYMBOL", "ENSEMBL", "ENTREZID", "CHR"))
select(Homo.sapiens, key="9575", keytype="ENTREZID", columns=c("SYMBOL", "ENSEMBL", "ENTREZID", "CHR"))

##We can list genes associated with certain biological processes or molecular functions if we know 
#the right vocabulary. One such vocabulary is the "Gene Ontology", often shortened "GO". 
#One term of interest in GO is "circadian rhythm". We use the following command to enumerate 
#genes by Entrez ID annotated to this process:

tab<- select(Homo.sapiens, key="circadian rhythm", keytype="TERM", columns=c("ENTREZID"))


### gene expression
library(devtools)
install_github("genomicsclass/tissuesGeneExpression")
library(tissuesGeneExpression)
data(tissuesGeneExpression)
head(e[,1:5])
table(tissue)

d<- dist(t(e))
hc<- hclust(d)
# the column names of e match the tissue names
plot(hc, cex=0.5, label=tissue)

abline(h=120)
cl<- cutree(hc, h=120)
table(true=tissue, cluster=cl)

as.numeric(e[rownames(e)== "209169_at",])
probe<- data.frame(value=as.numeric(e[rownames(e)== "209169_at",]), tissue=tissue)
ggplot(probe) + geom_boxplot(aes(x=tissue,y=value))

IDs<- c("201884_at", "209169_at", "206269_at", "207437_at", "219832_s_at", "212827_at")

probes<- as.data.frame(t(e[rownames(e) %in% IDs,]))
probes<- mutate(probes, tissue=tissue)
colnames(probes)<- c("A","B","C","D","E","F", "tissue")
par(mfrow=c(1,1))
ggplot(probes) + geom_boxplot(aes(x=tissue, y=B))

########## Another gene expression
library(Biobase)
data(sample.ExpressionSet)
sample.ExpressionSet
samp<- sample.ExpressionSet

# Note that you can access the information about the samples with:
pData(samp)

colnames(exprs(samp))
paste(colnames(exprs(samp)), rownames(pData(samp)), sep=".")
# match each other
head(exprs(samp))

Females<- rownames(pData(samp)[pData(samp)$sex=="Female",])
# expression values for only females 
head(exprs(samp)[, colnames(exprs(samp)) %in% Females])

experimentData(samp)
annotation(samp) 

cor(exprs(samp)[rownames(exprs(samp)) == "31489_at",], pData(samp)$score)

