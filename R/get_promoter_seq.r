
## get all the promoter sequences for human hg19 genome
## Author: Ming Tang (Tommy)
## Date: 04/30/2015

## load the libraries 
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
BSgenome.Hsapiens.UCSC.hg19
# or
Hsapiens
class(Hsapiens)

##BSgenome contains the DNA sequences

#############
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
Txdb<- TxDb.Hsapiens.UCSC.hg19.knownGene
Txdb
class(Txdb)
## Genomic features, contains the gene structure information

## all the genes, returns a GRange object 
hg19_genes<- genes(Txdb)
hg19_genes

## you can also get all the exons by gene, returns a GRangelist object
exbg<- exonsBy(Txdb, by="gene")
exbg
############### from another database
library(Homo.sapiens)
ghs<- genes(Homo.sapiens)
ghs
#############  hg19_genes and ghs are the same, the ENTREZID is the geneID

## get all the transcription start sites (TSSs)
tssgr<- resize( ghs, 1)

## how many bases you want to extract upstream and downstream of every TSSs
## in this case, I want the 1500 bases around TSS
upstream_pad<- 1000 
downstream_pad<- 500

promoters<- promoters(tssgr, upstream_pad, downstream_pad)
head(width(promoters))

## in fact, promoters(ghs, upstream_pad, downstream_pad) gives the same result
## because the default site for calcuation is the start(ghs), which is the TSS

### get sequences for all the promoters (1500bp sequences in this case)
## it takes seconds
promoter_seq<- getSeq(Hsapiens, promoters)

## we can annotate the promoter sequences with gene names, symbols in addition to the entrez id.

annotation<- select(Homo.sapiens, key=names(promoter_seq), keytype="ENTREZID", 
       columns=c("SYMBOL", "ENSEMBL", "ENTREZID"))

dim(annotation)

## more rows than the promoter_seq
## Be aware that there are some duplicated rows for the same ENTREZID, because there are multiple ENSEMBL id for the 
## same ENTREZID

## Let's only use the SYMBOL instead
annotation<- select(Homo.sapiens, key=names(promoter_seq), keytype="ENTREZID", 
                    columns=c("SYMBOL", "ENTREZID"))

## Now, it is the one to one mapping
## write to a fasta file
writeXStringSet(promoter_seq, "promoters_1500.fa", append=FALSE,
                compress=FALSE,compression_level=NA, format="fasta")

## I really do not want to dig in to find a way to write the name using the SYMBOL rather than the ENTREZID....
##

## I will prefer to prepare a bed file for all the promoters using bedtools slop (RefSeq table from UCSC, or from a GENECODE GTF file)
## then, use bedtools to extract DNA sequences using bedtools getfasta.
## To me, it is more flexible on the command lines.
## my previous post here http://crazyhottommy.blogspot.com/2015/02/fetch-genomic-sequences-from-coordinates.html
