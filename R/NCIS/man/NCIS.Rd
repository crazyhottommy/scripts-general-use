\name{NCIS}
\alias{NCIS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
NCIS
}
\description{
Estimates normalizing factor between a ChIP sample and a control/input sample}
\usage{
NCIS(chip.data, input.data, data.type = c("MCS", "BED", "AlignedRead"), frag.len = 200, min.binsize = 100, max.binsize = 20000, binsize.shift = 100, min.stop.binsize = 100, chr.vec = NULL, chr.len.vec = NULL, quant=0.75)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chip.data}{
ChIP data.
}
  \item{input.data}{
control data.
}
  \item{data.type}{
"MCS", "AlignedRead" or "BED".
}
  \item{frag.len}{
average fragment length. Default 200 bp.
}
  \item{min.binsize}{
minimum of binsize to search.
}
  \item{max.binsize}{
maximum of binsize to search.
}
  \item{binsize.shift}{
the threshold of binsize after which the normalization factor is computed as the average of two estimates, one on regular bins and the other on bins shifed half binsize.
}
  \item{min.stop.binsize}{
minimum of binsize to use (stop).
}
  \item{chr.vec}{
vector of chromosomes in the data. Only reads in chr.vec are considered for normalization purpose.
}
  \item{chr.len.vec}{
vector of chromosome lengths corresponding to chr.vec
}
  \item{quant}{
quantile to start the searching for total threshold. Default 0.75.
}
}
\details{
data.type
\itemize{
\item \code{MCS} Minimum ChIP-Seq format. data.frame with fields: chr (factor), pos (integer) and strand (factor, "+" and "-"). pos is 5' location. This is different from eland default which use 3' location for reverse strand.
\item \code{AlignedRead} from Bioconductor ShortRead package (with support of commonly used formats, including Eland, MAQ, Bowtie, SOAP and BAM).
\item \code{BED} with at least first 6 fields (chrom, start, end, name, score and strand), http://genome.ucsc.edu/FAQ/FAQformat.html#format1.
}
}
\value{
A list with following components:
\item{est }{the estimated normalizing factor.}
\item{binsize.est }{the estimated binsize.}
\item{r.seq.depth}{sequencing depth ratio.}
\item{pi0}{the estimated proportion of background reads among ChIP sample.}
}
\references{
Liang, K and Keles, S (2012). Normalization of ChIP-seq data with control, to appear.
}
\author{
Kun Liang, \email{kliang@stat.wisc.edu}
}

\examples{
data("SEG3")
res <- NCIS(chip.data, input.data, data.type="MCS")
res
}
