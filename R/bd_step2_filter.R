options(stringsAsFactors=F)
library(gdata)
library(parallel)
files = list.files(path='ctx/',pattern='*.bd$')
meta = read.csv("WGS.coverage.csv")

mclapply (files, function(f) {
	dat = read.delim(sprintf('ctx/%s', f),comment.char='#',header=F,as.is=T)[,-(12:14)]
	message(sprintf("File: %s, Dim: (%s)", f, paste(dim(dat), collapse=",")))
	
	## case
	case = gsub("bd", "TCGA", substr(f,0,10))
	dat = data.frame(Case.ID=case, dat)
	
	## tumor reads threshold
	idx = match(paste(case,"01",sep="-"),meta$Sample.ID)
	tt = floor(meta$Coverage[idx] * 0.25)
	tt = ifelse(tt < 1, 1, tt)
	
	## normal reads threshold
	idx = match(paste(case,"10",sep="-"), meta$Sample.ID)
	nt = ceiling(meta$Coverage[idx] * 0.05)
	nt = ifelse(nt < 1, 1, nt)
	
	## read pairs supporting
	rp = strsplit(dat[,12], ":")
	rp = lapply(rp, strsplit, "\\|")
	rp = lapply(rp, unlist)
	imax=length(rp)

	regex_t = "TCGA-[\\w]{2}-[\\w]{4}-01"
	regex_n = "TCGA-[\\w]{2}-[\\w]{4}-10"

	idx = lapply(rp, function(x) grep(regex_n, x, perl=T))
	idx = sapply(idx, function(x) ifelse(length(x)==0,NA,x))
	dat$NormalReads = as.numeric(sapply(1:imax, function(i) rp[[i]][idx[i]+1]))

	idx = lapply(rp, function(x) grep(regex_t, x, perl=T))
	idx = sapply(idx, function(x) ifelse(length(x)==0,NA,x))
	dat$TumorReads = as.numeric(sapply(1:imax, function(i) rp[[i]][idx[i]+1]))

	dat$PropNormal = round(dat$NormalReads / (dat$NormalReads + dat$TumorReads),2)
	
	dat = subset(dat, !is.na(dat$TumorReads) & dat$TumorReads >= tt)
	dat = subset(dat, is.na(dat$NormalReads) | dat$NormalReads <= nt)
	dat[,12] = NULL

	colnames(dat)=c('Case.ID','ChrA','PosA','OrientA','ChrB','PosB','OrientB','Type','Size','Score','reads','n_reads', 't_reads', 'prop_n')

	fout = sprintf("filtered/%s.filtered.bd", case)
	write.table(dat, file=fout, quote=F, sep="\t", row.names=F, col.names=T)
}, mc.cores=24)
