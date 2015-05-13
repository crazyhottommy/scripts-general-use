#NCIS (Normalization for ChIP-Seq) estimates normalizing factor between a ChIP sample and a control/input sample
#input parameters:
#chip.data	 ChIP data.
#input.data	 control data.
#data.type	 "MCS", "AlignedRead" or "BED".
#frag.len	 average fragment length. Default 200 bp.
#min.binsize	 minimum of binsize to search.
#max.binsize	 maximum of binsize to search.
#binsize.shift	 the threshold of binsize after which the normalization factor is computed as the average of two estimates, one on regular bins and the other on bins shifed half binsize.
#min.stop.binsize	 minimum of binsize to use (stop).
#chr.vec	 vector of chromosomes in the data. Only reads in chr.vec are considered for normalization purpose.
#chr.len.vec	 vector of chromosome lengths corresponding to chr.vec
#contact: kliang@stat.wisc.edu
#last modified: 2012.05.29

NCIS <- function(chip.data, input.data, data.type=c("MCS", "BED", "AlignedRead"), frag.len=200, min.binsize=100, max.binsize=20000, 
    binsize.shift=100, min.stop.binsize=100, chr.vec=NULL, chr.len.vec=NULL, quant=0.75){
    if(data.type=="MCS"){
        chip <- read.MCS(chip.data)
        input <- read.MCS(input.data)
    }else{
        if(data.type=="AlignedRead"){
            require(ShortRead)
            chip <- read.AlignedRead(chip.data)
            input <- read.AlignedRead(input.data)
        }else{
            if(data.type=="BED"){
                chip <- read.BED(chip.data)
                input <- read.BED(input.data)
            }else{
                stop("Unknown data format: type can only be 'NCIS', 'BED' or 'AlignedRead'")
            }
        }
    }
    shift.size <- round(frag.len/2)
    NCIS.internal(chip, input, shift.size=shift.size, min.binsize=min.binsize, max.binsize=max.binsize, 
    binsize.shift=binsize.shift, min.stop.binsize=min.stop.binsize, chr.vec=chr.vec, chr.len.vec=chr.len.vec, quant=quant)
}

#NCIS takes chip.pos and input.pos
NCIS.internal <- function(chip.pos, input.pos, shift.size=100, min.binsize=100, max.binsize=20000, 
    binsize.shift=100, min.stop.binsize=100, chr.vec=NULL, chr.len.vec=NULL, quant=0.75){

    if(is.null(chr.vec)){
        chip.name <- names(chip.pos)
        input.name <- names(input.pos)
        chr.vec <- intersect(chip.name, input.name)
        if(length(chr.vec)<length(chip.name)){
            cat("Control sample doesn't have chromosome", setdiff(chip.name, chr.vec), "which are in ChIP sample. These chromosomes are ignored.\n")
        }
        if(length(chr.vec)<length(input.name)){
            cat("ChIP sample doesn't have chromosome", setdiff(chip.name, chr.vec), "which are in control sample. These chromosomes are ignored.\n")
        }
    }
    
    nchip <- sum(sapply(chr.vec, function(x) sapply(chip.pos[[x]], length)))
    ninput <- sum(sapply(chr.vec, function(x) sapply(input.pos[[x]], length)))
    r.seq.depth <- nchip/ninput
    
    sizevec <- rep(c(1,2,5), times=3)*10^rep(2:4, each=3)
    sizevec <- sizevec[sizevec <= max.binsize]
    sizevec <- sizevec[sizevec >= min.binsize]
    norm.est <- rep(1000, length(sizevec))
    #names(norm.est) <- sizevec

    if(!is.null(chr.len.vec) & length(setdiff(chr.vec, names(chr.len.vec)))==0){
        chr.end.max <- chr.len.vec[chr.vec]
    }else{
        chr.end.max <- sapply(chr.vec, function(chr) max( c(max(chip.pos[[chr]][["+"]]), max(chip.pos[[chr]][["-"]]), 
                max(input.pos[[chr]][["+"]]), max(input.pos[[chr]][["-"]])) ))
        #names(chr.end.max) <- chr.vec
    }
    
    binsize.est <- -1
    
    for(si in 1:length(sizevec)){
        binsize <- sizevec[si]
        
        bindata <- bin.data(chip.pos, input.pos, binsize, shift.size=shift.size, chr.vec=chr.vec, chr.end.max=chr.end.max)
        
        res <- est.norm.quant.search(bindata$chip, bindata$input, quant=quant)
        
        if(binsize < binsize.shift){
            norm.est[si] <- res
        }else{
            #run it twice, 2nd time shift half binsize
            bindata <- bin.data(chip.pos, input.pos, binsize, shift.size=shift.size, shift.half.size=TRUE, chr.vec=chr.vec, chr.end.max=chr.end.max)
            
            res2 <- est.norm.quant.search(bindata$chip, bindata$input, quant=quant)

            norm.est[si] <- (res+res2)/2
        }#end else
        
        #stopping criteria
        if(si>1 & binsize.est<0){
            if(norm.est[si]>=norm.est[si-1]){
                est <- norm.est[si-1]
                binsize.est <- sizevec[si-1]
            }
            if(si==length(sizevec) & binsize.est<0){ #the end of binsize, no converge yet
                est <- norm.est[si]
                binsize.est <- sizevec[si]
            }
        } #end if(si>1 & binsize.est<0)
        if(binsize.est>0 & binsize>=min.stop.binsize){
            break
        }
    } #end for(si in 1:length(sizevec))

    return(list(est=est, binsize.est=binsize.est,
                r.seq.depth=r.seq.depth, pi0=est/r.seq.depth))
}


#data is dataframe with fields: chr(factor), pos(integer) and strand(factor, "+" and "-")
#pos is 5' location; this is different from eland default which use 3' location for reverse strand.
read.MCS <- function(data){
    if(is.data.frame(data)){
        #check data has required field
        if(!setequal(intersect(c("chr", "pos", "strand"), colnames(data)), c("chr", "pos", "strand"))) stop("MCS format need to be a data.frame with fields: chr, pos and strand")
        res <- split(data[, c("pos", "strand")], data$chr, drop = TRUE)
        res <- lapply(res, function(x) split(x[, "pos"], x$strand, drop = TRUE))
        return(res)
    }else{
        if(!setequal(names(data[[1]]), c("+", "-"))) stop("MCS list format need to be a list of chromosomes, each of which is a list of two strands.")
        return(data)
    }
}

#BED file should have at least first 6 fields (chrom, start, end, name, score and strand), see
#http://genome.ucsc.edu/FAQ/FAQformat.html#format1 for more details
#chip <- read.BED("BED.file")
read.BED <- function(bed.file){
    temp <- read.delim(bed.file, 
            header=FALSE, row.names=NULL)
    if(ncol(temp)<6) stop("BED file should has at least 6 fields!")
    colnames(temp)[1:6] <- c("chr", "pos", "posEnd", "t1", "t2", "strand")
    index <- (1:nrow(temp))[temp$strand=="-"]
    tt <- temp[index[1],]
    if(tt$pos < tt$posEnd){
        temp$pos[index] <- temp$posEnd[index]
    }
    res <- split(temp[, c("pos", "strand")], temp$chr, drop = TRUE)
    res <- lapply(res, function(x) split(x[, "pos"], x$strand, drop = TRUE))
    return(res)
}


#chip <- read.AlignedRead("AlignedRead.object")
read.AlignedRead <-
function(aln){
    if(is(aln, "AlignedRead")){
        res <- split(data.frame(pos=position(aln)+(strand(aln)=="-")*(width(aln)-1), strand=factor(strand(aln))), chromosome(aln), drop = TRUE)
        res <- res[names(res)!=""]
        res <- lapply(res, function(x) split(x[, "pos"], x$strand, drop = TRUE))
        return(res)
    }else{
        stop("Need an AlignedRead object.")
    }
}

bin.data <- function(chip.pos, input.pos, binsize,
    shift.size=100,
    shift.half.size=FALSE,
    zero.filter=TRUE,
    by.strand=FALSE,
    chr.vec=NULL, chr.end.max=NULL,
    by.chr=FALSE){

    if(is.null(chr.vec)){
        chip.name <- names(chip.pos)
        input.name <- names(input.pos)
        chr.vec <- intersect(chip.name, input.name)
        if(length(chr.vec)<length(chip.name)){
            cat("Control sample doesn't have chromosome", setdiff(chip.name, chr.vec), "which are in ChIP sample. These chromosomes are ignored.\n")
        }
        if(length(chr.vec)<length(input.name)){
            cat("ChIP sample doesn't have chromosome", setdiff(input.name, chr.vec), "which are in control sample. These chromosomes are ignored.\n")
        }
    }
    if(is.null(chr.end.max)){
        chr.end.max <- sapply(chr.vec, function(chr) max( c(max(chip.pos[[chr]][["+"]]), max(chip.pos[[chr]][["-"]]), 
                max(input.pos[[chr]][["+"]]), max(input.pos[[chr]][["-"]])) ))
    }
    chr.len <- ceiling((chr.end.max+shift.size)/binsize)

    chip.f <- list()
    chip.r <- list()
    input.f <- list()
    input.r <- list()
    if(!by.strand){
        chip <- list()
        input <- list()
    }

    for(chr in chr.vec){
        if(shift.half.size){
            bk <- c(0, seq(from=round(binsize/2), by=binsize, length.out=chr.len[[chr]]+1))
        }else{
            bk <- seq(from=0, by=binsize, length.out=chr.len[[chr]]+1)
        }
        bk.f <- bk-shift.size
        bk.r <- c(0, bk+shift.size)
        chip.f[[chr]] <- hist(chip.pos[[chr]][["+"]], breaks=bk.f, plot = FALSE)$counts
        chip.r[[chr]] <- hist(chip.pos[[chr]][["-"]], breaks=bk.r, plot = FALSE)$counts[-1]
        input.f[[chr]] <- hist(input.pos[[chr]][["+"]], breaks=bk.f, plot = FALSE)$counts
        input.r[[chr]] <- hist(input.pos[[chr]][["-"]], breaks=bk.r, plot = FALSE)$counts[-1]
        if(zero.filter){
            ind <- chip.f[[chr]]+chip.r[[chr]]+input.f[[chr]]+input.r[[chr]]>0
            chip.f[[chr]] <- chip.f[[chr]][ind]
            chip.r[[chr]] <- chip.r[[chr]][ind]
            input.f[[chr]] <- input.f[[chr]][ind]
            input.r[[chr]] <- input.r[[chr]][ind]
        }
        if(!by.strand){
            chip[[chr]] <- chip.f[[chr]]+chip.r[[chr]]
            input[[chr]] <- input.f[[chr]]+input.r[[chr]]
        }
    }
    if(by.strand){
        if(by.chr){
            return(list(chip.f=chip.f, chip.r=chip.r, input.f=input.f, input.r=input.r))
        }else{
            return(list(chip.f=unlist(chip.f, use.names = FALSE), chip.r=unlist(chip.r, use.names = FALSE), 
                input.f=unlist(input.f, use.names = FALSE), input.r=unlist(input.r, use.names = FALSE)))
        }
    }else{
        if(by.chr){
            return(list(chip=chip, input=input))
        }else{
            return(list(chip=unlist(chip, use.names = FALSE), input=unlist(input, use.names = FALSE)))
        }
    }    
}


est.norm.quant.search <- function(chip, input, quant=0.75){
   
    total <- chip+input
    tbl <- table(total)
    total.count <- as.integer(names(tbl))
    cum.bc <- cumsum(tbl)
    cum.prop <- cum.bc/length(total)

    if(cum.prop[1] >= quant){
        threshold <- 1
    }else{
        #largest total before quant
        threshold <- max(total.count[cum.prop < quant])
    }
    ind <- total <= threshold
    chip.sum.low <- sum(chip[ind])
    input.sum.low <- sum(input[ind])
    bin.count.low <- cum.bc[which(total.count==threshold)]

    chip <- chip[!ind]
    input <- input[!ind]
    od <- order(chip+input)
    chip <- chip[od]
    input <- input[od]

    #start after threshold
    cum.bc.high <- cum.bc[total.count>threshold]-bin.count.low
    
    all.cum.chip <- cumsum(chip)
    all.cum.input <- cumsum(input)
    cum.chip <- c(0, all.cum.chip[cum.bc.high])
    cum.input <- c(0, all.cum.input[cum.bc.high])
    cum.ratio <- (chip.sum.low+cum.chip)/(input.sum.low+cum.input)

    if(sum(diff(cum.ratio) >= 0) ==0) stop("No increase in cum.ratio, binding signal missing!")
    index <- min((2:length(cum.ratio))[diff(cum.ratio) >= 0])
    return(cum.ratio[index])

}
