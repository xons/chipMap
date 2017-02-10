#' bw2coverage
#'
#' read bigWig files into coverage RLElist.
#' @import GenomicRanges
#' @import GenomicAlignments
#' @import rtracklayer
#' @param fdir bigWig file path.
#' @return RLElist
#' @export

bw2coverage= function (fdir)
{
	message('Reading in BigWig file. Time consuming.')
	x=import.bw(fdir)
	coverage(x,weight=x$score)
}

#' bed2gr
#'
#' read bed files into GRange object.
#' @import GenomicRanges
#' @import GenomicAlignments
#' @importFrom IRanges IRanges
#' @param bedFile bed file path.
#' @return GRange object
#' @export

bed2gr=function(bedFile)
{
	peak=read.delim(bedFile,header=F)
	peakgr=GRanges(seqnames=peak[,1],ranges=IRanges(start=peak[,2],end=peak[,3]),strand=rep("*",nrow(peak)),peak_index=rownames(peak))
	peakgr
}
