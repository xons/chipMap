#' tssMap
#'
#' Calculate chip seq matrix around given tss. Flip starnd if strand is "-".
#' @import GenomicRanges
#' @import GenomicAlignments
#' @param query
#' A coverage RLEList. This can be get from a bigWig file using bw2coverage function in this package. BigWig file can be generated from Bam file using bamtools and bedtools. 
#' Or a peak GRanges. This can be get from import Bed file of peak into GRanges object. If the query is a peak GRanges, the output will be a enrichment profile instead of a coverage heatmap.
#' @param anndf
#' A data frame contains columns 'seqnames', 'start', 'end', 'strand'. The function will find TSS based on the strand information.
#' @param left
#' Specify the window around TSS to be calculate. Default is -5 kb to 5 kb.
#' @param right
#' Specify the window around TSS to be calculate. Default is -5 kb to 5 kb.
#' @param bin_width
#' Specify the bin size to smooth and compress the window. Default is 25.
#' @return A matrix with each line as a region around each TSS and each column is a bin.
#' @export

tssMap=function(query,anndf,left=-5000,right=5000,bin_width=25)
{
	ann=anndf
	Centres=ifelse(ann$strand==1|ann$strand=="+",ann$start,ann$end)
	ann$start=Centres+left
	ann$end=Centres+right

	strand=ann$strand

	#enrichment or coverage plot
	enrich <- class(query)=="GRanges"
	#if (class(query)=="GRanges") enrich=TRUE
	if (enrich) cvg=coverage(query) else cvg=query

	chr=as(ann$seqnames,"vector")
	nrows=length(strand)
	ncols=ifelse(enrich | bin_width==1,(right-left)/bin_width+1,(right-left)/bin_width)
	v = array( 0, c(nrows,ncols) )

	message("Calculating tssMap. Time consuming.")
	for( i in 1:nrow(ann) )
	{
		vtmp=rep(0,right-left+1) # prepare a empty matrix for saving data.
		chrcvg=cvg[[chr[i]]]
		# set to zeros if the window exit chromosome length.
		start=max(ann$start[i],1)
		end=min(ann$end[i],length(chrcvg))
		if (start<=length(chrcvg))
		{
			tmp2=ifelse(ann$start[i]<1,abs(ann$start[i])+2,1)
			tmp3=end-start+tmp2

			vtmp[tmp2:tmp3]= as(chrcvg[start:end],"vector")
			if (enrich | bin_width==1) v[i,]=bin(vtmp,bin_width) else v[i,]=bin_mean(vtmp,bin_width)

			if( strand[i]== -1| strand[i]=="-" ) v[i,] = rev( v[i,] )
		}
	}
	return(v)
}

#' bin
#'
#' Bin the vector by select the first element of the vector
#' @param v numeric vector
#' @param bin size of bin
#' @return A vector
bin=function(v,bin)
{
	bincol=seq(1,length(v),bin)
	out=v[bincol]
	return(out)
}


#' bin_mean
#'
#' Bin the vector by return the mean of the vector
#' @param v numeric vector
#' @param bin size of bin
#' @return A vector
bin_mean=function(v,bin)
{
	out=colMeans(array(v,c(bin,length(v)/bin)))
	return(out)
}

#' peakMap
#'
#' Calculate chip seq matrix around the centre of given peaks.
#' @import GenomicRanges
#' @import GenomicAlignments
#' @param query
#' A coverage RLE. This can be get from a bigWig file using bw2coverage function in this package. BigWig file can be generated from Bam file using bamtools and bedtools.
#' Or a peak GRanges. This can be get from Bed file of peak using bed2gr function in this package. If the query is a peak GRanges, the output will be a enrichment profile instead of a coverage heatmap.
#' @param anndf A data frame contains columns 'seqnames', 'peak_start', 'peak_end', 'strand'. The function will find peak center. 'strand' can be '*'. Flip strand if strand is "-".
#' @param left
#' Specify the window around TSS to be calculate. Default is -5 kb to 5 kb.
#' @param right
#' Specify the window around TSS to be calculate. Default is -5 kb to 5 kb.
#' @param bin_width
#' Specify the bin size to smooth and compress the window. Default is 25.
#' @param sortByWidth
#' If sort by the peak width or not. Default if FALSE.
#' @return A matrix with each line as a region around each peak centre and each column is a bin.
#' @export

peakMap=function(query,anndf,left=-5000,right=5000,bin_width=25,sortByWidth=FALSE)
{
	if (sortByWidth==TRUE)
	{
		ann=anndf
		ann=ann[order(ann$peak_end-ann$peak_start,decreasing=TRUE),]
	} else
	ann=anndf
	Centres=floor((ann$peak_start+ann$peak_end)/2)
	ann$start=Centres+left
	ann$end=Centres+right

	strand=ann$strand

	#enrichment or coverage plot
	enrich <- class(query)=="GRanges"
	#if (class(query)=="GRanges") enrich=TRUE
	if (enrich) cvg=coverage(query) else cvg=query

	chr=as(ann$seqnames,"vector")
	nrows=length(strand)
	ncols=ifelse(enrich | bin_width==1,(right-left)/bin_width+1,(right-left)/bin_width)
	v = array( 0, c(nrows,ncols) )

	message("Calculating peakMap. Time consuming.")
	for( i in 1:nrow(ann) )
	{
		vtmp=rep(0,right-left+1) # prepare a empty matrix for saving data.
		chrcvg=cvg[[chr[i]]]
		# set to zeros if the window exit chromosome length.
		start=max(ann$start[i],1)
		end=min(ann$end[i],length(chrcvg))
		if (start<=length(chrcvg))
		{
		 tmp2=ifelse(ann$start[i]<1,abs(ann$start[i])+2,1)
		 tmp3=end-start+tmp2

		 vtmp[tmp2:tmp3]= as(chrcvg[start:end],"vector")
		 if (enrich | bin_width==1) v[i,]=bin(vtmp,bin_width) else v[i,]=bin_mean(vtmp,bin_width)

 		if( strand[i]== -1| strand[i]=="-" ) v[i,] = rev( v[i,] )
		}
 	}
	v
}




