#' annPeak
#'
#' Annotate peak by it's nearest transcript/gene.
#' If peak overlaps any TSS, peak_to_feature will be 'overlap_tss'. If it overlaps multiple tss, random one will be selected.
#' The rest of the peak region will be assigned to it's nearest gene (not nearest tss). If it overlaps a gene body, the distance will be 0.
#' @import GenomicRanges
#' @param peak.gr
#' A peak GRanges. This can be get from Bed file of peak using bed2gr function in this package.
#' @param tx
#' A data frame of transcript/gene annotation, contains columns 'txid','seqnames', 'start', 'end', 'strand'.
#' @param tssup
#' Window size around TSS that will be defined as 'overlap_tss'. Default is 1000, means must overlap within 1 kb upstream of TSS.
#' @param tssdown
#' Window size around TSS that will be defined as 'overlap_tss'. Default is 1000, means must overlap within 1 kb downstream of TSS.
#' @return A dataframe of annotated peaks with peak and it's nearest gene information.
#' @export

annPeak=function(peak.gr,tx,tssup=1000,tssdown=1000)
{
	tx.gr=with(tx, GRanges(
			seqnames=seqnames,
			ranges=IRanges(start=start,end=end),
			strand=strand,
			txid=txid))
	##names(tx.gr)=tx$txid

	hits=distanceToNearest(peak.gr,tx.gr,ignore.strand = TRUE) ## select random gene if overlap multiple
	hits=as.data.frame(hits)
	ann=cbind(as.data.frame(peak.gr),as.data.frame(tx.gr)[hits$subjectHits,],distance=hits$distance)
	names(ann)[2:5]=c("peak_start", "peak_end", "peak_width", "peak_strand")
	
	ann$peak_to_feature="body"
	up=which((ann$strand=="+" & ann$start > ann$peak_end  ) |(ann$strand=="-" & ann$peak_start > ann$end ))
	if (length(up) != 0) ann[up,]$peak_to_feature="upstream"
	down=which((ann$strand=="+" & ann$peak_start > ann$end ) |(ann$strand=="-" & ann$start > ann$peak_end ))
	if (length(down) != 0) ann[down,]$peak_to_feature="downstream"

	## correct tss overlap
	s<- as(strand(tx.gr),"vector")=="+"
	tss=ifelse(s,start(tx.gr),end(tx.gr))
	tss.gr=GRanges(seqnames=seqnames(tx.gr),ranges=IRanges(start=ifelse(s,tss-tssup,tss-tssdown),end=ifelse(s,tss+tssdown,tss+tssup)),strand(tx.gr),tr.id=tx.gr$tr.id)

	hits_tss=distanceToNearest(peak.gr,tss.gr,ignore.strand = TRUE)## select random gene if overlap multiple
	hits_tss=as.data.frame(hits_tss)
	ann_tss=cbind(as.data.frame(peak.gr),as.data.frame(tx.gr)[hits_tss$subjectHits,],distance=hits_tss$distance)
	ann_tss$peak_to_feature="overlap_tss"

	oltss=which(ann_tss$distance==0)
	if (length(oltss)!=0) ann[oltss,]=ann_tss[oltss,]
	
	ann
}

#' annPeak_bytss
#'
#' Annotate peak by it's nearest transcription start site.
#' If peak overlaps multiple tss, random one will be selected.
#' @import GenomicRanges
#' @param peak.gr
#' A peak GRanges. This can be get from Bed file of peak using bed2gr function in this package.
#' @param tx
#' A data frame of transcript/gene annotation, contains columns 'txid','seqnames', 'start', 'end', 'strand'.
#' @param tssup
#' Window size around TSS that will be defined as overlap_tss, and 0 distance will be given. Default is 1000, means must overlap within 1 kb upstream of TSS.
#' @param tssdown
#' Window size around TSS that will be defined as overlap_tss, and 0 distance will be given. Default is 1000, means must overlap within 1 kb downstream of TSS.
#' @return A dataframe of annotated peaks with peak and it's nearest TSS information.
#' @export

annPeak_bytss=function(peak.gr,tx,tssup=1000,tssdown=1000)
{
	tx.gr=with(tx, GRanges(
			seqnames=seqnames,
			ranges=IRanges(start=start,end=end),
			strand=strand,
			txid=txid))
	##names(tx.gr)=tx$txid

	s<- as(strand(tx.gr),"vector")=="+"
	tss=ifelse(s,start(tx.gr),end(tx.gr))
	tss.gr=GRanges(seqnames=seqnames(tx.gr),ranges=IRanges(start=ifelse(s,tss-tssup,tss-tssdown),end=ifelse(s,tss+tssdown,tss+tssup)),strand(tx.gr),tr.id=tx.gr$tr.id)

	hits_tss=distanceToNearest(peak.gr,tss.gr,ignore.strand = TRUE)
	hits_tss=as.data.frame(hits_tss)
	ann_tss=cbind(as.data.frame(peak.gr),as.data.frame(tx.gr)[hits_tss$subjectHits,],distance=hits_tss$distance)
	names(ann_tss)[2:5]=c("peak_start", "peak_end", "peak_width", "peak_strand")

	## calculate distance_center_to_tss for sorting in the future.
	ann_tss$tss=ifelse(ann_tss$strand=="+",ann_tss$start,ann_tss$end)
	ann_tss$center=round((ann_tss$peak_start+ann_tss$peak_end)/2)
	ann_tss$distance_center_to_tss=ifelse(ann_tss$strand=="+", ann_tss$center-ann_tss$tss, ann_tss$tss - ann_tss$center)
	if (min(ann_tss$distance)== 0) ann_tss[which(ann_tss$distance==0),]$distance_center_to_tss= 0
	
	ann_tss
}

#' plotPie
#'
#' plot pie chart of 'peak_to_feature' distribution of annPeak results.
#' @param annp
#' A dataframe of annPeak output.
#' @param titles
#' Figure title and output file name
#' @param out.pdf
#' Save image as PDF or PNG. Default is PNG.
#' @param color
#' color in image.
#' @return A image of PNG or PDF.
#' @export
plotPie=function(annp,titles,out.pdf=F,color=c("yellow2","orange","red","red4"))
{
	N=nrow(annp)
	n= summary(as.factor(annp$peak_to_feature))[c(3,1,2,4)]

	if (out.pdf==T) pdf(paste0("pie.",titles,".pdf")) else
	png(paste0("pie.",titles,".png"))
		par(mar=c(0,0,4,0))
		pie(n,labels="",col=color,border=NA,clockwise = T,init.angle=180,main=titles)
		points(0,0,pch=16,cex=par("cex")*35,col="white")
		l=paste0(names(n)," = ",round(n/N*100),"%")
		l=c(paste0("N = ",N),"Peak to features:",l)
		legend(-0.4,0.3,legend=l,ncol=1,col=c("white","white",color),pch=15,bty="n",cex=1.3)
	dev.off()
}

#' plotPie_tss
#'
#' plot pie chart of 'distance_center_to_tss' distribution of annPeak_bytss results.
#' @import RColorBrewer
#' @param annp_bytss
#' A dataframe of annPeak_bytss output.
#' @param titles
#' Figure title and output file name
#' @param out.pdf
#' Save image as PDF or PNG. Default is PNG.
#' @param color
#' color in image.
#' @return A image of PNG or PDF.
#' @export

plotPie_tss=function(annp_bytss,titles,out.pdf=F,color=rev(brewer.pal(5,"Blues")))
{
	x=annp_bytss$distance_center_to_tss/1000
	g=cut(x,breaks=c(min(x),-50,-10,-5,-1,0,1,5,10,50,max(x)))
	xg=split(x,g)
	n=unlist(lapply(xg,length))
	dis=n/sum(n)
	dis2=dis[6:10]+dis[5:1]
	N=nrow(annp_bytss)
	
	if (out.pdf==T) pdf(paste0("pie.tss.",titles,".pdf")) else
	png(paste0("pie.tss",titles,".png"))
		par(mar=c(0,0,4,0))
		pie(dis2,labels="",col=color,border=NA,clockwise = T,init.angle=180,main=titles)
		points(0,0,pch=16,cex=par("cex")*35,col="white")
		l=paste0(c("<1 kb","1-5 kb","5-10 kb","10-50 kb",">50 kb")," : ",round(dis2*100),"%")
		l=c(paste0("N = ",N),"Peak center to TSS:",l)
		legend(-0.4,0.36,legend=l,ncol=1,col=c("white","white",color),pch=15,bty="n",cex=1.3)
	dev.off()
 }

#' nearest_peak
#'
#' Find nearest peak of given genes's transcription start site.
#' If tss overlaps multiple peaks, random one will be selected.
#' @import GenomicRanges
#' @param peak.gr
#' A peak GRanges. This can be get from Bed file of peak using bed2gr function in this package.
#' @param tx
#' A data frame of transcript/gene annotation.
#' @param tssup
#' Window size around TSS that will be defined as overlap_tss, and 0 distance will be given. Default is 1000, means must overlap within 1 kb upstream of TSS.
#' @param tssdown
#' Window size around TSS that will be defined as overlap_tss, and 0 distance will be given. Default is 1000, means must overlap within 1 kb downstream of TSS.
#' @return A dataframe with gene and it's nearest peak information.
#' @export

nearest_peak=function(peak.gr,tx,tssup=1000,tssdown=1000)
{
	tx.gr=with(tx, GRanges(
			seqnames=seqnames,
			ranges=IRanges(start=start,end=end),
			strand=strand,
			txid=txid))

	s<- as(strand(tx.gr),"vector")=="+"
	tss=ifelse(s,start(tx.gr),end(tx.gr))
	tss.gr=GRanges(seqnames=seqnames(tx.gr),ranges=IRanges(start=ifelse(s,tss-tssup,tss-tssdown),end=ifelse(s,tss+tssdown,tss+tssup)),strand(tx.gr),id=names(tx.gr))

	hits=distanceToNearest(tss.gr,peak.gr,ignore.strand = TRUE)
	hits=as.data.frame(hits)
	ann=cbind(as.data.frame(tx.gr)[hits$queryHits,],as.data.frame(peak.gr)[hits$subjectHits,],distance=hits$distance)
	colnames(ann)[7:11]=paste0('peak_',colnames(ann)[7:11])
	
	ann
}

