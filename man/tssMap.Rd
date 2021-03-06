% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chipMap.R
\name{tssMap}
\alias{tssMap}
\title{tssMap}
\usage{
tssMap(query, anndf, left = -5000, right = 5000, bin_width = 25)
}
\arguments{
\item{query}{A coverage RLEList. This can be get from a bigWig file using bw2coverage function in this package. BigWig file can be generated from Bam file using bamtools and bedtools. Or a peak GRanges. This can be get from import Bed file of peak into GRanges object. If the query is a peak GRanges, the output will be a enrichment profile instead of a coverage heatmap.}

\item{anndf}{A data frame contains columns 'seqnames', 'start', 'end', 'strand'. The function will find TSS based on the strand information.}

\item{left}{Specify the window around TSS to be calculate. Default is -5 kb to 5 kb.}

\item{right}{Specify the window around TSS to be calculate. Default is -5 kb to 5 kb.}

\item{bin_width}{Specify the bin size to smooth and compress the window. Default is 25.}
}
\value{
A matrix with each line as a region around each TSS and each column is a bin.
}
\description{
Calculate chip seq matrix around given tss. Flip starnd if strand is "-".
}

