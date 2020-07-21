saccer1<-read.delim(file.choose("SacCer1.fa"))
setwd("/Users/erinsmith/researchlab")
datafile<-read.delim(file.choose("nucdata3.txt"))
colnames(datafile)
colnames(datafile, do.NULL=FALSE)
colnames(datafile) <-c("chromosome","position","value")
head(datafile)


graphnucdata<- function(dat=datafile, chrom="2",startpos=1,endpos=1000) {
  chromSpecifDat<-subset(dat, chromosome==chrom)
  position<-as.numeric(chromSpecifDat$position)
  positionSpecifDat<-subset(chromSpecifDat, position>=startpos & position<=endpos)
  final<-plot(dat$position, dat$value, pch=19, main=x, ylab="data value", type="l")
}
graphnucdata()
summary(datafile)
library(GenomicRanges)
library(IRanges)
library(rtracklayer)
plotTracks(list(nucdatafinal))
nucdataplot <- function(data= nucdatafinal, seqnam="chr1", start=1, end=1000){
  sub_data <- data[data$seqnames == seqnam,]
  position <- as.numeric(sub_data$position)
  sub_data <- sub_data[position >= start & position <= end,]
  x <-seqnam
  with(sub_data, plot(position, score, pch = 19, main = x, ylab= "data value", type = "l", xlim=c(0,1000000), ylim=c(0,1000000)))
} 
nucdataplot(data=nucdatafinal, seqnam="chr2")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tracklayer")
library(Biostrings)

saccergenome<-import.2bit("sacCer1.2bit")
plot(consensusdata$ends, consensusdata$ByPos_MIndex)

consensusdata<-vmatchPattern("CACGTG",saccergenome,max.mismatch=0,min.mismatch=0,with.indels=FALSE,fixed=TRUE,algorithm="auto")
SCprom<-promoters(SacCerGFF)
colnames(mcols(consensusdata))
promoterconsensusdata<-vmatchPattern("CACGTG",SCprom,max.mismatch=0,min.mismatch=0,with.indels=FALSE,fixed=TRUE,algorithm="auto")
plot(consensusdata$MIndex, consensusdata$ends)

