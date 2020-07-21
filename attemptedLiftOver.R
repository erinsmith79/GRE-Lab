setwd("/Users/erinsmith/researchlab")

ScGFF <- import(file.path("/Users/erinsmith/researchlab", "GCF_000146045.2_R64_genomic.gff"))
wig_path<-file.path("/Users/erinsmith/researchlab", "nucdata.wig")
nucdatafinal<-import(con=wig_path, seqinfo=seqinfo(genome="sc"))
nucdatafinal
?system.file
SacCerGFF<-import("/Users/erinsmith/researchlab/GCF_000146045.2_R64_genomic.gff",format="GFF")

chain<-import.chain("/Users/erinsmith/researchlab/sacCer3tosacCer1.over.chain.gz")
merged<-as.data.frame(liftOver(SacCerGFF,chain))
?`liftOver,GenomicRanges,Chain-method`
?import.chain