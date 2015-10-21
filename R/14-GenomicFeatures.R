## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

## ----biocLite, eval=FALSE------------------------------------------------
#  source("http://www.bioconductor.org/biocLite.R")
#  biocLite(c("GenomicFeatures", "TxDb.Hsapiens.UCSC.hg19.knownGene"))

## ----txdb----------------------------------------------------------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
txdb

## ----gr------------------------------------------------------------------
gr <- GRanges(seqnames = "chr1", strand = "+", ranges = IRanges(start = 11874, end = 14409))
subsetByOverlaps(genes(txdb), gr)
subsetByOverlaps(genes(txdb), gr, ignore.strand = TRUE)

## ----transcripts---------------------------------------------------------
subsetByOverlaps(transcripts(txdb), gr)

## ----exons---------------------------------------------------------------
subsetByOverlaps(exons(txdb), gr)

## ----exonsBy-------------------------------------------------------------
subsetByOverlaps(exonsBy(txdb, by = "tx"), gr)

# coding sequences 
## ----cds-----------------------------------------------------------------
subsetByOverlaps(cds(txdb), gr)
# what we are really interested in:
# coding sequences overlapping our little fragment
subsetByOverlaps(cdsBy(txdb, by = "tx"), gr)

## ----transcriptLengths---------------------------------------------------
# in terms of spliced RNA. gene id from above
subset(transcriptLengths(txdb, with.cds_len = TRUE), gene_id == "100287102")

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

