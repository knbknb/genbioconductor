## ----dependencies, warning=FALSE, message=FALSE--------------------------
#
library(GenomeInfoDb)
library(GenomicRanges)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("GenomeInfoDb", "GenomicRanges"))

## ----seqlevelsForce------------------------------------------------------
gr <- GRanges(seqnames = c("chr1", "chr2"),
              ranges = IRanges(start = 1:2, end = 4:5))
# tow chromosomes
gr
# drop everything except chr1
seqlevels(gr, force=TRUE) <- "chr1"
gr

# same feature
## ----dropSeqlevels-------------------------------------------------------
gr <- GRanges(seqnames = c("chr1", "chr2"),
              ranges = IRanges(start = 1:2, end = 4:5))
dropSeqlevels(gr, "chr1")
keepSeqlevels(gr, "chr2")

## ----keepStandard--------------------------------------------------------
gr <- GRanges(seqnames = c("chr1", "chrU345"),
              ranges = IRanges(start = 1:2, end = 4:5))
keepStandardChromosomes(gr)

## ----GRanges-------------------------------------------------------------
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = 1:2, width = 2))
# there are inconsistent naming conventions for
# chromosomes: chr1, chr01, 1, I, chrI

## ----seqStyle------------------------------------------------------------
# assigns new styles / naming conventions
newStyle <- mapSeqlevels(seqlevels(gr), "NCBI")
gr <- renameSeqlevels(gr, newStyle)

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

# ?? seqinfo(bsgenomeName())
