## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(BSgenome)
library(BSgenome.Scerevisiae.UCSC.sacCer2)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("BSgenome", "BSgenome.Scerevisiae.UCSC.sacCer2"))

## ----BSgenome------------------------------------------------------------
available.genomes()
library("BSgenome.Scerevisiae.UCSC.sacCer2")
Scerevisiae

## ----BSgenomeLength------------------------------------------------------
seqlengths(Scerevisiae)
seqnames(Scerevisiae)

## ----BSgenomeLoad--------------------------------------------------------
Scerevisiae$chrI

## ----gcChrI--------------------------------------------------------------
letterFrequency(Scerevisiae$chrI, "CG", as.prob = TRUE)

## ----gcGenome------------------------------------------------------------
# set up a new BSparam object
param <- new("BSParams", X = Scerevisiae, FUN = letterFrequency)
unlist(bsapply(param, letters = "GC"))
unlist(bsapply(param, letters = "GC", as.prob=TRUE))
## ----gcGenome2-----------------------------------------------------------
param <- new("BSParams", X = Scerevisiae, FUN = letterFrequency, simplify = TRUE)
bsapply(param, letters = "GC")

## ----gcGenome3-----------------------------------------------------------
sum(bsapply(param, letters = "GC")) / sum(seqlengths(Scerevisiae))

        
## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

