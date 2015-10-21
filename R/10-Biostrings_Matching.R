## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(Biostrings)
library(BSgenome)
library(BSgenome.Scerevisiae.UCSC.sacCer2)


## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("Biostrings", "BSgenome",
##            "BSgenome.Scerevisiae.UCSC.sacCer2", "AnnotationHub"))

## ----mmatchPattern-------------------------------------------------------
dnaseq <- DNAString("ACGTACGT")

# this returns a views object
matchPattern(dnaseq, Scerevisiae$chrI) # single string to single string

countPattern(dnaseq, Scerevisiae$chrI)

# one sequence againsta  set of sequences
# forward strand and reverse strand
vmatchPattern(dnaseq, Scerevisiae)
head(vcountPattern(dnaseq, Scerevisiae))

## ----revCompCheck--------------------------------------------------------
dnaseq == reverseComplement(dnaseq)


# matchPWN - aka sequence logo
# pairwiseAlignment() - for small genome seuqnces usuable in R

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

