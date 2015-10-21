## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("GenomicRanges", "rtracklayer", "AnnotationHub"))
#biocLite(c("rtracklayer"))

## Question 1
# Use the AnnotationHub package to obtain data on "CpG Islands" in the human genome.

## ----ahub_species--------------------------------------------------------
ah <- AnnotationHub()
ah <- subset(ah, species == "Homo sapiens" & genome=="hg19")

## ----ahub_histone--------------------------------------------------------

qcpg <- query(ah, "CpG")
qcpg

head(qcpg)


# Question: How many islands exists on the autosomes?

# no 1 is the one we want, track to download is in field 1
cpg <- qcpg[1][[1]]

# there are not only autosomes in the dataset, but also X, Y chromosomes u M
#refseq <-
cpg <- keepStandardChromosomes(cpg)
seqlevels(cpg)
cpg <- dropSeqlevels(cpg, c("chrX", "chrY", "chrM"))
unique(seqnames(cpg[grep(seqnames(cpg),pattern='^Chr\\d\\d?', ignore.case=TRUE, perl=TRUE)]))

# ANswer
length(cpg)

newStyle <- mapSeqlevels(seqlevels(cpg), "NCBI")
cpg <- renameSeqlevels(cpg, newStyle)
seqlevels(cpg)


# Question 2
# Question: How many CpG Islands exists on chromosome 4.

length(keepSeqlevels(cpg, "chr4"))


## Question 3

# Obtain the data for the H3K4me3 histone modification for the H1 cell line
# from Epigenomics Roadmap, using AnnotationHub. 
# Subset these regions to only keep regions mapped to the autosomes (chromosomes 1 to 22).

# Question: How many bases does these regions cover?


q3 <- query(ah, "H3K4me3")
q3 <- query(q3, "H1")
q3 <- query(q3, "narrowPeak")
q3 <- query(q3, "E003")


gr3 <- subset(q3, title == "E003-H3K4me3.narrowPeak.gz")[[1]]
gr3 <- ah[["AH29884"]]
gr3 <- keepStandardChromosomes(gr3)
gr3 <- dropSeqlevels(gr3, c("chrX", "chrY"))
seqlevels(gr3)

newStyle <- mapSeqlevels(seqlevels(gr3), "NCBI")
gr3 <- renameSeqlevels(gr3, newStyle)
seqlevels(gr3)

genome(gr3) <- "hg19"
# ANswer:
sum(width(gr3))



### Question 4

# Obtain the data for the H3K27me3 histone modification for the H1 cell line from
# Epigenomics Roadmap, using the AnnotationHub package. 
# Subset these regions to only keep regions mapped to the autosomes. 
# In the return data, each region has an associated "signalValue".
# Question: What is the mean signalValue 
#           $ranges)across all regions on the standard chromosomes?

q4 <- query(ah, "H3K27me3")
q4 <- query(q4, "H1")
#q4 <- query(q4, "E003")
gr4 <- subset(q4, title == "E003-H3K27me3.narrowPeak.gz")[[1]]

gr4 <- keepStandardChromosomes(gr4)
seqlevels(gr4)
gr4 <- dropSeqlevels(gr4, c("chrX", "chrY"))
seqlevels(gr4)

gr4
head(gr4, 1)
head(end(gr4), 1)
#Answer: 4.771
mean(gr4$signalValue)


## Question 5

# Bivalent regions are bound by both H3K4me3 and H3K27me3.

# Question: Using the regions we have obtained above, 
# how many bases on the standard chromosomes are bivalently marked?
newStyle <- mapSeqlevels(seqlevels(gr4), "NCBI")
gr4 <- renameSeqlevels(gr4, newStyle)
seqlevels(gr4)

is <- intersect(gr3, gr4)
#Answer:
(sumis <- sum(width(is)))



# Question 6

# We will examine the extent to which bivalent regions overlap CpG Islands.

# Question: how big a fraction (expressed as a number between 0 and 1) 
# of the bivalent regions, overlap one or more CpG Islands?

ov <- findOverlaps(is, cpg)
# ANswer
length(unique(queryHits(ov)))/length(is)


## Question 7

# Question: How big a fraction (expressed as a number between 0 and 1) 
# of the bases which are part of CpG Islands, are also bivalent marked.

# The way I calculated was by dividing  the
# 
#Total number of bases of CpG intersecting k4 and k27
is2 <- intersect(intersect(gr4, gr3), cpg)

#by the
#Total number of bases from CpG islands
(sumcpg <- sum(width(cpg)))

#Answer: 0.2417
sumis2/sumcpg


# 
## Question 8

# Question: How many bases are bivalently marked within 10kb of CpG Islands?
rsz <- resize(cpg, width = 20000+ width(cpg), fix ="center" )
#Answer: 9782086
sum(width(intersect(is, rsz)))

# Tip: consider using the "resize()"" function.


# Question 9
# 
# Question: How big a fraction (expressed as a number between 0 and 1) 
# of the human genome is contained in a CpG Island?
#
#Answer: 0.007
sum(width(cpg))/sum(as.numeric(seqlengths(cpg)))

# Question 10
# 
# Question: Compute an odds-ratio for the overlap of bivalent marks
# with CpG islands.
overlapMat <- matrix(0, nrow=2, ncol=2)
overlapMat[1,1] <- sumis2

overlapMat[1,2] <- sum(width(setdiff(cpg,is2)))
overlapMat[2,1] <- sum(width(setdiff(is,is2)))
overlapMat[2,2] <- sum(as.numeric(seqlengths(cpg))) - sum(overlapMat)

overlapMat
# Answer: 169.1
(oddsRatio <- overlapMat[1,1] * overlapMat[2,2] / (overlapMat[2,1] * overlapMat[1,2]))


# overlapMat[1,1] <- sum(width(both))
# overlapMat[1,2] <- sum(width(only.prom))
# overlapMat[2,1] <- sum(width(only.peaks))
# overlapMat[2,2] <- 3*10^9 - sum(overlapMat)

# 160.4022
# 169.0962
# 138.4391
# 211.0553




