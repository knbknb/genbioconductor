## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(IRanges)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("IRanges"))

## ----iranges1------------------------------------------------------------
# give two arguments, last 3rd can be inferred.
ir1 <- IRanges(start = c(1,3,5), end = c(3,5,7))
ir1
ir2 <- IRanges(start = c(1,3,5), width = 3)
all.equal(ir1, ir2)

## access elements
## ----ir_width------------------------------------------------------------
start(ir1)
width(ir2) <- 1
ir2

## ----ir_names------------------------------------------------------------
names(ir1) <- paste("A", 1:3, sep = "")
ir1

## ----ir_dim--------------------------------------------------------------
dim(ir1)
length(ir1)

## ----ir_subset-----------------------------------------------------------
ir1[1]
ir1["A1"]

## ----concatenate---------------------------------------------------------
c(ir1, ir2)

## "Normal" IRanges
## ----irNormal1, echo=FALSE-----------------------------------------------
ir <- IRanges(start = c(1,3,7,9), end = c(4,4,8,10))
source("plotRanges.R")
## ----irNormal2, echo=FALSE, fig.height=2, small.mar=TRUE-----------------
plotRanges(ir)

## ----irNormal3, echo=FALSE, fig.height=1.75, small.mar=TRUE--------------
# minimal representation of the ranges as a set
plotRanges(reduce(ir))

## ----irNormal4-----------------------------------------------------------
ir
reduce(ir)

## ----irDisjoin1, eval=FALSE----------------------------------------------
# handy when you need it
# disjoin(ir1)

## ----irDisjoin2, echo=FALSE, fig.height=2, small.mar=TRUE----------------
plotRanges(ir)

## ----irDisjoin3, echo=FALSE, fig.height=1.75, small.mar=TRUE-------------
plotRanges(disjoin(ir))

## ----ir_resize-----------------------------------------------------------
ir <- IRanges(start = c(1,3,7,9), end = c(4,4,8,10))
plotRanges(ir)
rsz1 <-resize(ir, width = 1, fix = "start")
plotRanges(rsz1)


rsz2 <- resize(ir, width = 1, fix = "center")
plotRanges(rsz2)
## ----ir_sets-------------------------------------------------------------
ir1 <- IRanges(start = c(1, 3, 5), width = 1)
ir2 <- IRanges(start = c(4, 5, 6), width = 1)
plotRanges(ir1)
plotRanges(ir2)

ir3 <- union(ir1, ir2)
ir4 <-intersect(ir1, ir2)

plotRanges(ir3)
plotRanges(ir4)

## ----union2--------------------------------------------------------------
reduce(c(ir1, ir2))

## ----findOverlaps--------------------------------------------------------
ir1 <- IRanges(start = c(1,4,8), end = c(3,7,10))
ir2 <- IRanges(start = c(3,4), width = 3, )
par(mfrow=c(1,1))
dev.new()
plotRanges(ir1)
plotRanges(ir2)
ov <- findOverlaps(ir1, ir2)
ov

#Hits object with 3 hits and 0 metadata columns:
#        queryHits subjectHits
#        <integer>   <integer>
#[1]         1           1
#[2]         2           1
#[3]         2           2
# Index 1 So the first row, or
# the first element of the Hit object,
# means that range number one, in the query,
# overlaps range number one in the subject.

## ----findOverlaps_ill----------------------------------------------------
intersect(ir1[subjectHits(ov)[1]],
          ir2[queryHits(ov)[2]])

## ----subjectHits---------------------------------------------------------
queryHits(ov)
unique(queryHits(ov))

## ----argsFindOverlaps, tidy=TRUE-----------------------------------------
args(findOverlaps)

## ----countOverlaps-------------------------------------------------------
countOverlaps(ir1, ir2)
# countOverlaps  returns a vector.
# In this case, it means that range
# number one in the ir1 overlaps ir2[1].
# Element two overlaps two elements of ir2.
# [1] 1 2 0

## ----nearest-------------------------------------------------------------
ir1
ir2
nearest(ir1, ir2)

# Finally we can also relate IRanges in
# a different way than through the overlaps.
# We can look at which ones
# are close to each other.
# So again, we take our two IRanges and
# we can ask,
# which of these IRanges in ir2
# are closer to the ones in ir1?

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

