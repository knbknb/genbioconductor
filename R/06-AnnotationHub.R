## ----dependencies, warning=FALSE, message=FALSE--------------------------
# interfaace to a lot of online databases
# relatively new package, powerful
library(AnnotationHub)

## ----biocLite, eval=FALSE------------------------------------------------
## source("http://www.bioconductor.org/biocLite.R")
## biocLite(c("AnnotationHub"))

## ----annoHub,results="hide"----------------------------------------------
ah <- AnnotationHub()
ah

## ----no1-----------------------------------------------------------------
ah[1]

## ----look----------------------------------------------------------------
unique(ah$dataprovider)
unique(ah$rdataclass)

## ----subset--------------------------------------------------------------
ah <- subset(ah, species == "Homo sapiens")
ah

## ----query---------------------------------------------------------------
query(ah, "H3K4me3")

## ----display, eval=FALSE-------------------------------------------------
## hist <- display(ah)

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

