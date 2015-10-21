## ----dependencies, warning=FALSE, message=FALSE--------------------------
library(GenomicRanges)
library(airway)

## ----biocLite, eval=FALSE------------------------------------------------
source("http://www.bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
#biocLite("GenomicRanges")

#c('abind', 'acepack', 'affy', 'affyio', 'ALL', 'annotate', 'AnnotationDbi', 'AnnotationHub', 'ape', 'assertthat', 'base64enc', 'BatchJobs', 'BBmisc', 'BH', 'Biobase', 'BiocGenerics', 'BiocInstaller', 'BiocParallel', 'biomaRt', 'Biostrings', 'bitops', 'boot', 'brew', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Scerevisiae.UCSC.sacCer2', 'car', 'caTools', 'checkmate', 'class', 'cluster', 'codetools', 'colorspace', 'curl', 'DBI', 'DESeq', 'devtools', 'dichromat', 'digest', 'dplyr', 'e1071', 'evaluate', 'fail', 'fastcluster', 'foreign', 'formatR', 'Formula', 'futile.logger', 'futile.options', 'gdata', 'genefilter', 'geneplotter', 'GenomeInfoDb', 'GenomicRanges', 'GEOquery', 'ggplot2', 'ggvis', 'git2r', 'gridExtra', 'gtable', 'gtools', 'hexbin', 'hgu95av2.db', 'highr', 'HistData', 'Hmisc', 'htmltools', 'httpuv', 'httr', 'hwriter', 'interactiveDisplayBase', 'IRanges', 'jsonlite', 'KernSmooth', 'knitr', 'labeling', 'lambda.r', 'lattice', 'latticeExtra', 'lazyeval', 'limma', 'lme4', 'locfit', 'lubridate', 'magrittr', 'maps', 'markdown', 'MASS', 'Matrix', 'MatrixModels', 'matrixStats', 'memoise', 'mgcv', 'mime', 'minqa', 'munsell', 'nlme', 'nloptr', 'NLP', 'nnet', 'nutshell', 'nutshell.audioscrobbler', 'nutshell.bbdb', 'org.Hs.eg.db', 'pbkrtest', 'pheatmap', 'plyr', 'preprocessCore', 'proto', 'quantreg', 'R6', 'rafalib', 'Rbowtie', 'RcmdrMisc', 'RColorBrewer', 'Rcpp', 'RcppArmadillo', 'RcppEigen', 'RCurl', 'readxl', 'reshape2', 'rgl', 'rJava', 'rmarkdown', 'RMySQL', 'roxygen2', 'rpart', 'RSQLite', 'rstudioapi', 'rversions', 'S4Vectors', S4Vectors-'old', 'sandwich', 'scales', 'sendmailR', 'shiny', 'slam', 'snow', 'sp', 'SparseM', 'spatial', 'statmod', 'stringi', 'stringr', 'survival', 'tm', 'whisker', 'XML', 'xml2', 'xtable', 'XVector', 'yaml', 'zlibbioc', 'zoo')
## ----airway--------------------------------------------------------------
library(airway)
data(airway)
airway
n
## ----colData-------------------------------------------------------------
colData(airway)

## ----getColumn-----------------------------------------------------------
airway$cell

## ----exptData------------------------------------------------------------
exptData(airway)

## ----names---------------------------------------------------------------
colnames(airway)
head(rownames(airway))

## ----assay---------------------------------------------------------------
airway
assayNames(airway)
assays(airway)
head(assay(airway, "counts"))

## ----rowRanges-----------------------------------------------------------
length(rowRanges(airway))
dim(airway)
rowRanges(airway)

## ----numberOfExons-------------------------------------------------------
length(rowRanges(airway))
sum(elementLengths(rowRanges(airway)))

## ----start---------------------------------------------------------------
start(rowRanges(airway))
start(airway)

## ----subsetByOverlaps----------------------------------------------------
gr <- GRanges(seqnames = "1", ranges = IRanges(start = 1, end = 10^7))
subsetByOverlaps(airway, gr)

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

