#!/usr/bin/env Rscript

## R-script for adding q-values

library(qvalue)

FILEXT = '.tsv.gz'



## main
sessionInfo()
.Machine
args.v <- commandArgs(trailingOnly = TRUE)
cat("invoked R like this:\n")
cat(paste(args.v, collapse=" "), "\n")
args.n = length(args.v)

if (args.n == 2) {
  infnam <- args.v[1]
  oufnam <- args.v[2]
} else {
  stop("Call with <tensortqtl table [TSV.GZ]> <output file [TSV.GZ]\n", call.=FALSE)
}

# assume file name is gzipped TSV
dat <- read.delim(gzfile(infnam,'rt'), row.names = NULL, quote = "")

qval <- qvalue(dat$pval_beta)

dat$qvalueR <- qval$qvalue

oufh = gzfile(oufnam, 'wt')
write.table(dat, file = oufh, sep = "\t", row.names = T, quote = FALSE)
flush(oufh)
close(oufh)
