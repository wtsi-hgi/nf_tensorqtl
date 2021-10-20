#!/usr/bin/env Rscript

## R-script for normalizing pseudo-bulk RNA-seq counts (dSum, TMM)
## using edgeR::calcNormFactors()

## assume there are donors to be excluded who have a number of cells below a specified threshold

DEFAULT_MIN_CELLS_PER_DONOR <- 5 # remove samples with fewer than this number of cells

library(edgeR)

write.singleBED <- function(df, oufn = "test.bed")
{
  names(df)[1] <- "#chr"
  write.table(df, file = oufn, sep = "\t", row.names = FALSE, quote = FALSE)
  cat ("Wrote normalised counts to file '", oufn, "' ...\n", sep="")
}

write.perChrBED <- function(df, dirnam)
{
  chrnams <- factor(df$chr)
  dl <- split(df, chrnams)
  for (cn in levels(chrnams)) {
    cat("chromosome: ", cn, "\n")
    oufpath = file.path(dirnam, paste0(dirnam, "_", cn, ".bed"))
    cat(oufpath, "\n")
    dd <- dl[[cn]]
    write.singleBED(dl[[cn]])
  }
  dl
}
write.principal.components <- function(df, dirnam, n.pcs = 20)
{
  oufn <- paste0(dirnam, "_PC", n.pcs, ".tsv")
  n_donors <- dim(df)[2]
  m <- t(df[,5:n_donors])
  pca <- prcomp(m, scale. = TRUE, rank. = n.pcs)
  write.table(pca$x, file = oufn, sep = "\t", row.names = TRUE, quote = FALSE)
  pca
}

normalise.dSumTMM <- function(phenotype.bed.file, donor.tsv.file, min.cells.per.donor = 5)
{
  #oufnprfx <- sub("\\.bed$", "", phenotype.bed.file)

  # read table with the number of cells per donor
  df.donors <- read.delim(donor.tsv.file, row.names = NULL, quote = "")

sample.names.included <- df.donors$sample_name[
  df.donors$donor_included == "True" &
  df.donors$n_cells_pass_qc >= min.cells.per.donor
  ]

  dat.counts <- read.delim(phenotype.bed.file, row.names = NULL, quote = "")
  row.names(dat.counts) <- dat.counts[,4]

  ## subset matrix to remove donors which fewer than MIN_CELLS_PER_DONOR cells
  n_cols = dim(dat.counts)[2]
  m.sum <- as.matrix(dat.counts[,5:n_cols])[, as.vector(sample.names.included)]

  # TMM normalization
  dge <- edgeR::DGEList(counts=m.sum)
  dge <- edgeR::calcNormFactors(dge)

  # defaults:
  # calcNormFactors(
  #  m.sum,
  #  lib.size = None, method = "TMM", refColumn = NULL,
  #  logratioTrim = 0.3, sumTrim = 0.5, doWeighting = TRUE,
  #  Acutoff = -1e10
  #  )

  # log2(CPM) as output
  m.logcpm <- edgeR::cpm(dge, log=TRUE)

  df <- cbind(dat.counts[,1:4], as.data.frame(m.logcpm))
}

## main
## infnam.bed <- "count_all.bed.gz"
## infnam.tsv <- "ncells_all.tsv"
args.v <- commandArgs(trailingOnly = TRUE)
cat("invoked R like this:\n");
cat(paste(args.v, collapse=" "), "\n")
args.n = length(args.v)
if (args.n > 1 & args.n < 5) {
  infnam.bed <- args.v[1]
  infnam.tsv <- args.v[2]
  if (args.n == 3) {
    n.cells.min <- args.v[3]
  } else {
    n.cells.min <- DEFAULT_MIN_CELLS_PER_DONOR
  }
} else {
  stop("Call with <input BED file> <input TSV table> <min_cell_mum>\n", call.=FALSE)
}

cat ("bed file:", infnam.bed, ", tsv table:", infnam.tsv, "n.cells.min =", n.cells.min, "\n")
dirnam <- paste0("norm", sub("\\.bed.gz$", basename(infnam.bed), replacement=""))
if (!file.exists(dirnam)) {
  dir.create(dirnam)
}

df <- normalise.dSumTMM(
  phenotype.bed.file = infnam.bed,
  donor.tsv.file = infnam.tsv,
  min.cells.per.donor = n.cells.min
  )
write.principal.components(df, dirnam, n.pcs = 20)

if (args.n == 4 & args.v[4] == 'perChr') {
  write.perChrBED(df, dirnam)
} else {
  write.singleBED(df, oufn = paste0(dirnam, "_chrAll.bed"))
}
cat ("Wrote", length(df), "files.\n")
