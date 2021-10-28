#!/usr/bin/env Rscript

## R-script for normalizing pseudo-bulk RNA-seq counts (dSum, TMM)
## using edgeR::calcNormFactors()

## assume there are donors to be excluded who have a number of cells below a specified threshold

DEFAULT_MIN_CELLS_PER_DONOR <- 5 # remove samples with fewer than this number of cells
DEFAULT_NUM_PCS_WRITTEN <- 20 # number of principal components to write out
COLUMN_GENE_ID <- 4
VARIANCE_ZERO <- 10e-7

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
}

write.principal.components <- function(df, dirnam, n.pcs = DEFAULT_NUM_PCS_WRITTEN)
{
  oufn <- paste0(dirnam, "_PC", n.pcs, ".tsv")
  cat("writing file", oufn, "...\n")
  n_donors <- dim(df)[2]
  m <- t(df[,(COLUMN_GENE_ID+1):n_donors])
  # remove zero-variance columns
  vv <- apply(m,2,var)
  z.v = vv < VARIANCE_ZERO
  nz = sum(z.v)
  cat("There are ", nz, "columns with zero variance.\n")
  if (nz > 0) {
    m <- m[, which(z.v)]
    }
  pca <- prcomp(m, center = TRUE, scale. = TRUE)
  nc <- dim(pca$x)[1]
  if (nc > n.pcs) {
    nc <- n.pcs
  }
  write.table(pca$x[,1:nc], file = oufn, sep = "\t", row.names = TRUE, quote = FALSE)
  pca
}

normalise.dSumTMM <- function(phenotype.bed.file, donor.tsv.file, min.cells.per.donor = 5)
{
  #oufnprfx <- sub("\\.bed$", "", phenotype.bed.file)

  # read table with the number of cells per donor
  df.donors <- read.delim(donor.tsv.file, row.names = NULL, quote = "")

  sample.names.included <- df.donors$donor_id[
    df.donors$is_included == "True" &
    df.donors$n_cells_qc_pass >= min.cells.per.donor
  ]

  dat.counts <- read.delim(phenotype.bed.file, row.names = NULL, quote = "")
  row.names(dat.counts) <- dat.counts[,COLUMN_GENE_ID]


  n_cols = dim(dat.counts)[2]
  dm <- as.matrix(dat.counts[,(COLUMN_GENE_ID+1):n_cols])
  ## subset matrix to remove donors which fewer than MIN_CELLS_PER_DONOR cells
  dm <- dm[,as.vector(sample.names.included)]
  ## subset matrix to remove unobserved genes
  sv.nonzero.genes = rowSums(dm) != 0
  dm <- dm[sv.nonzero.genes,]

  # TMM normalization
  dge <- edgeR::DGEList(counts=dm)
  dge <- edgeR::calcNormFactors(dge)

  # defaults:
  # calcNormFactors(
  #  dm,
  #  lib.size = None, method = "TMM", refColumn = NULL,
  #  logratioTrim = 0.3, sumTrim = 0.5, doWeighting = TRUE,
  #  Acutoff = -1e10
  #  )

  # log2(CPM) as output
  m.logcpm <- edgeR::cpm(dge, log=TRUE)

  df <- cbind(dat.counts[sv.nonzero.genes,1:COLUMN_GENE_ID], as.data.frame(m.logcpm))
}

## main
## infnam.bed <- "aggrsum_all_counts.bed.gz"
## infnam.tsv <- "aggrsum_all_ncells.tsv"
args.v <- commandArgs(trailingOnly = TRUE)
cat("invoked R like this:\n");
cat(paste(args.v, collapse=" "), "\n")
args.n = length(args.v)
if (args.n > 1 & args.n < 6) {
  infnam.bed <- args.v[1]
  infnam.tsv <- args.v[2]
  n.cells.min <- DEFAULT_MIN_CELLS_PER_DONOR
  n.pcs <- 0
  output_is_perChr <- FALSE
  if (args.n > 2) {
    n.cells.min <- args.v[3]
    if (args.n > 3) {
      n.pcs <- args.v[4]
      if (args.n > 4) {
        output_is_perChr <- args.v[5] == 'perChr'
      }
    }
  }
} else {
  stop("Call with <input BED file> <input TSV table> <min_cell_mum> <n_PCs> [perChr]\n", call.=FALSE)
}

cat ("bed file:", infnam.bed, ", tsv table:", infnam.tsv, "n.cells.min =", n.cells.min, "n.pcs =", n.pcs, "\n")
dirnam <- paste0("norm", sub("\\.bed.gz$", basename(infnam.bed), replacement=""))
if (output_is_perChr && !file.exists(dirnam)) {
  dir.create(dirnam)
}

df <- normalise.dSumTMM(
  phenotype.bed.file = infnam.bed,
  donor.tsv.file = infnam.tsv,
  min.cells.per.donor = n.cells.min
  )
if (n.pcs > 0) {
  cat("write", n.pcs, "principal components ...\n")
  pca <- write.principal.components(df, dirnam, n.pcs = n.pcs)
}

if (output_is_perChr) {
  write.perChrBED(df, dirnam)
} else {
  write.singleBED(df, oufn = paste0(dirnam, "_chrAll.bed"))
}
cat("finished.")
