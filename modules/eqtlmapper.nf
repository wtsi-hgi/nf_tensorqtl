#!/usr/bin/env nextflow

process run_tensorqtl {
  label 'gpu'

  input:
    path(aggrnorm_counts_bed)
    path(genotype_pcs_tsv)
    path(plink_files_dir)
    val(plink_files_prefix)

  """
  oufnprfx = "${aggrnorm_counts_bed}".minus("_counts_chrAll.bed.gz")
  mapqtl.py \
    --gene-expression-bed ${aggrnorm_counts_bed} \
    --genotype-plink-prefix ${plink_files_dir}/${plink_files_prefix} \
    --genotype-pcs-plink ${genotype_pcs_tsv}
    --output-prefix mapqtl_${oufnprfx}
  """
}

workflow map_eqtl {
  take:
    aggrnormcounts_bedfile
    genotpcs_tsvfile
    plink_dir
    plink_prefix

  main:
    run_tensorqtl(
      aggrnormcounts_bedfile,
      genotpcs_tsvfile,
      plink_dir,
      plink_prefix
    )
}
