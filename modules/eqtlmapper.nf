#!/usr/bin/env nextflow

process run_tensorqtl {
  label 'gpu'

  publishDir  path: "${outdir}",
              overwrite: "true"

  input:
    path(aggrnorm_counts_bed)
    path(gene_annotation_tsv)
    path(genotype_pcs_tsv)
    path(plink_files_dir)
    val(plink_files_prefix)

  output:
    path("mapqtl_${oufnprfx}.cis_eqtl.tsv.gz"), emit: qtl_tsv
    path("mapqtl_${oufnprfx}.cis_eqtl_dropped.tsv.gz"), emit: dropped_tsv
    path("mapqtl_${oufnprfx}.cis_eqtl_qval.tsv.gz"), emit: qval_tsv
    //path("mapqtl_${oufnprfx}.bin.gz"), emit: qtl_bin

  script:
  outdir = "${launchDir}/" + params.output_dir
  oufnprfx = "${aggrnorm_counts_bed}".minus("_counts_chrAll.bed.gz")
  """
  mapqtl.py \
    --gene-expression-bed ${aggrnorm_counts_bed} \
    --gene-annotation ${gene_annotation_tsv} \
    --genotype-plink-prefix ${plink_files_dir}/${plink_files_prefix} \
    --genotype-pcs-plink ${genotype_pcs_tsv} \
    --output-prefix mapqtl_${oufnprfx}
  """
}

process add_qvalues {
  publishDir  path: "${outdir}",
              overwrite: "true"

  input:
    path(cis_qtl_tsv)

  output:
    path("${oufnam}", emit:cis_qtl_qval_tsv)

  script:
    outdir = "${launchDir}/" + params.output_dir
    oufnam = "${cis_qtl_tsv}".minus(".tsv.gz") + "_qvalR.tsv.gz"
  """
  qval.R ${cis_qtl_tsv} ${oufnam}
  """
}

workflow map_eqtl {
  take:
    aggrnormcounts_bedfile
    gene_annotation_tsvfile
    genotpcs_tsvfile
    plink_dir
    plink_prefix

  main:
    run_tensorqtl(
      aggrnormcounts_bedfile,
      gene_annotation_tsvfile,
      genotpcs_tsvfile,
      plink_dir,
      plink_prefix
    )
    add_qvalues(
      run_tensorqtl.out.qval_tsv
    )
  emit:
    cis_qtl_tsvfile = add_qvalues.out.cis_qtl_qval_tsv
}
