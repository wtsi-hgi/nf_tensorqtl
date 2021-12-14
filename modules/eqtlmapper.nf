#!/usr/bin/env nextflow

process prep_plink_files_from_genotype_vcf {
  input:
    path(genotype_vcf)

  output:
    path "${plink_bed_fnprfx}.*", emit: plink_bed_files

  script:
    plink_bed_fnprfx = "${genotype_vcf}".minus('.vcf.gz').plus('_plink_bed')
    // plink2 --memory argument interpreted as MiB
    """
      plink2 --make-pgen --max-alleles 2 --vcf ${genotype_vcf} --out tmp_gt_plink
      plink2 --make-pgen --sort-vars --threads ${task.cpus} --pfile tmp_gt_plink --out tmp_gt_plink_srt
      plink2 --make-bed --output-chr chrM --pfile tmp_gt_plink_srt --snps-only --out ${plink_bed_fnprfx}
    """
}

process calc_genotype_pca_plink {
  input:
    path plink_bed_files

  output:
    path "${plink_gtpca_fnprfx}.eigenvec", emit: plink_gtpca_eigenvec_tsv

  script:
    plink_bed_fnprfx = plink_bed_files[0].getBaseName()
    plink_gtpca_fnprfx = plink_bed_fnprfx.minus('_bed').plus('_gtpca')
    """
    plink2 --freq counts --bfile ${plink_bed_fnprfx} --out tmp_gt_plink_freq
    plink2 --pca --read-freq tmp_gt_plink_freq.acount --bfile ${plink_bed_fnprfx} --out ${plink_gtpca_fnprfx}
    """
}

process run_tensorqtl {
  label 'gpu'

  publishDir  path: "${outdir}",
              overwrite: "true"

  input:
    path(aggrnorm_counts_bed)
    path(gene_annotation_tsv)
    path(genotype_pcs_tsv)
    path(plink_files)

  output:
    path("mapqtl_${oufnprfx}.cis_eqtl.tsv.gz"), emit: qtl_tsv
    path("mapqtl_${oufnprfx}.cis_eqtl_dropped.tsv.gz"), emit: dropped_tsv
    path("mapqtl_${oufnprfx}.cis_eqtl_qval.tsv.gz"), emit: qval_tsv
    //path("mapqtl_${oufnprfx}.bin.gz"), emit: qtl_bin

  script:
  outdir = "${launchDir}/" + params.output_dir
  oufnprfx = "${aggrnorm_counts_bed}".minus("_counts_chrAll.bed.gz")
  plink_files_prefix = plink_files[0].getBaseName()
  """
  mapqtl.py \
    --gene-expression-bed ${aggrnorm_counts_bed} \
    --gene-annotation ${gene_annotation_tsv} \
    --genotype-plink-prefix ${plink_files_prefix} \
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

process plot_rank_comparison {

  publishDir  path: "${outdir}",
              overwrite: "true"

  input:
    path(cis_qtl_qval_tsv)
    path(cis_qtl_ref_csv)

  output:
    path("${plotfnout}*", emit: plotfiles)

  script:
    outdir = "${launchDir}/" + params.output_dir
    plotfnout = "${cis_qtl_qval_tsv}".minus(".tsv.gz").plus("plot")
    """
    plotQval.py ${cis_qtl_ref_csv} ${cis_qtl_qval_tsv} ${plotfnout}
    """
}

workflow prep_genotypes {
  take:
    genotype_vcf

  main:
    prep_plink_files_from_genotype_vcf(
      genotype_vcf
      )
    calc_genotype_pca_plink(
      prep_plink_files_from_genotype_vcf.out.plink_bed_files
      )
  emit:
    plink_genotype_pcs_tsv = calc_genotype_pca_plink.out.plink_gtpca_eigenvec_tsv
    plink_files = prep_plink_files_from_genotype_vcf.out.plink_bed_files
}

workflow map_eqtl {
  take:
    aggrnormcounts_bedfile
    gene_annotation_tsvfile
    genotpcs_tsvfile
    plink_files
    eqtl_reference_csvfile

  main:
    run_tensorqtl(
      aggrnormcounts_bedfile,
      gene_annotation_tsvfile,
      genotpcs_tsvfile,
      plink_files
    )
    add_qvalues(
      run_tensorqtl.out.qval_tsv
    )
    add_qvalues.out.cis_qtl_qval_tsv
      .filter( ~/\S+_pseudobulk.\S+/ )
      .set {ch_pseudobulk_eqtls}

    plot_rank_comparison(
      ch_pseudobulk_eqtls,
      eqtl_reference_csvfile
    )
  emit:
    cis_qtl_tsvfile = add_qvalues.out.cis_qtl_qval_tsv
}
