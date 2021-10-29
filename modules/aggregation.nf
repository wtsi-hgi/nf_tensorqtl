#!/usr/bin/env nextflow

process aggregate_UMI_counts_total_sum {
  input:
    path(donor_file_dir) // lists input files per donor
    path(gene_annotation_file)
    val(n_cells_min)
    val(n_donors_min)
    val(cis_window_pos)

  output:
    path("aggrsum_files.csv", emit:file_list_csv)
    path("aggrsum_output", emit:output_dir)

  script:
  """
    aggrsce.py --handover-data-dir ${donor_file_dir} \
      --gene-annotations ${gene_annotation_file} \
      --output-dir aggrsum_output \
      --output-file-prefix aggrsum \
      --output-file-list aggrsum_files.csv \
      --minium-cell-number ${n_cells_min} \
      --minimum-donor-number ${n_donors_min} \
      --bed-cis-position-in-gene ${cis_window_pos}
  """
}

process normalize_counts_TMM {
  input:
    tuple val(celltype), val(counts_bed), val(donor_tsv)
    path(data_dir)
    val(n_cell_min)
    val(n_prcmp)

  output:
    path("${outprfx}_chrAll.bed.gz", emit:norm_bed)
    path("${outprfx}_PC${npcs}.tsv", emit:pc_tsv) optional true

  script:
    outprfx = "norm" + "${counts_bed}".minus(".bed.gz")

    //only calculate PCs for pseudo-bulk sample
    if ("${celltype}" == "all") {
      npcs = "${n_prcmp}"
    } else {
      npcs = "0" // signals not to attempt PC caculations
    }
  """
    dSumTMM.R ${data_dir}/${counts_bed} ${data_dir}/${donor_tsv} ${n_cell_min} ${npcs}
    gzip ${outprfx}_chrAll.bed
  """
}

workflow aggregate_normalize_dSum {
  take:
    donor_file_dir
    gene_annotation_file
    n_cell_min
    n_donor_min
    cis_window_pos
    n_expression_pcs

  main:
    aggregate_UMI_counts_total_sum (
      donor_file_dir,
      gene_annotation_file,
      n_cell_min,
      n_donor_min,
      cis_window_pos
    )
    aggregate_UMI_counts_total_sum.out.file_list_csv
      .splitCsv(header: true)
      .view()
      .set { ch_celltype_files }

    normalize_counts_TMM(
      ch_celltype_files,
      aggregate_UMI_counts_total_sum.out.output_dir,
      n_cell_min,
      n_expression_pcs
    )
  emit:
    aggrnorm_bed = normalize_counts_TMM.out.norm_bed
    expression_pcs_tsv = normalize_counts_TMM.out.pc_tsv
}
