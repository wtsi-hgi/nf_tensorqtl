#!/usr/bin/env nextflow

process aggregate_UMI_counts_total_sum {
  input:
    path(donor_file_dir) // lists input files per donor
    path(gene_annotation_file)
    val(n_cells_min)
    val(n_donors_min)

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
      --minimum-donor-number ${n_donors_min}
  """
}

process normalize_counts_TMM {
  input:
    tuple val(celltype), val(counts_bed), val(donor_tsv)
    path(data_dir)
    val(n_cell_min)

  output:
    path("${outfn}.gz", emit:norm_bed)

  script:
    outfn = "norm" + "${counts_bed}".minus(".bed.gz") + "_chrAll.bed"
  """
    dSumTMM.R ${data_dir}/${counts_bed} ${data_dir}/${donor_tsv} ${n_cell_min}
    gzip ${outfn}
  """
}

workflow dSUM_aggregation {
  take:
    donor_file_dir
    gene_annotation_file
    n_cell_min
    n_cell_donor

  main:
    aggregate_UMI_counts_total_sum (
      donor_file_dir,
      gene_annotation_file,
      n_cell_min,
      n_cell_donor
    )
    aggregate_UMI_counts_total_sum.out.file_list_csv
      .splitCsv(header: true)
      .view()
      .set { ch_celltype_files }

    normalize_counts_TMM(
      ch_celltype_files,
      aggregate_UMI_counts_total_sum.out.output_dir,
      n_cell_min
    )
  emit:
    aggrnorm_bed = normalize_counts_TMM.out.norm_bed
}
