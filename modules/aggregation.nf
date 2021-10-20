#!/usr/bin/env nextflow

process aggregate_UMI_counts_dSUM {
  input:
    path(donor_file_dir) // lists input files per donor
    path(gene_annotation_file)
    val(minimum_number_of_cells)

  output:
    path("count_*.bed.gz", emit: bed_files)
    path("ncells_*.tsv", emit:tsv_files)

  script:
  """
    aggrsce.py --handover-data-dir ${donor_file_dir} --gene-annotations ${gene_annotation_file} --min-cell-num ${minimum_number_of_cells}
  """
}

workflow dSUM_aggregation {
  take:
    donor_file_dir
    gene_annotation_file
    minimum_number_of_cells

  main:
    aggregate_UMI_counts_dSUM (
      donor_file_dir,
      gene_annotation_file,
      minimum_number_of_cells
    )
    aggregate_UMI_counts_dSUM.out.bed_files
    .flatMap()
    .toSortedList()
    .set{ch_bed_files}

    aggregate_UMI_counts_dSUM.out.tsv_files
    .flatMap()
    .toSortedList()
    .set{ch_tsv_files}

    //ch_bed_files.merge(ch_tsv_files)
    //.subscribe onNext {  println "ch_bed_files: $it"  }, onComplete { println "Done."}

    emit:
      bed_and_tsv_files = ch_bed_files
}
