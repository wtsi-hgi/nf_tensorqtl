#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

VERSION = "0.0.1" // Do not edit, controlled by bumpversion.

include {
    dSUM_aggregation;
} from "./modules/aggregation.nf"

// default parameters

params.input_dir = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/handover/franke_data4"
// directory of input files with a file 'files.tsv' listing donors and input files
// ----
// experiment_id	donor_id	filename_h5ad	filename_annotation_tsv
// franke_Pilot_3_lane_1	s3	franke_Pilot_3_lane_1.donor0.h5ad	franke_Pilot_3_lane_1.donor0.tsv
// franke_Pilot_3_lane_1	s1	franke_Pilot_3_lane_1.donor1.h5ad	franke_Pilot_3_lane_1.donor1.tsv
// franke_Pilot_3_lane_1	s5	franke_Pilot_3_lane_1.donor2.h5ad	franke_Pilot_3_lane_1.donor2.tsv

// params.ensembl_url = "ensembldb-mirror.internal.sanger.ac.uk"
// needed for automated download/generation of gene annotation file

params.gene_annotation = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/handover/franke_data4/gene_annot.tsv"
// from EnsEMBL
// prepared using .bin/lookupEnsEMBLGeneIds.pl
// ----
// query_id	ensembl_id	display_id	chromosome	start	end	strand	hgnc_symbols
// ENSG00000243485	ENSG00000243485	ENSG00000243485	1	29554	31109	1	MIR1302-2HG
// ENSG00000237613	ENSG00000237613	ENSG00000237613	1	34554	36081	-1	FAM138A	F379
// ENSG00000186092	ENSG00000186092	ENSG00000186092	1	65419	71585	1	OR4F5
// ...

params.minimum_cell_number = 5
// minimum number of cells of a given cell-type donors must have to be included in eQTL analysis

log.info """
============================================================================
  sceQTL analysis using linear models with tensorQTL ~ v${VERSION}
============================================================================
input directory                                : ${params.input_dir}
gene annotation file                           : ${params.gene_annotation}
minimum number of cells per donor              : ${params.minimum_cell_number}
""".stripIndent()

workflow {
    main:
    log.info "Running main workflow."
    log.info """\n
      test
    """.stripIndent()

    dSUM_aggregation(
      params.input_dir,
      params.gene_annotation,
      params.minimum_cell_number
    )
    dSUM_aggregation.out.bed_and_tsv_files
    .view()
}

workflow.onComplete {
    log.info """\n
    ----------------------------------------------------------------------------
     pipeline execution summary
    ----------------------------------------------------------------------------
    Completed         : ${workflow.complete}
    Duration          : ${workflow.duration}
    Success           : ${workflow.success}
    Work directory    : ${workflow.workDir}
    Exit status       : ${workflow.exitStatus}
    """.stripIndent()
}