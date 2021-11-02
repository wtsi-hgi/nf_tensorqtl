#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

VERSION = "0.0.1" // Do not edit, controlled by bumpversion.

include {
    aggregate_normalize_dSum;
} from "./modules/aggregation.nf"

include {
  map_eqtl;
} from "./modules/eqtlmapper.nf"

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

params.gene_annotation = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/handover/franke_data4/gene_annot_chr.tsv"
// from EnsEMBL
// prepared using .bin/lookupEnsEMBLGeneIds.pl
// IMPORTANT: chromosome name need 'chr' prefix as in PLINK files
// ----
// query_id	ensembl_id	display_id	chromosome	start	end	strand	hgnc_symbols
// ENSG00000243485	ENSG00000243485	ENSG00000243485	chr1	29554	31109	1	MIR1302-2HG
// ENSG00000237613	ENSG00000237613	ENSG00000237613	chr1	34554	36081	-1	FAM138A	F379
// ENSG00000186092	ENSG00000186092	ENSG00000186092	chr1	65419	71585	1	OR4F5
// ...

params.plink_dir = "/lustre/scratch123/hgi/projects/ukbb_scrna/eval/groningen/"
params.plink_filnam_prefix = "franke_gt_plink_bed"
params.genotype_pcs_tsvfile = "/lustre/scratch123/hgi/projects/ukbb_scrna/eval/groningen/franke_gtpca.tsv"

params.minimum_cell_number = 5
// minimum number of cells of a given cell-type donors must have to be included in eQTL analysis

params.minimum_donor_number = 30
// minimum number of donors, with at least params.minimum_cell_number each, a cell type must have to be included in eQTL analysis

params.cis_window_pos = "mid"
// position of CIS window
//'mid': middle of annotated gene start and end point
//'start': start of annotated gene taking strand into account (i.e. approx. transcription start site)

params.expression_prcmp_num = 10
// number of principal components to use for expression vectors

params.output_dir = "results"
// output directory relative to ${launchDir}
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

    aggregate_normalize_dSum(
      params.input_dir,
      params.gene_annotation,
      params.minimum_cell_number,
      params.minimum_donor_number,
      params.cis_window_pos,
      params.expression_prcmp_num
    )
    //aggregate_normalize_dSum.out.aggrnorm_bed
    //.view()
    aggregate_normalize_dSum.out.expression_pcs_tsv
    .view()

    map_eqtl(
      aggregate_normalize_dSum.out.aggrnorm_bed,
      params.genotype_pcs_tsvfile,
      params.plink_dir,
      params.plink_filnam_prefix
    )
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
