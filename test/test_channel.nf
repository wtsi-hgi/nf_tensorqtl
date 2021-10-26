#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//include {
//    file_list_to_channel;
//} from "../modules/aggregation.nf"

datadir = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/handover/franke_data4/"
ch_input_file_table = Channel.fromPath(
  datadir + "aggr_output_files.csv"
  )
ch_data_dir = Channel.fromPath(datadir + "sumaggr_out")

process dummy {
  input:
    tuple val(celltype), val(bedfile), val(tsvfile)
    path(data_dir)

  script:
    """
    echo "data_dir = ${data_dir}"
    echo "celltype = ${celltype}"
    echo "bedfile = ${bedfile}"
    """
}

workflow {
    main:
      ch_input_file_table
      .splitCsv(header: true)
      .view()
      .set { my_channel }

//      my_channel.subscribe onNext: {println it}, onComplete: {println 'Done.'}

    dummy(
     my_channel,
     ch_data_dir
   )
}
