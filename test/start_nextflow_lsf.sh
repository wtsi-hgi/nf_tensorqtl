#!/bin/sh

# activate Nextflow conda env
conda init bash
eval "$(conda shell.bash hook)"
conda activate nextflow

# run nextflow main.nf with inputs and lsf config:
export NXF_OPTS="-Xms2G -Xmx2G"
#nextflow run ./nextflow_ci/pipelines/main.nf \

nextflow run ./test/test_channel.nf \
  -c ./nextflow.config \
  -profile lsf
