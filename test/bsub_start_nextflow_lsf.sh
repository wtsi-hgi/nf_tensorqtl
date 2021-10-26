#!/bin/sh

mem=2400
queue=yesterday

# clean up previous run files
rm -f *.log
rm -f bsub.o
rm -f bsub.e
rm -f bjob.id

# OBS: nextflow expects the bin directory for software and custom scripts in projectDir (directory of the main script)
# set a link if not yet there

if [ ! -L ./test/bin ]; then
  echo "Nextflow expects the bin directory for custom scripts in projectDir (directory of the main script)."
  echo "Setting a symlink in ./test/bin."
  cd ./test
  ln -s ../bin ./bin
  cd ..
fi

# start Nextflow via bsub:
bsub -G team151 \
     -R"select[mem>${mem}] rusage[mem=${mem}]" \
     -M ${mem} \
     -o bsub.o -e bsub.e \
     -q ${queue} \
     bash ./test/start_nextflow_lsf.sh > bjob.id

# get process PID
echo "Nextflow Bjob ID saved in file bjob.id"
echo kill with \"bkill ID_number\" command
echo "check logs files bsub.o, bsub.e and .nextflow.log"
