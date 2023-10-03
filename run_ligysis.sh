#!/bin/bash

#$ -cwd
#$ -V
#$ -N LIGYSIS_BioLiP
#$ -M 2394007@dundee.ac.uk
#$ -m a
#$ -t 1-100
#$ -tc 1
#$ -o logs/r2/o
#$ -e logs/r2/e

# Lookup task-array and cd to job directory to avoid possible tempfile collisions
UP_ACC=$(sed "${SGE_TASK_ID}q;d" ./input/biolip_up_ids_filt.txt)

dest=./output/${UP_ACC}

# Check if the directory exists; if not, create it
if [ ! -d "${dest}" ]; then
    mkdir -p "${dest}"
fi

cd ${dest}

# Run ligysis for each protein within its own directory
python3.6 ./../../ligysis.py ${UP_ACC} --override
