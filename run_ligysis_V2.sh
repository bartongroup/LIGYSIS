#!/bin/bash

#$ -cwd
#$ -V
#$ -N LIGYSIS_BioLiP_V2
#$ -M 2394007@dundee.ac.uk
#$ -m a
#$ -t 1501-2000
#$ -tc 10
#$ -o logs/V2/r1/o
#$ -e logs/V2/r1/e
# -jc long
#$ -mods l_hard mfree 16G
# -adds l_hard hostname '!(m910-*|c6*) #|m910-6*)


# Lookup task-array and cd to job directory to avoid possible tempfile collisions
UP_ACC=$(sed "${SGE_TASK_ID}q;d" ./input/biolip_up_ids_filt.txt)

dest=./output_V2/${UP_ACC}

# Check if the directory exists; if not, create it
if [ ! -d "${dest}" ]; then
    mkdir -p "${dest}"
fi

cd ${dest}

# Run ligysis for each protein within its own directory
python3.6 ./../../ligysis.py ${UP_ACC} --override
