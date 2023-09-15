#!/bin/bash

#$ -cwd
#$ -V
#$ -N LIGYSIS_arr
#$ -M 2394007@dundee.ac.uk
# -m a
#$ -t 1-5
#$ -tc 1
#$ -o logs/ligysis_pdb/o
#$ -e logs/ligysis_pdb/e

# Lookup task-array and cd to job directory to avoid possible tempfile collisions
UP_ACC=$(sed "${SGE_TASK_ID}q;d" ./input/up_ids.txt)

dest=./output/${UP_ACC}

# Check if the directory exists; if not, create it
if [ ! -d "${dest}" ]; then
    mkdir -p "${dest}"
fi

cd ${dest}

# Run ligysis for each protein within its own directory
python3.6 ./scripts/ligysis.py --transform --experimental --variants --override --up_acc ${UP_ACC}
