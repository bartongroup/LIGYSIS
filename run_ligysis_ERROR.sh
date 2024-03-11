#!/bin/bash

#$ -cwd
#$ -V
#$ -N LIGYSIS_BioLiP_V2
#$ -M 2394007@dundee.ac.uk
#$ -m a
#$ -t 1-479
#$ -tc 5
#$ -o logs/V2/r1/o
#$ -e logs/V2/r1/e
# -jc long
#$ -mods l_hard mfree 64G
# -adds l_hard hostname '!(m910-*|c6*|gpu-33*) #|m910-6*)

# Look up UP_ACC and original_task_id
UP_ACC=$(awk "NR==$SGE_TASK_ID{print \$1}" LIGYSIS_V3_unproc_accs_1.txt)
og_task_id=$(awk "NR==$SGE_TASK_ID{print \$2}" LIGYSIS_V3_unproc_accs_1.txt)

dest=./output_V2/${UP_ACC}

# Check if the directory exists; if not, create it
if [ ! -d "${dest}" ]; then
    mkdir -p "${dest}"
fi

cd ${dest}

#echo $(pwd)

python3.6 ./../../ligysis.py ${UP_ACC} --override --override_trans

#python3.6 ./../../ligysis.py ${UP_ACC}

# Rename output and error files to include original_task_id
cd ../..
mv logs/V2/r1/e/LIGYSIS_BioLiP_V2.e${JOB_ID}.${SGE_TASK_ID} logs/V2/r1/e/LIGYSIS_BioLiP_V2.e${JOB_ID}.${og_task_id}
mv logs/V2/r1/o/LIGYSIS_BioLiP_V2.o${JOB_ID}.${SGE_TASK_ID} logs/V2/r1/o/LIGYSIS_BioLiP_V2.o${JOB_ID}.${og_task_id}
