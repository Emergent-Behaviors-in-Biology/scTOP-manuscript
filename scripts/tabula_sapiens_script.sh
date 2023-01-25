#!/bin/bash -l
module load python3/3.8.10

# Set SCC project
#$ -P biophys

# Request at least 256 GB of ram
#$ -pe omp 16
#$ -l mem_per_core=16G

# Load packages
source /projectnb/biophys/mariay/venvs/myenv/bin/activate
export PYTHONPATH=/projectnb/biophys/mariay/venvs/myenv/lib/python3.8/site-packages/:$PYTHONPATH
export PATH=/projectnb/biophys/mariay/venvs/myenv/bin:$PATH

python /project/biophys/mouse\ scRNA-seq/scTOP\ explorations/scTOP/tutorial/tabula_sapiens_accuracies.py
