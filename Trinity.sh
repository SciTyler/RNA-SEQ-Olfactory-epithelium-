#!/bin/bash
#SBATCH --time=7-00:00:00
#SBATCH --account=def-dogfish
#SBATCH --mem=0
#SBATCH --mail-user=edward32@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --propagate=STACK

module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3
module load trinity/2.14.0
module load samtools/1.17
module load jellyfish/2.3.0 salmon/1.4.0
module load scipy-stack
pip list | grep -a numpy

Trinity --seqType fq --max_memory 160G --CPU 8 --min_kmer_cov 2 --min_contig_length 300 --normalize_by_read_set --verbose  --samples_file /home/edward32/scratch/sturgeon_OE/trinity_samples.txt \
--output /home/edward32/scratch/sturgeon_OE/trinity_outputs
