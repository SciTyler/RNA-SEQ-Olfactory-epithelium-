#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --account=def-dogfish
#SBATCH --mem=20000M
#SBATCH --mail-user=edward32@myumanitoba.ca
#SBATCH --mail-type=ALL

module load StdEnv/2023  gcc/12.3  openmpi/4.1.5
module load salmon/1.10.2

salmon index -t /home/edward32/scratch/sturgeon_OE/Transcriptome/trinity_outputs.Trinity.fasta -i /home/edward32/scratch/sturgeon_OE/Transcriptome/salmon_index -k 31
