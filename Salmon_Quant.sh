#!/bin/bash
#SBATCH --time=18:00:00
#SBATCH --account=def-dogfish
#SBATCH --mem=20000M
#SBATCH --mail-user=edward32@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --array=1-57

module load StdEnv/2023  gcc/12.3  openmpi/4.1.5
module load salmon/1.10.2

raw_read_list=/home/edward32/scratch/sturgeon_OE/raw_reads_paired.txt
string="sed -n "$SLURM_ARRAY_TASK_ID"p ${raw_read_list}"
str=$($string)

var=$(echo $str | awk -F"\t" '{print $1, $2}')
set -- $var

forward=$1
reverse=$2

echo ${forward}
echo ${reverse}

sample=$(basename ${forward} _R1.fastq.gz | cut -d "." -f5)

echo ${sample}

salmon quant -i /home/edward32/scratch/sturgeon_OE/Transcriptome/salmon_index -l IU -1 /home/edward32/scratch/sturgeon_OE/trimmed_data/${sample}.R1.fq.gz -2 /home/edward32/scratch/sturgeon_OE/trimmed_data/${sample}.R2.fq.gz --validateMappings --seqBias --gcBias -o /home/edward32/scratch/sturgeon_OE/sturgeon_quant/${sample}
