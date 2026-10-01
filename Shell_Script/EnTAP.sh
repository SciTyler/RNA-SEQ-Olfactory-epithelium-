#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-dogfish
#SBATCH --mem=0
#SBATCH --mail-user=edward32@myumanitoba.ca
#SBATCH --mail-type=ALL

module load interproscan_data/5.64-96.0 (needs stdenv2020)
module load sqlite/3.42.0
module load StdEnv/2023 
module load diamond/2.1.8
module load python/3.10.13
pip install --user pysqlite3 biopython tqdm numpy

export PATH=/home/edward32/scratch/sturgeon_OE/EnTAP/EnTAP-v2.0.0/libs/eggnog-mapper-2.1.12/build/lib/eggnogmapper/bin:$PATH

#configuration
EnTAP --config --run-ini /home/edward32/scratch/sturgeon_OE/EnTAP/EnTAP-v2.0.0/entap_run.params --entap-ini /home/edward32/scratch/sturgeon_OE/EnTAP/EnTAP-v2.0.0/entap_config.ini

#running 
EnTAP --runP --run-ini /home/edward32/scratch/sturgeon_OE/EnTAP/EnTAP-v2.0.0/entap_run.params --entap-ini /home/edward32/scratch/sturgeon_OE/EnTAP/EnTAP-v2.0.0/entap_config.ini
