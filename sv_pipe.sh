#!/bin/bash

#SBATCH --job-name=
#SBATCH --time=03:00:00
#SBATCH --account=marth-rw
#SBATCH --partition=marth-shared-rw
#SBATCH -o {path}/slurm_out/%j.out
#SBATCH -e {path}/slurm_out/%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=

#CRAM Paths
PBDCRAM=$1
MOMCRAM=$2
DADCRAM=$3
DOC_VCF_OUTPUT_FILE=$4
REF_FASTA=$5

#run the the manta script run_manta_trio.sh (needs the trio's urls CRAM format) -> new joint called vcf

#modify the vcf with his doctor_manta.py script (needs input vcf, and output) -> new doc_vcf

#run DupHold on the manta vcf -> dhanno_vcf
#per tom: user the smoove duphold module (will need to start singularity pointing to smoove docker container)

#now we have a manta with duphold annotations, smoove we will already have

#Run svafotate on both files (.8 ol threshold) -> filtered svaf_vcf
#we will either create a conda env with smoove or activate it
