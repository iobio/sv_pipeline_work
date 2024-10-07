#!/bin/bash

#SBATCH --job-name=
#SBATCH --time=04:00:00
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

# set up scratch dir in scratch
SCRDIR=/scratch/ucgd/lustre-labs/marth/scratch/$USER/$SLURM_JOB_ID
mkdir -p $SCRDIR

# copy the scripts (not going to copy the data files themselves)
cp run_manta_trio.sh doctor_manta.py sv_pipe_env.yml $SCRDIR

# cd into that dir
cd $SCRDIR

# un the the manta script run_manta_trio.sh (needs the trio's urls CRAM format) -> new joint called vcf
./run_manta_trio.sh $PBDCRAM $MOMCRAM $DADCRAM 

# load the miniconda3 module will be needed for the doctor_manta.py and for the svafotate run at the end
module load miniconda3/23.11.0

ENV_NAME="sv_pipe_conda_env"

# Check if the environment already exists
if conda info --envs | grep -q "$ENV_NAME"; then
    conda activate "$ENV_NAME"
else
    conda env create -f sv_pipe_env.yml
    conda activate "$ENV_NAME"
fi

#the script is going to put a simlink to the manta vcf at ./diploidSV.vcf.gz
#modify the vcf with his doctor_manta.py script (needs input vcf, and output) -> new doc_vcf

#need to have a conda environment that will have 

#point to the singularity the singularity file will be in that same place
#run DupHold on the manta vcf -> dhanno_vcf
#per tom: user the smoove duphold module (will need to start singularity pointing to smoove docker container)

#now we have a manta with duphold annotations, smoove we will already have

#Run svafotate on both files (.8 ol threshold) -> filtered svaf_vcf
#we will either create a conda env with smoove or activate it
