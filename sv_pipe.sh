#!/bin/bash

#SBATCH --job-name=
#SBATCH --time=05:00:00
#SBATCH --account=marth-rw
#SBATCH --partition=marth-shared-rw
#SBATCH -o {path}/slurm_out/%j.out
#SBATCH -e {path}/slurm_out/%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=

#Inputs to the script
PBDCRAM=$1
MOMCRAM=$2
DADCRAM=$3
OG_SMOOVE_VCF=$4
REF_FASTA=$5

DOCTORED_MANTA_OUTPUT="doc_manta_vcf"
DUPHOLD_MANTA_OUTPUT="dhanno_manta_vcf"
SVAF_SMOOVE_OUTPUT="svaf_smoove_vcf"
SVAF_MANTA_OUTPUT="svaf_manta_vcf"

# set up scratch directory
SCRDIR=/scratch/ucgd/lustre-labs/marth/scratch/$USER/$SLURM_JOB_ID
mkdir -p $SCRDIR
# copy the scripts and reference files to the scratch dir
cp \
    run_manta_trio.sh \
    doctor_manta.py \
    ./smoove/bp_smoove.sif \
    sv_pipe.yml \
    ./ref_files/recommended_exclusions_37.bed \
    ./ref_files/recommended_exclusions_38.bed \
    SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz \
    $SCRDIR

cd $SCRDIR

# run the the manta script run_manta_trio.sh (needs the trio's urls CRAM format) -> new joint called vcf
./run_manta_trio.sh $PBDCRAM $MOMCRAM $DADCRAM 

# load the miniconda3 module will be needed for the doctor_manta.py and for the svafotate run at the end
module load \
    miniconda3/23.11.0 \
    singularity/4.1.1

# Check if the environment already exists, if not create it from the yml
ENV_NAME="sv_pipe_conda_env"
if conda info --envs | grep -q "$ENV_NAME"; then
    conda activate "$ENV_NAME"
else
    conda env create -f sv_pipe.yml
    conda activate "$ENV_NAME"
fi

# the script is going to put a simlink to the manta vcf at ./diploidSV.vcf.gz
# modify the vcf with his doctor_manta.py script (needs input vcf, and output) -> new doc_vcf
python doctor_manta.py diploidSV.vcf.gz $DOCTORED_MANTA_OUTPUT
conda deactivate "$ENV_NAME"

# TODO: run smoove duphold on the manta vcf -> dhanno_manta_vcf
singularity exec bp_smoove.sif duphold -v $DOCTORED_MANTA_OUTPUT

#Run svafotate on both files (.8 ol threshold) -> filtered svaf_vcf
#we will either create a conda env or activate it
svafotate annotate -v $DUPHOLD_MANTA_OUTPUT -b SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz -o $SVAF_MANTA_OUTPUT -f 0.8
svafotate annotate -v $OG_SMOOVE_VCF -b SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz -o $SVAF_SMOOVE_OUTPUT -f 0.8
