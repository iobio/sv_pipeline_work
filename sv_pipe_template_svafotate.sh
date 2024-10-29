#!/bin/bash

# THIS SCRIPT IS MINIMIZED FOR 

#SBATCH --job-name=<any_name_you_want>
#SBATCH --time=10:00:00

# Your Account and Partition as Applicable
#SBATCH --account=marth-rw
#SBATCH --partition=marth-shared-rw

# Output in your directory in the scratch space
#SBATCH -o /scratch/ucgd/lustre-labs/marth/scratch/<your_uid>/slurm-%j.out-%N
#SBATCH -e /scratch/ucgd/lustre-labs/marth/scratch/<your_uid>/slurm-%j.err-%N

# Since this isnt calling manta cutting down the mem and cpu usage
#SBATCH --ntasks=10
#SBATCH --mem=60G
#SBATCH --mail-type=ALL

# The email to send system messages to
#SBATCH --mail-user=<your_email>

#Inputs to the script 
# THIS YML ONLY NEEDS THE 'probandId' AND THE 'smoove' FILE PATH
YMLDEF=$1

# Grab the neccessary items from the YMLDEF passed via the command
PROBANDID=$(grep 'probandId:' "$YMLDEF" | awk '{print $2}')

O_SMOOVE_VCF=$(grep 'smoove:' "$YMLDEF" | awk '{print $2}')

REF_FASTA='human_g1k_v38_decoy_phix.fasta'
FASTA_FOLDER='/scratch/ucgd/lustre/common/data/Reference/GRCh38'

SVAF_SMOOVE_OUTPUT="svaf_smoove.vcf"
DSMOOVEFILE=deno_smoove.vcf.gz

# Set up scratch directory
FOLDERNAME=${PROBANDID}_${SLURM_JOB_ID}
SCRDIR=/scratch/ucgd/lustre-labs/marth/scratch/$USER/$FOLDERNAME
mkdir -p $SCRDIR

# Copy the scripts and reference files to the scratch dir
cp \
    sv_pipe.yml \
    ./ref_files/SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz \
    $SCRDIR

cd $SCRDIR

# Load the miniconda3 module will be needed for the doctor_manta.py and for the svafotate run
module load \
    miniconda3/23.11.0 \
    bcftools/1.16

# Check if the environment already exists, if not create it from the yml
ENV_NAME="sv_pipe_conda_env"
if conda info --envs | grep -q "$ENV_NAME"; then
    conda activate "$ENV_NAME"
else
    conda env create -f sv_pipe.yml
    conda activate "$ENV_NAME"
    
    conda install --file https://raw.githubusercontent.com/fakedrtom/SVAFotate/master/requirements.txt
    pip install git+https://github.com/fakedrtom/SVAFotate.git
fi

echo "running svafotate on smoove"
svafotate annotate -v $O_SMOOVE_VCF -b SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz -o $SVAF_SMOOVE_OUTPUT -f 0.8 --cpu 10
bgzip $SVAF_SMOOVE_OUTPUT
echo "smoove svafotate complete"

conda deactivate

#These assume that your proband is at the 0 position if youre called in a different order then you'd need to use sample name to id where your proband is
bcftools view -i '((INFO/SVTYPE="DEL" && FMT/DHFFC[0]<0.7) || (INFO/SVTYPE="DUP" && FMT/DHBFC[0]>1.3) || (INFO/SVTYPE!="DEL" && INFO/SVTYPE!="DUP")) && INFO/Max_AF<0.05 && ((GT[0]="1/1" && GT[1]!="1/1" && GT[2]!="1/1") || ((GT[0]="1/0" || GT[0]="0/1") && GT[1]="0/0" && GT[2]="0/0") || (GT[0]="0/0"))' $SVAF_SMOOVE_OUTPUT -Oz -o $DSMOOVEFILE

# Cleanup
rm run_manta_trio.sh \
    sv_pipe.yml \
    SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz