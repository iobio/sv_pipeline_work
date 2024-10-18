#!/bin/bash

#SBATCH --job-name=<any_name_you_want>
#SBATCH --time=14:00:00

# Your Account and Partition as Applicable
#SBATCH --account=marth-rw
#SBATCH --partition=marth-shared-rw

# Output in your directory in the scratch space
#SBATCH -o /scratch/ucgd/lustre-labs/marth/scratch/<your_uid>/slurm-%j.out-%N
#SBATCH -e /scratch/ucgd/lustre-labs/marth/scratch/<your_uid>/slurm-%j.err-%N

# We need at least 4 cpus because smoove will try to use 4 threads by default more is okay but remember each cpu may increase the mem needed
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --mail-type=ALL

# The email to send system messages to
#SBATCH --mail-user=<your_email>

#Inputs to the script
YMLDEF=$1

# Grab the neccessary items from the YMLDEF passed via the command
CRAMSPATH=$(grep 'cramsPath:' "$YMLDEF" | awk '{print $2}')
PROBANDID=$(grep 'probandId:' "$YMLDEF" | awk '{print $2}')
PARENT1ID=$(grep 'parent1Id:' "$YMLDEF" | awk '{print $2}')
PARENT2ID=$(grep 'parent2Id:' "$YMLDEF" | awk '{print $2}')

O_SMOOVE_VCF=$(grep 'smoove:' "$YMLDEF" | awk '{print $2}')

REF_FASTA='human_g1k_v38_decoy_phix.fasta'
FASTA_FOLDER='/scratch/ucgd/lustre/common/data/Reference/GRCh38'

DOCTORED_MANTA_OUTPUT="doctor_manta.vcf.gz"
DUPHOLD_MANTA_OUTPUT="duphold_manta.vcf.gz"
SVAF_SMOOVE_OUTPUT="svaf_smoove.vcf"
SVAF_MANTA_OUTPUT="svaf_manta.vcf"

# Set up scratch directory
FOLDERNAME=${PROBANDID}_${SLURM_JOB_ID}
SCRDIR=/scratch/ucgd/lustre-labs/marth/scratch/$USER/$FOLDERNAME
mkdir -p $SCRDIR

# Copy the scripts and reference files to the scratch dir
cp \
    run_manta_trio.sh \
    doctor_manta.py \
    ./smoove/bp_smoove.sif \
    sv_pipe.yml \
    ./ref_files/SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz \
    $SCRDIR

cd $SCRDIR
mkdir "duphold_run"

# Run the the manta script run_manta_trio.sh (needs the trio's urls CRAM format) -> new joint called vcf
echo "run manta trio"
./run_manta_trio.sh $CRAMSPATH/$PROBANDID.cram $CRAMSPATH/$PARENT1ID.cram $CRAMSPATH/$PARENT2ID.cram
echo "manta complete"

# Load the miniconda3 module will be needed for the doctor_manta.py and for the svafotate run
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
    
    conda install --file https://raw.githubusercontent.com/fakedrtom/SVAFotate/master/requirements.txt
    pip install git+https://github.com/fakedrtom/SVAFotate.git
fi

# Modify the manta vcf (./diploidSV.vcf.gz) with his doctor_manta.py (needs input vcf, and output) -> new doc_vcf
python doctor_manta.py diploidSV.vcf.gz $DOCTORED_MANTA_OUTPUT

BIND_DIR=$(pwd -P)
cp $DOCTORED_MANTA_OUTPUT ./duphold_run/$DOCTORED_MANTA_OUTPUT

# Run smoove duphold on the manta vcf -> dhanno_manta_vcf
echo "Run Duphold"
singularity exec \
    --bind $BIND_DIR/duphold_run:/output \
    --bind $CRAMSPATH:/crams \
    --bind $FASTA_FOLDER:/fastas \
    bp_smoove.sif \
    smoove duphold \
    -f /fastas/$REF_FASTA  \
    -v /output/$DOCTORED_MANTA_OUTPUT \
    -o /output/$DUPHOLD_MANTA_OUTPUT \
    /crams/$PROBANDID.cram /crams/$PARENT1ID.cram /crams/$PARENT2ID.cram
echo "duphold complete"

# Move the duphold output to the main folder for use in svafotate
mv ./duphold_run/$DUPHOLD_MANTA_OUTPUT ./$DUPHOLD_MANTA_OUTPUT

# Run svafotate on both files (.8 ol threshold) -> filtered svaf_vcf
echo "Run svafotate on MANTA"
svafotate annotate -v ./$DUPHOLD_MANTA_OUTPUT -b SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz -o $SVAF_MANTA_OUTPUT -f 0.8 --cpu 4
bgzip $SVAF_MANTA_OUTPUT
echo "manta svafotate complete"

echo "running svafotate on smoove"
svafotate annotate -v $O_SMOOVE_VCF -b SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz -o $SVAF_SMOOVE_OUTPUT -f 0.8 --cpu 4
bgzip $SVAF_SMOOVE_OUTPUT
echo "smoove svafotate complete"

conda deactivate

# Cleanup
rm run_manta_trio.sh \
    runWorkflow.py \
    runWorkflow.py.config.pickle \
    doctor_manta.py \
    bp_smoove.sif \
    sv_pipe.yml \
    SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz

rm -r ./duphold_run