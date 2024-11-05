#!/bin/bash

#SBATCH --job-name=<sv_pipeline>
#SBATCH --time=14:00:00

#SBATCH --account=marth-rw
#SBATCH --partition=marth-rw

#SBATCH -o ./slurm-%j.out-%N
#SBATCH -e ./slurm-%j.err-%N

#SBATCH --nodes=1
#SBATCH --mail-type=ALL

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
DMANTAFILE=filtered_manta.vcf.gz
DSMOOVEFILE=filtered_smoove.vcf.gz

# Set up scratch directory
FOLDERNAME=${PROBANDID}_${SLURM_JOB_ID}
SCRDIR=./$FOLDERNAME
mkdir $SCRDIR
cd $SCRDIR

# Copy the scripts and reference files to the scratch dir
cp \
    run_manta_trio.sh \
    doctor_manta.py \
    ./smoove/bp_smoove.sif \
    sv_pipe.yml \
    ./ref_files/SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz \
    .

mkdir "duphold_run"

# Run the the manta script run_manta_trio.sh (needs the trio's urls CRAM format) -> new joint called vcf
echo "run manta trio"
./run_manta_trio.sh $CRAMSPATH/$PROBANDID.cram $CRAMSPATH/$PARENT1ID.cram $CRAMSPATH/$PARENT2ID.cram
echo "manta complete"

# Load the miniconda3 module will be needed for the doctor_manta.py and for the svafotate run
module load \
    miniconda3/23.11.0 \
    singularity/4.1.1 \
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
svafotate annotate -v ./$DUPHOLD_MANTA_OUTPUT -b SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz -o $SVAF_MANTA_OUTPUT -f 0.8 --cpu 12
bgzip $SVAF_MANTA_OUTPUT
echo "manta svafotate complete"

echo "running svafotate on smoove"
svafotate annotate -v $O_SMOOVE_VCF -b SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz -o $SVAF_SMOOVE_OUTPUT -f 0.8 --cpu 12
bgzip $SVAF_SMOOVE_OUTPUT
echo "smoove svafotate complete"

conda deactivate

bcftools view -i '((INFO/SVTYPE="DEL" && FMT/DHFFC[0]<0.7) || (INFO/SVTYPE="DUP" && FMT/DHBFC[0]>1.3) || (INFO/SVTYPE!="DEL" && INFO/SVTYPE!="DUP")) && INFO/Max_AF<0.05 && ((GT[0]="1/1" && GT[1]!="1/1" && GT[2]!="1/1") || ((GT[0]="1/0" || GT[0]="0/1") && GT[1]="0/0" && GT[2]="0/0") || (GT[0]="0/0"))' $SVAF_MANTA_OUTPUT -Oz -o $DMANTAFILE
bcftools view -i '((INFO/SVTYPE="DEL" && FMT/DHFFC[0]<0.7) || (INFO/SVTYPE="DUP" && FMT/DHBFC[0]>1.3) || (INFO/SVTYPE!="DEL" && INFO/SVTYPE!="DUP")) && INFO/Max_AF<0.05 && ((GT[0]="1/1" && GT[1]!="1/1" && GT[2]!="1/1") || ((GT[0]="1/0" || GT[0]="0/1") && GT[1]="0/0" && GT[2]="0/0") || (GT[0]="0/0"))' $SVAF_SMOOVE_OUTPUT -Oz -o $DSMOOVEFILE

# Cleanup
rm run_manta_trio.sh \
    runWorkflow.py \
    runWorkflow.py.config.pickle \
    doctor_manta.py \
    bp_smoove.sif \
    sv_pipe.yml \
    SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz

rm -r ./duphold_run