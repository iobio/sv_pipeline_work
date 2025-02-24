#!/bin/bash

#SBATCH --job-name=sv_pipeline
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

REF_FASTA='human_g1k_v37_decoy_phix.fasta'
FASTA_FOLDER='/scratch/ucgd/lustre-core/common/data/Reference/GRCh37'

DOCTORED_MANTA_OUTPUT="doctor_manta.vcf.gz"
DUPHOLD_MANTA_OUTPUT="duphold_manta.vcf.gz"
SVAF_SMOOVE_OUTPUT="svaf_smoove.vcf"
SVAF_MANTA_OUTPUT="svaf_manta.vcf"
DMANTAFILE="filtered_manta.vcf.gz"
DSMOOVEFILE="filtered_smoove.vcf.gz"

# Set up scratch directory
FOLDERNAME=${PROBANDID}_${SLURM_JOB_ID}
SCRDIR=./$FOLDERNAME
mkdir $SCRDIR
cd $SCRDIR

# Copy the scripts and reference files to the scratch dir
cp \
    /scratch/ucgd/lustre-labs/marth/scratch/u1069837/data/sv_pipeline_work/run_manta_trio_37.sh \
    /scratch/ucgd/lustre-labs/marth/scratch/u1069837/data/sv_pipeline_work/doctor_manta.py \
    /scratch/ucgd/lustre-labs/marth/scratch/u1069837/data/sv_pipeline_work/smoove/bp_smoove.sif \
    /scratch/ucgd/lustre-labs/marth/scratch/u1069837/data/sv_pipeline_work/sv_pipe.yml \
    /scratch/ucgd/lustre-labs/marth/scratch/u1069837/data/sv_pipeline_work/ref_files/SVAFotate_core_SV_popAFs.GRCh37.v4.1.bed.gz \
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
#TODO: Change both references to the 37 version
echo "Run svafotate on MANTA"
svafotate annotate -v ./$DUPHOLD_MANTA_OUTPUT -b SVAFotate_core_SV_popAFs.GRCh37.v4.1.bed.gz -o $SVAF_MANTA_OUTPUT -f 0.8 --cpu 12
bgzip $SVAF_MANTA_OUTPUT
echo "manta svafotate complete"

echo "running svafotate on smoove"
svafotate annotate -v $O_SMOOVE_VCF -b SVAFotate_core_SV_popAFs.GRCh37.v4.1.bed.gz -o $SVAF_SMOOVE_OUTPUT -f 0.8 --cpu 12
bgzip $SVAF_SMOOVE_OUTPUT
echo "smoove svafotate complete"

conda deactivate

# Load sample names from the VCF file into the SAMPLES array
SAMPLES=($(bcftools query -l "$SVAF_SMOOVE_OUTPUT.gz"))
PROB_INDEX=-1

# Find the index of the proband sample in the SAMPLES array
for i in "${!SAMPLES[@]}"; do
    if [[ "${SAMPLES[$i]}" == "$PROBANDID" ]]; then
        PROB_INDEX=$i
        break
    fi
done

# Check if the proband index was found
if [[ "$PROB_INDEX" -lt 0 ]]; then
    echo "Error: Proband sample not found in VCF file."
    exit 1
fi

# Construct the bcftools command using the dynamically determined indices
bcftools view -i '((INFO/SVTYPE="DEL" && FMT/DHFFC['"$PROB_INDEX"']<0.7) || \
                   (INFO/SVTYPE="DUP" && FMT/DHBFC['"$PROB_INDEX"']>1.3) || \
                   (INFO/SVTYPE!="DEL" && INFO/SVTYPE!="DUP")) && \
                   INFO/Max_AF<0.05' "$SVAF_SMOOVE_OUTPUT.gz" -Oz -o "$DSMOOVEFILE"

# Construct the bcftools command using the dynamically determined indices
bcftools view -i '((INFO/SVTYPE="DEL" && FMT/DHFFC['"$PROB_INDEX"']<0.7) || \
                   (INFO/SVTYPE="DUP" && FMT/DHBFC['"$PROB_INDEX"']>1.3) || \
                   (INFO/SVTYPE!="DEL" && INFO/SVTYPE!="DUP")) && \
                   INFO/Max_AF<0.05' "$SVAF_MANTA_OUTPUT.gz" -Oz -o "$DMANTAFILE"

# Cleanup
rm run_manta_trio.sh \
    runWorkflow.py \
    runWorkflow.py.config.pickle \
    diploidSV.vcf.gz \
    doctor_manta.py \
    doctor_manta.vcf.gz \
    duphold_manta.vcf.gz \
    bp_smoove.sif \
    sv_pipe.yml \
    SVAFotate_core_SV_popAFs.GRCh37.v4.1.bed.gz \
    workflow.error.log.txt \
    workflow.exitcode.txt \
    workflow.warning.log.txt


rm -r ./duphold_run
rm -r ./results
rm -r ./workspace