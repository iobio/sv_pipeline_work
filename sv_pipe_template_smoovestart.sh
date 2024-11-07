#!/bin/bash

# THIS SCRIPT IS MINIMIZED FOR SMOOVE RUN ONLY
#SBATCH --job-name=<any_name_you_want>
#SBATCH --time=10:00:00

# Your Account and Partition as Applicable
#SBATCH --account=marth-rw
#SBATCH --partition=marth-shared-rw

# Output in your directory in the scratch space
#SBATCH -o ./slurm-%j.out-%N
#SBATCH -e ./slurm-%j.err-%N

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
DSMOOVEFILE=filtered_smoove.vcf.gz

# Set up scratch directory
FOLDERNAME=${PROBANDID}_${SLURM_JOB_ID}
SCRDIR=./$FOLDERNAME
mkdir $SCRDIR
cd $SCRDIR

# Copy the scripts and reference files to the scratch dir
cp \
    /scratch/ucgd/lustre-labs/marth/scratch/u1069837/data/sv_pipeline_work/sv_pipe.yml \
    /scratch/ucgd/lustre-labs/marth/scratch/u1069837/data/sv_pipeline_work/ref_files/SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz \
    .

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

# Assign the other two indices dynamically
OTHER_INDICES=($(seq 0 $((${#SAMPLES[@]} - 1)) | grep -v "$PROB_INDEX"))

# Ensure we only have two other indices to work with
if [[ ${#OTHER_INDICES[@]} -ne 2 ]]; then
    echo "Error: Expected exactly two other samples, found ${#OTHER_INDICES[@]}."
    exit 1
fi

OTHER_INDEX1=${OTHER_INDICES[0]}
OTHER_INDEX2=${OTHER_INDICES[1]}

# Construct the bcftools command using the dynamically determined indices
bcftools view -i '((INFO/SVTYPE="DEL" && FMT/DHFFC['"$PROB_INDEX"']<0.7) || \
                   (INFO/SVTYPE="DUP" && FMT/DHBFC['"$PROB_INDEX"']>1.3) || \
                   (INFO/SVTYPE!="DEL" && INFO/SVTYPE!="DUP")) && \
                   INFO/Max_AF<0.05 && \
                   ((GT['"$PROB_INDEX"']="1/1" && GT['"$OTHER_INDEX1"']!="1/1" && GT['"$OTHER_INDEX2"']!="1/1") || \
                   ((GT['"$PROB_INDEX"']="1/0" || GT['"$PROB_INDEX"']="0/1") && GT['"$OTHER_INDEX1"']="0/0" && GT['"$OTHER_INDEX2"']="0/0") || \
                   (GT['"$PROB_INDEX"']="0/0"))' "$SVAF_SMOOVE_OUTPUT.gz" -Oz -o "$DSMOOVEFILE"

# Construct the bcftools command using the dynamically determined indices
bcftools view -i '((INFO/SVTYPE="DEL" && FMT/DHFFC['"$PROB_INDEX"']<0.7) || \
                   (INFO/SVTYPE="DUP" && FMT/DHBFC['"$PROB_INDEX"']>1.3) || \
                   (INFO/SVTYPE!="DEL" && INFO/SVTYPE!="DUP")) && \
                   INFO/Max_AF<0.05 && \
                   ((GT['"$PROB_INDEX"']="1/1" && GT['"$OTHER_INDEX1"']!="1/1" && GT['"$OTHER_INDEX2"']!="1/1") || \
                   ((GT['"$PROB_INDEX"']="1/0" || GT['"$PROB_INDEX"']="0/1") && GT['"$OTHER_INDEX1"']="0/0" && GT['"$OTHER_INDEX2"']="0/0") || \
                   (GT['"$PROB_INDEX"']="0/0"))' "$SVAF_MANTA_OUTPUT.gz" -Oz -o "$DMANTAFILE"
# Cleanup
rm run_manta_trio.sh \
    sv_pipe.yml \
    SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz