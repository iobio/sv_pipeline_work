#!/bin/bash

#SBATCH --job-name=<any_name_you_want>
#SBATCH --time=14:00:00

# Your Account and Partition as Applicable
#SBATCH --account=marth-rw
#SBATCH --partition=marth-shared-rw

# Output in your directory from where the script is launched
#SBATCH -o ./slurm-%j.out-%N
#SBATCH -e ./slurm-%j.err-%N

#SBATCH --ntasks=10
#SBATCH --mem=60G
#SBATCH --mail-type=ALL

# The email to send system messages to
#SBATCH --mail-user=<your_email>

# ===================================== Help message
if [[ "$1" == "--help" ]]; then
    echo "Usage: $0 [options] PROBAND_ID INPUT_VCF VCF_TYPE BUILD [CRAMS_PATH PARENT1_ID PARENT2_ID]"
    echo ""
    echo "Options:"
    echo "  --help                Show this help message"
    echo "  --config FILE         Use a configuration yaml file for input variables"

    echo "Notes:"
    echo ""
    echo " * If using --config, all other options are ignored"
    echo " * If VCF_TYPE is 'smoove' then CRAMS_PATH, PARENT1_ID, and PARENT2_ID are not required"
    echo " * BUILD must be either 'GRCh37' or 'GRCh38'"
    echo " * VCF_TYPE must be either 'smoove', 'dragen', or 'manta'"
    exit 0
fi

# ===================================== Configuration file handling
if [[ "$1" == "--config" ]]; then
    CONFIG_FILE=$2

    if [[ ! -f "$CONFIG_FILE" ]]; then
        echo "Error: Configuration file not found!"
        exit 1
    fi

    # Set default value to null if not found
    PROBAND_ID=${PROBAND_ID:-""}
    INPUT_VCF=${INPUT_VCF:-""}
    VCF_TYPE=${VCF_TYPE:-""}
    BUILD=${BUILD:-""}
    CRAMS_PATH=${CRAMS_PATH:-""}
    PARENT1_ID=${PARENT1_ID:-""}
    PARENT2_ID=${PARENT2_ID:-""}

    PROBAND_ID=$(grep 'PROBAND_ID:' "$CONFIG_FILE" | awk '{print $2}')
    INPUT_VCF=$(grep 'INPUT_VCF:' "$CONFIG_FILE" | awk '{print $2}')
    VCF_TYPE=$(grep 'VCF_TYPE:' "$CONFIG_FILE" | awk '{print $2}')
    BUILD=$(grep 'BUILD:' "$CONFIG_FILE" | awk '{print $2}')
    CRAMS_PATH=$(grep 'CRAMS_PATH:' "$CONFIG_FILE" | awk '{print $2}')
    PARENT1_ID=$(grep 'PARENT1_ID:' "$CONFIG_FILE" | awk '{print $2}')
    PARENT2_ID=$(grep 'PARENT2_ID:' "$CONFIG_FILE" | awk '{print $2}')
else
    # The ID we should look for in the VCF
    PROBAND_ID=$1
    # The path to the input VCF file
    INPUT_VCF=$2
    # Allowed: smoove, dragen, manta
    VCF_TYPE=$3
    # Allowed: GRCh37 or GRCh38
    BUILD=$4
    if [[ "$VCF_TYPE" != "smoove" ]]; then
        CRAMS_PATH=$5
        PARENT1_ID=$6
        PARENT2_ID=$7
    else 
        CRAMS_PATH=""
        PARENT1_ID=""
        PARENT2_ID=""
    fi
fi

# Check if all required variables are set
if [[ -z "$PROBAND_ID" || -z "$INPUT_VCF" || -z "$VCF_TYPE" || -z "$BUILD" ]]; then
    echo "Error: Missing required input variables."
    echo "Use --help for usage information."
    exit 1
fi

if [[ "$VCF_TYPE" != "smoove" && ( -z "$CRAMS_PATH" || -z "$PARENT1_ID" || -z "$PARENT2_ID" ) ]]; then
    echo "Error: Missing required CRAM files for non-smoove VCF type."
    echo "Use --help for usage information."
    exit 1
fi

# Assign the ref fasta and fasta folder depending on the build
# This will only be used ultimately if we run duphold (i.e. we don't have a smoove file)
if [[ $BUILD == "GRCh37" ]]; then
    REF_FASTA="human_g1k_v37_decoy_phix.fasta"
    FASTA_FOLDER="/scratch/ucgd/lustre-core/common/data/Reference/GRCh37"
elif [[ $BUILD == "GRCh38" ]]; then
    REF_FASTA="human_g1k_v38_decoy_phix.fasta"
    FASTA_FOLDER="/scratch/ucgd/lustre/common/data/Reference/GRCh38"
else
    echo "Error: Build is not valid (GRCh37 or GRCh38 accepted)"
    exit 1
fi

# Variables for the output file names
# Duphold (dh) files only used if we have something other than a smoove vcf
DH_DOCTORED_OUTPUT="dh_doctored.vcf.gz"
DUPHOLD_OUTPUT="duphold_annotated.vcf.gz"
SVAF_OUTPUT="svaf_annot.vcf"
FILTERED_VCF="calypso_filtered_svs.vcf.gz"

# ===================================== Begin Setup
FOLDERNAME=${PROBAND_ID}_calypso_sv
SCRDIR=./$FOLDERNAME
mkdir $SCRDIR
cd $SCRDIR

# Copy the scripts and reference files to the scratch dir
cp \
    /scratch/ucgd/lustre-labs/marth/scratch/u1069837/data/sv_pipeline_work/doctor_manta.py \
    /scratch/ucgd/lustre-labs/marth/scratch/u1069837/data/sv_pipeline_work/smoove/bp_smoove.sif \
    /scratch/ucgd/lustre-labs/marth/scratch/u1069837/data/sv_pipeline_work/ref_files/SVAFotate_core_SV_popAFs.GRCh37.v4.1.bed.gz \
    /scratch/ucgd/lustre-labs/marth/scratch/u1069837/data/sv_pipeline_work/ref_files/SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz \
    .

# Creating this here because even if it is empty and ultimately isn't used it is part of the cleanup so I want to make sure it is created anyway
mkdir "duphold_run"

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
    if [[ ! -f "sv_pipe.yml" ]]; then
        echo "Error: sv_pipe.yml file not found!"
        exit 1
    fi
    conda env create -f sv_pipe.yml
    conda activate "$ENV_NAME"
    
    conda install --file https://raw.githubusercontent.com/fakedrtom/SVAFotate/master/requirements.txt
    pip install git+https://github.com/fakedrtom/SVAFotate.git
fi

# If we have something other than a smoove vcf then we need to doctor the VCF and then run duphold before svafotate
if [[ $VCF_TYPE != "smoove" ]]; then
    python doctor_manta.py $INPUT_VCF $DH_DOCTORED_OUTPUT

    BIND_DIR=$(pwd -P)
    cp $DH_DOCTORED_OUTPUT ./duphold_run/$DH_DOCTORED_OUTPUT

    # IMPORTANT: This is the main reason we need to think about the input structure of the cram files
    echo "Run Duphold"
    singularity exec \
        --bind $BIND_DIR/duphold_run:/output \
        --bind $CRAMSPATH:/crams \
        --bind $FASTA_FOLDER:/fastas \
        bp_smoove.sif \
        smoove duphold \
        -f /fastas/$REF_FASTA  \
        -v /output/$DH_DOCTORED_OUTPUT \
        -o /output/$DUPHOLD_OUTPUT \
        /crams/$PROBAND_ID.cram /crams/$PARENT1_ID.cram /crams/$PARENT2_ID.cram
    echo "duphold complete"

    # Move the duphold output to the main folder for use in svafotate
    mv ./duphold_run/$DUPHOLD_OUTPUT ./$DUPHOLD_OUTPUT

    # Run svafotate
    echo "running svafotate"
    svafotate annotate -v $DUPHOLD_OUTPUT -b SVAFotate_core_SV_popAFs.GRCh37.v4.1.bed.gz -o $SVAF_OUTPUT -f 0.8 --cpu 10
    bgzip $SVAF_OUTPUT
    echo "svafotate complete"

    conda deactivate

# If it is smoove duphold has already been run just run svafotate
else
    echo "running svafotate"
    svafotate annotate -v $INPUT_VCF -b SVAFotate_core_SV_popAFs.GRCh37.v4.1.bed.gz -o $SVAF_OUTPUT -f 0.8 --cpu 10
    bgzip $SVAF_OUTPUT
    echo "svafotate complete"

    conda deactivate
fi

# Load sample names from the VCF file into the SAMPLES array
SAMPLES=($(bcftools query -l "$SVAF_OUTPUT.gz"))
PROB_INDEX=-1

# Find the index of the proband sample in the SAMPLES array
for i in "${!SAMPLES[@]}"; do
    if [[ "${SAMPLES[$i]}" == "$PROBAND_ID" ]]; then
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
bcftools view -i "((INFO/SVTYPE=\"DEL\" && FMT/DHFFC[$PROB_INDEX]<0.7) || \
                   (INFO/SVTYPE=\"DUP\" && FMT/DHBFC[$PROB_INDEX]>1.3) || \
                   (INFO/SVTYPE!=\"DEL\" && INFO/SVTYPE!=\"DUP\")) && \
                   INFO/Max_AF<0.05" "$SVAF_OUTPUT.gz" -Oz -o "$FILTERED_VCF"

# Cleanup
rm bp_smoove.sif \
    SVAFotate_core_SV_popAFs.GRCh37.v4.1.bed.gz \
    SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz \
    doctor_manta.py \
    dh_doctored.vcf.gz \
    svaf_annot.vcf.gz \
    duphold_annotated.vcf.gz 

rm -r ./duphold_run

# FINAL FILE: calypso_filtered_svs.vcf.gz

# This final file will have svafotate annotations, duphold annotations, and will be filtered based on the following
#   svafotate max af < .05
#   deletions - duphold dhffc < .7
#   dupolications - duphold dhbfc > 1.3
#   insertions or other types will be allowed any duphold value but filtered on svafotate af