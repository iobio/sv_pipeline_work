#!/bin/bash
#Tom's *fakedrtom* Script for running manta jointly on a trio

#This will obviously only work on our environment as it requires some files located in our computing server system

#Set some variables
configManta='/scratch/ucgd/lustre-work/quinlan/u0055382/src/Manta/manta-1.6.0.centos6_x86_64/bin/configManta.py'
MEM=384
CPU=48
FASTA='/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta'
BED='/scratch/ucgd/lustre-work/quinlan/u0055382/src/Manta/canonical_chroms.bed.gz'

#Feed it CRAM paths from command line
PBDCRAM=$1
MOMCRAM=$2
DADCRAM=$3

$configManta \
--bam $PBDCRAM \
--bam $MOMCRAM \
--bam $DADCRAM \
--referenceFasta $FASTA \
--generateEvidenceBam \
--callRegions $BED \
--runDir ./

./runWorkflow.py -j $CPU -g $MEM

ln -s results/variants/diploidSV.vcf.gz ./
