# SV Pipeline Notes
**Formalized Documentation of the SV QC Pipeline**
- This pipeline is being created in order to pre-filter and annotate structural variants in a resonable
way for use within mosaic and downstream with the SV.iobio application.
- As a first pass the pipeline will be conservative, filtering out only those variants which we are most sure are not relevant or trustworthy.

## Goals
- Document the tools used in the pipeline
- Document and formalize the steps in the pipeline

## Callers
By default the mosaic pipeline will expect VCFs from the smoove pipeline as well as manta called 
VCFs.

### Smoove
> Smoove Git: [link-to-smoove-git](https://github.com/brentp/smoove)

### Manta
> Manta Git: [link-to-manta-git](https://github.com/Illumina/manta/tree/master)

## Tools

### Duphold
This tool is used for annotating SV vcfs with scores related to depth information. These can be used
to support or create skepticism around SV calls with discordant signals/call types.

> Duphold Git: [link-to-duphold-git](https://github.com/brentp/duphold/tree/master)

### Svafotate
This tool is used for reviewing the 'commonness' of particular variants. Variants that are commonly 
found in population data bases can be thought of as being of less relevance in a rare disease context.

> Svafotate Git: [link-to-svafotate-git](https://github.com/fakedrtom/SVAFotate)

## Procedures
- Get Calls
    - Manta 
        - UCGD calls manta seperately but Tom calls jointly with trio or quad whatever he has.
    - Smoove
        - is already going to include duphold
    - Svafotate
        - requires a few things (namely a big .bed file to go over)

### Inputs
- BAM, BAI, Fasta (for each samples)
- Bed (excluded regions) - Maybe Package Ours by Default? `investigate`
- Bed (svafotate) - For the population overlaps but maybe we can package this? It may make the docker image too large? `investigate`
