# SV QC & Filtering Pipeline

This repo contains a few bioinformatics tools chained together with the goal of creating a conservative structural variant quality and annotation pipeline to complement the mosaic platform and sv.iobio application.

---
### `smoove` folder

The **smoove** folder contains a `.def` file and `.sh` script that can build an apptainer/singularity container for [Smoove](https://github.com/brentp/smoove/tree/master). 

*I had a lot of trouble using the docker container on the HPC environment with apptainer so ultimately had to build it natively with apptainer/singularity.*

I have a `.sif` that can be used in our computing environment already but if you needed to build it again you could from this.

---

### scripts

The `sv_pipe.sh` is the main entry point which will run a few things. It expects a trio of `CRAM` files, an already joint called smoove `vcf` for this same trio, and the 

1. `run_manta_trio.sh` courtesy of [fakedrtom](https://github.com/fakedrtom)
    - runs a trio with the manta tool on our HPC environment (expects you to be in that environment to access some files and tools).
2. `doctor_manta.py` also provided to me by fakedrtom
    - Puts the neccessary confidence intervals into the joint-called manta file.
3. It will use the `bp_smoove.sif` and run the duphold module on this doctored manta file
    - Annotates the SVs with the depth change information
    - `bp_smoove.sif` is expected to be in the `/smoove` folder in the working directory from which the entry script `sv_pipe.sh` is called
4. Finally, it will annotate the vcfs with population AF information using [svafotate](https://github.com/fakedrtom/SVAFotate)
    - The tool looks at multiple SV cohorts and aggrigates them to annotate the variants with various allele frequency metrics.
    - This requires a bed file as defined in the svafotate repo the script expects it to be located in the ref_files directory with the name `SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz`