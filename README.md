# SV QC & Filtering Pipeline

This repo contains a few bioinformatics tools chained together with the goal of creating a conservative structural variant quality and annotation pipeline to complement the mosaic platform and sv.iobio application.

**NOTE**

As written this is really only intended to be used by our team in our own internal environment. However parts could be modified for other purposes in theory. 

## Background Information

### The `/smoove` folder

The **smoove** folder contains a `.def` file and `.sh` script that can build an apptainer/singularity container for [Smoove](https://github.com/brentp/smoove/tree/master). 

*I had a lot of trouble using the docker container on the HPC environment with apptainer so ultimately had to build it natively with apptainer/singularity.*

I have a `.sif` that can be used in our computing environment already but if you needed to build it again you could from this. **Ask me for the `.sif` if you need and I'll provide the shared path**

---

### Scripts

The `sv_pipe.sh` is the **main entry point** which will run a few things. It is intended to be run as a `slurm` script because it is resource intensive and runs for quite a while with the current settings (10+ hours). **The input for the `sv_pipe.sh` is a `.yml` file.**

The pipeline uses a trio of `CRAM` files in one folder with their associated indexes, an already joint called smoove `vcf` for this same trio **AND** a `.sif` file of the smoove tool (See #3 Below). 

The cram folder path, sample names, and path to the smoove vcf are defined from within the yaml file (*see example in next section*)

Example Call:
```
sbatch sv_pipe.sh my_sample_yaml.yml
```

#### The sv_pipe script runs the following additional scripts in this order

1. `run_manta_trio.sh` courtesy of [fakedrtom](https://github.com/fakedrtom)
    - runs a trio with the manta tool on our HPC environment (expects you to be in that environment to access some files and tools).
2. `doctor_manta.py` also provided to me by fakedrtom
    - Puts the neccessary confidence intervals into the joint-called manta file.
3. It will use the `bp_smoove.sif` and run the duphold module on this doctored manta file
    - Annotates the SVs with the depth change information
    - `bp_smoove.sif` is expected to be in the `/smoove` folder in the working directory from which the entry script `sv_pipe.sh` is called
    - Either build the `.sif` yourself or ask me for it. Then copy it into your local repo's `/smoove` folder.
4. Finally, it will annotate the vcfs with population AF information using [svafotate](https://github.com/fakedrtom/SVAFotate)
    - The tool looks at multiple SV cohorts and aggrigates them to annotate the variants with various allele frequency metrics.
    - This requires a bed file as defined in the svafotate repo the script expects it to be located in the ref_files directory with the name `SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz`

---

### Yaml Example*

```

# NO TRAILING SLASH ON PATHS

smoove: /folder/folder1/my_samples_smoove.vcf.gz
cramsPath: /folder/folder1/folder_with_crams_and_indexes

# NO FILE EXT ON IDS
probandId: 1234-1
parent1Id: 1234-2
parent2Id: 1234-3

```

## Running Instructions In Order

