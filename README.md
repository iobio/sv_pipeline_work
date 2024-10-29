# SV QC & Filtering Pipeline

This repo contains a few bioinformatics tools chained together with the goal of creating a conservative structural variant quality and annotation pipeline to complement the mosaic platform and sv.iobio application.

## Outline

- [Background](#background-information)
- [Scripts Information](#scripts)
- [Yaml Examples](#yaml-example)
- [Running Instructions](#running-instructions-in-order)

**NOTE**

As written this is really only intended to be used by our team in our own internal environment. However parts could be modified for other purposes in theory. 

## Background Information

### The `/smoove` folder

The **smoove** folder contains a `.def` file and `.sh` script that can build an apptainer/singularity container for [Smoove](https://github.com/brentp/smoove/tree/master). 

*I had a lot of trouble using the docker container on the HPC environment with apptainer so ultimately had to build it natively with apptainer/singularity.*

I have a `.sif` that can be used in our computing environment already but if you needed to build it again you could from this. **Ask me for the `.sif` if you need and I'll provide the shared path**

---

### Scripts

There are currently two options for template entry scripts `sv_pipe_template.sh` and `sv_pipe_template_smoovestart.sh`. The first is the full pipeline including all tools. The second expects a joint called smoove file and only runs #4 and #5 below on the smoove file.

The **sv_pipe.sh** (whichever you are using) is the **main entry point** which will run a few things. It is intended to be run as a `slurm` script because it is resource intensive and runs for quite a while if running the full shell with the current settings (10+ hours). **The input for the `sv_pipe.sh` is a `.yml` file.**

The full pipeline uses a trio of `CRAM` files in one folder with their associated indexes, an already joint called smoove `vcf` for this same trio, a svafotate population bed file in the `/ref_files`, **AND** a `.sif` file of the smoove tool (See #3 Below). 

- *The smoovestart shell uses a probandId an already joint called smoove `vcf`, and a svafotate population bed file in the `/ref_files`.* Smoovestart should not take nearly as long as the full pipeline.

The cram folder path, sample names, and path to the smoove vcf are defined from within the yaml file (*see example in next section*). 

Example Call:
```
sbatch sv_pipe_personal.sh my_sample_yaml.yml
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
4. It will annotate the vcfs with population AF information using [svafotate](https://github.com/fakedrtom/SVAFotate)
    - The tool looks at multiple SV cohorts and aggrigates them to annotate the variants with various allele frequency metrics.
    - This requires a bed file as defined in the svafotate repo the script expects it to be located in the ref_files directory with the name `SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz`
5. Finally, it will run filtering based on conservative duphold, svafotate criteria, and denovo status (conservative exact match only) and return the filtered down vcfs as `filtered_<smoove/manta>.vcf.gz`

---

### Yaml Example*

```
# FULL PIPE EXAMPLE

# NO TRAILING SLASH ON PATHS

smoove: /folder/folder1/my_samples_smoove.vcf.gz
cramsPath: /folder/folder1/folder_with_crams_and_indexes

# NO FILE EXT ON IDS
probandId: 1234-1
parent1Id: 1234-2
parent2Id: 1234-3

```
```
# SMOOVE START EXAMPLE

# NO TRAILING SLASH ON PATHS

smoove: /folder/folder1/my_samples_smoove.vcf.gz

# NO FILE EXT ON IDS
probandId: 1234-1

```
The `sv_pipe.yml` within the repo is internal to the script and is used to define the python environment for other tools, not as the input to the pipeline.

## Running Instructions In Order

1. Go into your folder within the shared workspace something like: 
```
cd /scratch/ucgd/lustre-labs/marth/scratch/your_uid
```

2. Clone this repository then `cd` into it 
```
git clone https://github.com/iobio/sv_pipeline_work.git

cd sv_pipeline_work
```

3. *smoovestart SKIP --* Get `.sif` location and copy it into the `./smoove` folder 
```
cp /example_path_to_sif.sif ./smoove
```

4. Get the svafotate bed reference and copy it into `./ref_files` with the required name `SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz`
```
wget -O ./ref_files SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz https://zenodo.org/records/11642574/files/SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz
```
**Please note this is the current (18 OCT 2024) location of the recommended bed but this could change take a look at the current recommended bed location on the [svafotate](https://github.com/fakedrtom/SVAFotate?tab=readme-ov-file) repository.**

5. Make your yaml and add your file paths to it using the example in the previous section (do this in the yml_defs folder as this folder is on the .gitignore list)
```
touch ./yml_defs/my_sample_yaml.yml
```
*edit file with information needed*

6. Copy the `sv_pipe_template.sh` into a new file named `sv_pipe_personal.sh` so that you have your own untracked version of the script.
```
cp sv_pipe_template.sh sv_pipe_personal.sh
```

7. Edit all the **< >** fields in the #SBATCH directives in your `sv_pipe_personal.sh`, the list below should be the ones that need editing
```
#SBATCH --job-name=<any_name_you_want>
#SBATCH -o /scratch/ucgd/lustre-labs/marth/scratch/<your_uid>/slurm-%j.out-%N
#SBATCH -e /scratch/ucgd/lustre-labs/marth/scratch/<your_uid>/slurm-%j.err-%N
#SBATCH --mail-user=<your_email>
```
8. Run your personal script
```
sbatch sv_pipe_personal.sh ./yml_defs/my_sample_yaml.yml
```

**If this runs you should see the slurm ouputs appear in your own scratch folder, and a folder should appear that corresponds to your sample analysis, work will take place in that folder aside from the slurm logging/outputs and your final files should be present within that folder as well.**

The final results of this pipeline are two files the `svaf_smoove.vcf.gz` & `svaf_manta.vcf.gz` as well as the corresponding filltered files `filtered_smoove.vcf.gz` and `filtered_manta.vcf.gz`