# MC-spike: Assessment of bias and reproducibility of viral metagenomics methods

The pipeline consists of a Nextflow script which runs the bioinformatics tools that perform the bulk of the analysis. A collection of R and R Markdown scripts is provided to perform the final part of the analysis, combining output files from the Nextflow pipeline to generating figures and statistics.

## Requirements

The Nextflow pipeline has been tested on the following software:

- Nextflow 24.04.3
- Ubuntu 22.04.4 running Conda 24.3.0
- AlmaLinux 9.4 running Singularity 3.8.7-1.el9.centos
- CheckV 1.0.3 scripts

The R scripts have been tested with R 4.1.1 on:

- Windows 10 and 11
- Ubuntu 22.04.4

## Installation

1. Install Nextflow by following the instructions on https://www.nextflow.io/
2. Install R version 4.4.1 from https://www.r-project.org/
3. Download the MC Spike pipeline from the git repository using `git clone https://github.com/RHaagmans/mc-spike/mc-spike.git`
4. In the MC Spike project folder, run `./setup_renv.sh` to set up the project R environment. 
5. Download the reads to the `./reads/raw` folder of the mc-spike project.
6. Download the CheckV `anicalc.py` and `aniclyst.py` scripts from `https://bitbucket.org/berkeleylab/checkv/src/master/scripts/`, and place them in the `./scripts/` folder. 
7. Optional: obtain a Singularity or Docker image (see below)

### Singularity and Docker images

When running Nextflow with Singularity or Docker, the relevant image can be pulled automatically from Singularity Library or Docker Hub. Alternatively, you can obtain a Singularity image manually and  move into the MC Spike `./envs` directory. 

For Singularity, either:
    1. Pull the image using `singularity pull --arch amd64 library://rhaagmans/mc-spike/mc-spike:v1.0`
    2. Build the image using `sudo singularity build mc-spike.img ubuntu-mc.def`

For Docker, either:
    1. Pull the image using `docker pull rhaagmans/mc-spike:1.0`
    2. Build the image using `docker build . rhaagmans/mc-spike:1.0`
    
### Conda

The pipeline can also be run in a Conda environment without the use of Docker or Singularity. If the server running the pipeline is connected to the internet, Nextflow can handle setting up the environment through the provided `./envs/ubuntu-mc-environment.yml` environment file.  It is also possible to set up the environment yourself by moving into the `./envs` folder and running `conda env create -n mc -f ubuntu-mc-environment.yml`. If this fails for some reason, you could try `ubuntu-mc-environment-history.yml` instead, which only lists the required packages without dependencies. The `conda` parameter in the `conda` profile of `main.config` needs to be edited to point to the location of the environment, e.g. `conda = "/opt/miniforge3/envs/mc"`.

### Databases

The pipeline uses CheckV and geNomad, which each depend on their databases. The first step of the pipeline consists of downloading these databases, and if your computing setup does not allow access to online resources, these can be downloaded manually to a folder of your choice. If the files are downloaded manually, their location can be specified in the configuration file (see below). Or, the script can be run using the `--download_only` option on a computer with internet access to download these online resources: `nextflow run ${mc_spike}/main.nf -c ${mc_spike}/config/main.config -profile vm,conda --download_only` (select the appropriate profile for your situation). By default, the pipeline downloads the latest database versions. If you already have downloaded these databases, you can place them in a folder and add this to the configuration (see below). For the paper, geNomad database version 1.7 and CheckV database version 1.5 were used. 

## Configuring the pipeline

There are two main files that can be edited to configure the pipeline, including for analysis of your own data. The first is `./config/main.config`, in which the nextflow script is configured. This includes changing paths to various files and folders, and setting the computing resources for the various steps of the pipeline. These parameters are configured separately for the `vm` and `hpc` profiles. The other file is `./config/sample_sheet.csv`, which defines the sample metadata and raw read file names.

### main.config

The `samples` parameter points to the sample metadata sheet `sample_sheet.csv` and can be changed to a user-supplied metadata sheet.

The `reads` parameter defines the location of the raw reads. By default, the raw read files are assumed to be located in the `./reads/raw` folder in the main project folder.

The `data_out` parameter defines the location of the pipeline output. By default, the pipeline output will be stored in the `mc-spike` project folder that also contains the code. However, a different folder can be specified by changing the `data_out` parameter in the `vm` and/or `hpc` profile. All output will be stored within the folder specified there:

- Reads are stored in `${data_out}/reads`.
- Output from bioinformatics tools are stored in `${data_out}/data`.  
- Figures from the R scripts stored in `${data_out}/figures`.
- R Markdown reports stored in `${data_out}/reports`. 

The `db_dir` parameter defines the download location of CheckV and geNomad databases. By default, the pipeline will download CheckV and geNomad databases to `./dbs` in the `mc-spike` project folder. If you've already downloaded these databases (for the manuscript, CheckV DB version 1.5 and geNomad DB version 1.7 were used), you can point this parameter to a folder containing (symbolic links to) these databases. The CheckV database files should be in `${db_dir}/checkv_db/checkv-db-v*`. Make sure there is only one database folder matching that pattern in the `${db_dir}/checkv_db/` folder. The geNomad database files should be in `${db_dir}/genomad_db`. 

The `ictv_vmr` parameter defines the path to the ICTV virus metadata resource, which is downloaded to `${data_out}/data/taxonomy` by default. This can be supplied by the user separately.

The `assembly_groups` parameter defines the columns in the sample sheet that define the assembly groups. By default, individual samples are assembled, and all samples belonging to the same base faecal sample, processed using the same method, are co-assembled. These default assembly groups are defined by the `assembly` and `coassembly` columns in `sample_sheet.csv`. The sample sheet can be edited to define alternative assembly groupings, and are configured by specifying the column names in the `assembly_groups` parameter. Only use alphanumeric characters, underscores (`_`), and dashes (`-`) in the to define assembly groups.

The `adapters` parameter can be defined by the user to supply custom adapter sequences to fastp in the read trimming and filtering step.

The `megahit_tmp` parameter can be used to specify an alternative location for the storage of temporary files for MEGAHIT.

Additional pipeline configurations, including CPU and memory requirements for processes, nextflow executors (e.g, SLURM), total available CPUs and memory, etc., can be changes in `main.config` as well. For more information, see the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html). 

### sample_sheet.csv

The sample sheet consists of the following columns:

- `sample_name`: unique identifier for each sample (including identical replicates).
- `base_sample`: identifier for the base sample to which a mock community was added.
- `spike`: identifier for the mock community added to the sample
- `rt_method`: identifier for the method used to process the sample
- `replicate`: identifier for the replicate, in case multiple replicates of base_sample + spike + rt_method are used
- `assembly`: specifies single assemblies. 
- `coassembly`: specifies groups of samples that are co-assembled.
- `fw_readfile` and `rv_readfile`: the base filenames of the forward and reverse read files in the reads folder specified in `main.config` (excluding the path).

For the identifiers, please only use alphanumeric characters, underscores (`_`), and dashes (`-`). 

Samples can be co-assembled based on user-defined criteria, either by changing the groups in the `coassembly` column, or by adding (a) additional column(s) that specifiy additional assembly groups. For example, if you would like to co-assemble both all sample that share a base sample and method, and all samples that share a base sample (thus co-assemlbing samples from multiple methods), an additional column can be added, e.g. `coasembly_by_base`. Samples can be excluded from an assembly grouping by either leaving the cell blank or entering a single dash (`-`).

If the assembly groups are changed, or multiple replicates of identical treatments (base sample + spike + method) are used, the R scripts will need to be modified.

### Mock Community particle counts and reference names

When analysing your own data, there are two additional files that the R reports depend on. The file `mock_community_particles_added.csv` contains three columns with info on the number of particles added to the Mock Community for each virus, and the number of strands (single our double). The file `mock_community_ref_to_virus.csv` contains a table translating the sequence ID in the mock community fasta file to an individual virus. The `--vlp_counts` and `--mc_refnames` options can be used to change the files.

## Running the pipeline

The nextflow pipeline is run in a temporary folder, for example on a scratch disk, which requires ~150 Gb of disk space. Most ouput files are copied to the `${data_out}` folder, although some larger files, particularly `*.sam` files, are instead symlinked to save space. The process parameter `publishDir` option `mode` can be changed from `"symlink"` to `"copy"` for those files to include the full file in the output folder. The pipeline can be rerun using the Nextflow `-resume` option from the same location to resume after an error, or to copy files to the `${data_out}`. The Nextflow run folder is not required for subsequent analysis and can be deleted. 

Two main profiles are specified in `main.config`, `vm` and `hpc`. By default, the `vm` profile is set up to run on a single computer. The `hpc` profile is set up to run on a high performance computing cluster running the SLURM workload manager.

Additionally, container profiles are provied for Docker (`docker`), Singularity (`singularity_local` and `singularity_online`), and Conda (`conda`). For Docker and Singularity, Nextflow can handle downloading the images (in the case of Singularity, use profile `singularity_online`). Alternatively, you can pull or build the images as described above. In that case, use `singularity_local` profile and make sure the configuration file points to the singularity image using the `container` parameter. The `docker` profile works in both cases. For Conda, an environment file is provided in the `conda` profile and Nextflow will handle setting up the environment. The `conda` parameter can be edited to point to a different local conda environment. The profiles can be selected using the Nextflow `-profile` option and should always consist of `vm` or `hpc` and one of the container profiles, separated by a comma (e.g., `-profile vm,docker` or `-profile hpc,singularity_local`).

To run the pipeline:

1. Make sure the configuration file points to the correct file paths. 
2. Make sure the sample sheet contains the correct metadata and file names.
3. Optionally, pre-download the Singularity or Docker image or create a Conda environment.
4. On your server, move into a temporary folder with sufficient disk space.
5. Run the nextflow pipeline using `nextflow run ${mc_spike}/main.nf -config ${mcspike}/config/main.config -profile hpc,singularity_local` (or another `-profile`, as described above). Here, `$mcspike` points to the MC Spike project folder and the profile can be selected using the `-profile` option. A different configuration file can be used with the Netflow `-config` option.
6. Once Nextflow has finished, run `./initiate_r_analysis.sh`, or run `Rscript ./R/process_data.R`, followed by `Rscript ./R/generate_reports.R` to generate figures and reports.

## Citation

A manuscript is in preparation, meanwhile, a preprint is available at:

> Haagmans, R.; Charity, O. J.; Baker, D.; Telatin, A.; Savva, G. M.; Adriaenssens, E. M.; Powell, P.; Carding, S. R. Assessing Bias and | Reproducibility of Viral Metagenomics Methods for the Combined Detection of Faecal RNA and DNA Viruses. Preprints 2024, 2024112016. https://doi.org/10.20944/preprints202411.2016.v1
