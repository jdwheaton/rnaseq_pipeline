# A Snakemake pipeline for RNA-seq analysis

This repository contains everything you need to go from raw FASTQ files to summarized counts ready for differential expression analysis. All of the major software is containerized via Singularity, so all you need is singularity, anaconda and snakemake.

## Dependencies
1. Singularity (requires admin privileges for installation)
2. Anaconda 3
3. Snakemake

## Quick-start guide

1. Install anaconda. If you are installing at the user-level on a remote cluster, I recommend installing [miniconda3](https://docs.conda.io/en/latest/miniconda.html) using `wget`.

2. Use conda to install [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html):

```
$ conda install -c bioconda -c conda-forge snakemake
```

3. Clone this repository:

```
$ git clone https://github.com/jdwheaton/rnaseq_pipeline
```

4. Copy or symlink your raw data (`*.fastq.gz` files) into the `raw_data` folder.

5. Download the latest genome FASTA file and .gtf annotation for your organism (i.e. [Gencode](http://gencodegenes.org)

6. Build an index for STAR using the following command (executed as a batch script if using a cluster):

```
singularity exec -H $PWD shub://jdwheaton/
```

7. Edit the cluster.json file to suit your particular cluster (i.e. partition name, memory, etc.)

8. Edit the config.yaml file for your experimental setup
    - Path to STAR index
    - Path to transcriptome reference (.gtf file)
    - Name of final .counts file
    - Boolean value (True or False) indicating whether your data are single- or paired-end.
