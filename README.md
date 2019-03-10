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
singularity exec -H $PWD -B <path_to_index_directory> shub://jdwheaton/singularity_ngs:star_htseq \
    STAR \
		--runMode genomeGenerate \
		--genomeDir <path_to_index_directory> \
		--genomeFastaFiles <path_to_genome_fasta> \
		--runThreadN <num_threads> \
		--sjdbGTFfile <path_to_genome_annotation> \
		--sjdbOverhang 49
```

7. Edit the cluster.json file to suit your particular cluster (i.e. partition name, memory, etc.)

8. Edit the config.yaml file for your experimental setup
    - Path to STAR index
    - Path to transcriptome reference (.gtf file)
    - Name of final .counts file
    - Boolean value (True or False) indicating whether your data are single- or paired-end.
    
9. Edit the `snakemake.sh` submission script to reflect your cluster configuration. Right now, the command works for SLURM systems, but you will need to modify the `-B <path>` to allow singularity access to directories other than your home directory. For example, if you create your STAR index outside of your home directory, you will need to provide either the STAR index directory or a location in its parent tree.

10. Run the `snakemake.sh` script as a batch script using your cluster's scheduler. For SLURM, this would be:
```
sbatch snakemake.sh
```
That's it! Snakemake will automagically perform QC, trimming, alignment, and read-counting steps for all of the FASTQ files provided in the `raw_data` folder.
