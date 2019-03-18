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

5. Download the latest genome FASTA file and .gtf annotation for your organism (i.e. [Gencode](http://gencodegenes.org))

6. Build an index for STAR using the following command (executed as a batch script if using a cluster):

```
$ singularity exec -H $PWD -B <path_to_index_directory> shub://jdwheaton/singularity_ngs:star_htseq \
    STAR \
		--runMode genomeGenerate \
		--genomeDir <path_to_index_directory> \
		--genomeFastaFiles <path_to_genome_fasta> \
		--runThreadN <num_threads> \
		--sjdbGTFfile <path_to_genome_annotation> \
		--sjdbOverhang 49
```
More information on STAR can be found [here](https://github.com/alexdobin/STAR).

7. Edit the cluster.json file to suit your particular cluster (i.e. partition name, memory, etc.)

8. Edit the config.yaml file for your experimental setup
    - Path to STAR index
    - Path to transcriptome reference (.gtf file)
    - Name of final .counts file
    - Boolean value (true or false) indicating whether your data are single- or paired-end.
    
9. Edit the `snakemake.sh` submission script to reflect your cluster configuration. Right now, the command works for SLURM systems, but you will need to modify the `-B <path>` to allow singularity access to directories other than your home directory. For example, if you create your STAR index outside of your home directory, you will need to provide either the STAR index directory or a location in its parent tree.

If you get lost, the documenation for Singularity is [here](https://www.sylabs.io/docs/).

10. Run the `snakemake.sh` script as a batch script using your cluster's scheduler. For SLURM, this would be:
```
sbatch snakemake.sh
```

> **NOTE:** If you are using SLURM, you **must** create the `logs/slurm/output` and `logs/slurm/error` directories prior to running snakemake, otherwise your jobs will instantly fail *with no error message*. This can be circumvented by installing a SLURM-specific snakemake [profile](https://github.com/Snakemake-Profiles/slurm) and adjusting the submit script accordingly.


That's it! Snakemake will automagically perform QC, trimming, alignment, and read-counting steps for all of the FASTQ files provided in the `raw_data` folder.

## Workflow description

Your files will be processed as follows:

1. QC reports will be generated for all FASTQ files using [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
2. Adapters are automatically detected and removed from reads using [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
3. Trimmed reads are aligned to the genome using [STAR](https://github.com/alexdobin/STAR).
4. Reads aligning to genes are counted using [featureCounts](http://bioinf.wehi.edu.au/featureCounts/), resulting in a matrix of counts that can be used as input to edgeR or DESeq2 for differential expression analysis.
5. Coverage tracks for visualization (bigWig) are generated using [deepTools](https://deeptools.readthedocs.io/en/develop/) bamCoverage
