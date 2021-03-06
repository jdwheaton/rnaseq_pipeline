configfile: "config.yaml"

SAMPLES, = glob_wildcards("raw_data/{smp}.fastq.gz")
INDEX_DIR = config["star_index"]
ANNOTATION = config["annotation"]
STRANDED = config["stranded"]

COUNT_FILENAME = [config["count_filename"]]

ALL_FASTQC = expand("fastqc_out/{smp}_fastqc.zip", smp=SAMPLES)
ALL_BAMCOV = expand("bigwigs/{sample}.rpkm.bw", sample=SAMPLES)

rule all:
    input: ALL_BAMCOV + ALL_FASTQC + COUNT_FILENAME + ["multiqc_report.html"]

rule fastqc:
    input:
        "raw_data/{sample}.fastq.gz"
    output:
        "fastqc_out/{sample}_fastqc.html",
        "fastqc_out/{sample}_fastqc.zip"
    threads: 1
    singularity:
        "shub://jdwheaton/singularity-ngs:qc_trim"
    log:
        "logs/{sample}.fastqc.out"
    shell:
        "fastqc {input} -o fastqc_out/ &> {log}"

rule trim_galore:
    input:
        "raw_data/{sample}.fastq.gz"
    output:
        "trimmed_fastq/{sample}_trimmed.fq"
    threads: 1
    singularity:
        "shub://jdwheaton/singularity-ngs:qc_trim"
    log:
        "logs/{sample}.trimgalore.out"
    shell:
        "trim_galore --dont_gzip --fastqc -o trimmed_fastq/ {input} &> {log}"

rule star_align:
    input:
        "trimmed_fastq/{sample}_trimmed.fq"
    output:
        temp("alignment/{sample}.Aligned.out.bam")
    params:
        index = INDEX_DIR,
        prefix = "results/{sample}."
    threads: 6
    singularity:
        "shub://jdwheaton/singularity-ngs:star_htseq"
    log:
        "logs/{sample}.star.out"
    shell:
        '''
        STAR \
        --genomeDir {params.index} \
        --runThreadN {threads} \
        --outFilterMultimapNmax 1 \
        --outSAMattributes Standard \
        --outSAMstrandField intronMotif \
        --outSAMtype BAM Unsorted \
        --alignSJoverhangMin 500 \
        --readFilesIn {input} \
        --outFileNamePrefix {params.prefix} \
        &> {log}
        '''

rule samtools_sort:
    input:
        "alignment/{sample}.Aligned.out.bam"
    output:
        "alignment/{sample}.Aligned.out.sorted.bam"
    threads: 3
    singularity:
        "shub://jdwheaton/singularity-ngs:latest"
    log:
        "logs/{sample}.sort.out"
    shell:
        "samtools sort --threads {threads} -o {output} {input} &> {log}"

rule samtools_index:
    input:
        "alignment/{sample}.Aligned.out.sorted.bam"
    output:
        "alignment/{sample}.Aligned.out.sorted.bai"
    threads: 1
    singularity:
        "shub://jdwheaton/singularity-ngs:latest"
    log:
        "logs/{sample}.index.out"
    shell:
        "samtools index {input} {output} &> {log}"

rule featureCounts:
    input:
        expand("alignment/{sample}.Aligned.out.sorted.bam", sample=SAMPLES)
    output:
        COUNT_FILENAME
    threads: 4
    params:
        strand = [1 if STRANDED == 'forward' else 2 if STRANDED == 'reverse' else 0]
    singularity:
        "shub://jdwheaton/singularity-ngs:chip_atac_post"
    log:
        "logs/featureCounts.log"
    shell:
        "featureCounts -s {params.strand} -T {threads} \
        -a {ANNOTATION} \
        -o {output} \
        {input} \
        &> {log}"

rule bam_coverage:
    input:
        BAM = "results/{sample}.Aligned.out.sorted.bam",
        BAI = "results/{sample}.Aligned.out.sorted.bai"
    output:
        "bigwigs/{sample}.rpkm.bw"
    singularity:
        "shub://jdwheaton/singularity-ngs:chip_atac_post"
    threads: 4
    log:
        "logs/{sample}_bamcoverage.log"
    shell: """
            bamCoverage -p {threads} -b {input.BAM} -o {output} \
                --binSize 10 \
                --normalizeUsing RPKM \
                --ignoreForNormalization chrX chrY chrM \
                &> {log}
                """

rule mutltiqc:
    input:
        ALL_BAMCOV + ALL_FASTQC + COUNT_FILENAME
    output:
        "multiqc_report.html"
    singularity:
        "docker://ewels/multiqc"
    threads: 1
    shell:
        "multiqc ."