configfile: "config.yaml"

SAMPLES, = glob_wildcards("raw_data/{smp}_R1.fastq.gz")
INDEX_DIR = config["star_index"]
ANNOTATION = config["annotation"]
STRANDED = config["stranded"]

COUNT_FILENAME = [config["count_filename"]]

ALL_FASTQC = expand("fastqc_out/{sample}_R1_fastqc.zip", sample=SAMPLES)
ALL_BAMCOV = expand("results/{sample}.rpkm.bw", sample=SAMPLES)

rule all:
    input: ALL_BAMCOV + ALL_FASTQC + COUNT_FILENAME + ["multiqc_report.html"]

rule fastqc:
    input:
        "raw_data/{sample}_R1.fastq.gz",
        "raw_data/{sample}_R2.fastq.gz"
    output:
        "fastqc_out/{sample}_R1_fastqc.html",
        "fastqc_out/{sample}_R1_fastqc.zip",
        "fastqc_out/{sample}_R2_fastqc.html",
        "fastqc_out/{sample}_R2_fastqc.zip"
    log:
        "logs/{sample}_fastqc.log"
    singularity:
        "docker://biocontainers/fastqc"
    threads: 1
    shell:
        "fastqc -o fastqc_out/ {input} &> {log}"

rule trim_galore:
    input:
        "raw_data/{sample}_R1.fastq.gz",
        "raw_data/{sample}_R2.fastq.gz"
    output:
        "trimmed_fastq/{sample}_R1_val_1.fq",
        "trimmed_fastq/{sample}_R2_val_2.fq"
    log:
        "logs/{sample}_trimgalore.log"
    threads: 1
    singularity:
        "docker://dukegcb/trim-galore"
    shell:
        "trim_galore --dont_gzip --fastqc -o trimmed_fastq/ --paired {input} &> {log}"

rule star_align:
    input:
        "trimmed_fastq/{sample}_R1_val_1.fq",
        "trimmed_fastq/{sample}_R2_val_2.fq"
    output:
        temp("results/{sample}.Aligned.out.bam")
    params:
        index = INDEX_DIR,
        prefix = "results/{sample}."
    threads: 6
    singularity:
        "docker://alexdobin/star"
    log:
        "logs/{sample}.star.log"
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
        "results/{sample}.Aligned.out.bam"
    output:
        "results/{sample}.Aligned.out.sorted.bam"
    threads: 3
    singularity:
        "docker://biocontainers/samtools"
    log:
        "logs/{sample}.sort.out"
    shell:
        "samtools sort --threads {threads} -o {output} {input} &> {log}"

rule samtools_index:
    input:
        "results/{sample}.Aligned.out.sorted.bam"
    output:
        "results/{sample}.Aligned.out.sorted.bai"
    threads: 1
    singularity:
        "docker://biocontainers/samtools"
    log:
        "logs/{sample}.index.out"
    shell:
        "samtools index {input} {output} &> {log}"

rule featureCounts:
    input:
        expand("results/{sample}.Aligned.out.sorted.bam", sample=SAMPLES)
    output:
        COUNT_FILENAME
    threads: 4
    params:
        strand = [1 if STRANDED == 'forward' else 2 if STRANDED == 'reverse' else 0]
    singularity:
        "docker://genomicpariscentre/deeptools"
    log:
        "logs/featureCounts.log"
    shell:
        "featureCounts -p -s {params.strand} -T {threads} \
        -a {ANNOTATION} \
        -o {output} \
        {input} \
        &> {log}"

rule bam_coverage:
    input:
        BAM = "results/{sample}.Aligned.out.sorted.bam",
        BAI = "results/{sample}.Aligned.out.sorted.bai"
    output:
        "results/{sample}.rpkm.bw"
    singularity:
        "docker://dukegcb/deeptools"
    threads: 4
    log:
        "logs/{sample}_bamcoverage.log"
    shell: """
            bamCoverage -p {threads} -b {input.BAM} -o {output} \
                --binSize 10 \
                --normalizeUsing RPKM \
                --ignoreForNormalization chrX chrY chrM \
                --extendReads \
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
