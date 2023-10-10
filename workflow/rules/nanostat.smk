
rule raw_nanostat:
    input:
        get_raw_input_ont,
    output:
        dir=directory("samples/{sample}/raw_nanostat/"),
        stats="samples/{sample}/raw_nanostat/{sample}_stats.txt",
    conda:
        "../envs/nanostat.yaml"
    resources:
        mem_mb=memory["raw_nanostat"],
    log:
        "logs/samples/{sample}/raw_nanostat/{sample}.log",
    benchmark:
        "benchmark/samples/{sample}/raw_nanostat/{sample}.tsv"
    threads: config["nanostat_threads"]
    shell:
        "mkdir {output.dir} -p && NanoStat --fastq {input} --outdir {output.dir} -n {wildcards.sample}_stats.txt -t {threads}"


rule trimmed_nanostat:
    input:
        "samples/{sample}/trimmed_reads/nano_{sample}.trimmed.fastq.gz",
    output:
        dir=directory("samples/{sample}/nano_trimmed_fastqc/"),
        stats="samples/{sample}/nano_trimmed_fastqc/{sample}_stats.txt",
    conda:
        "../envs/nanostat.yaml"
    resources:
        mem_mb=memory["trimmed_nanostat"],
    log:
        "logs/samples/{sample}/nano_trimmed_fastqc/{sample}.log",
    benchmark:
        "benchmark/samples/{sample}/nano_trimmed_fastqc/{sample}.tsv"
    threads: config["nanostat_threads"]
    shell:
        "NanoStat --fastq {input} --outdir {output.dir} -n {wildcards.sample}_stats.txt -t {threads}"


rule raw_rabbit_qc:
    input:
        get_raw_input_ont,
    output:
        dir=directory("samples/{sample}/raw_rabbit_qc/"),
        stats="samples/{sample}/raw_rabbit_qc/{sample}_stats.html",
    conda:
        "../envs/rabbitqc.yaml"
    resources:
        mem_mb=memory["raw_nanostat"],
    log:
        "logs/samples/{sample}/raw_rabbit_qc/{sample}.log",
    benchmark:
        "benchmark/samples/{sample}/raw_rabbit_qc/{sample}.tsv"
    params:
        "-w 3 -D"
    shell:
       "../workflow/software/RabbitQC/rabbit_qc {params} -i {input} -h {output.stats}"


rule processed_rabbit_qc:
    input:
        "samples/{sample}/trimmed_reads/nano_{sample}.trimmed.fastq.gz",
    output:
        dir=directory("samples/{sample}/processed_rabbit_qc/"),
        stats="samples/{sample}/processed_rabbit_qc/{sample}_stats.html",
    conda:
        "../envs/rabbitqc.yaml"
    resources:
        mem_mb=memory["trimmed_nanostat"],
    log:
        "logs/samples/{sample}/nano_trimmed_fastqc/{sample}.log",
    benchmark:
        "benchmark/samples/{sample}/nano_trimmed_fastqc/{sample}.tsv"
    params:
        "-w 3 -D"
    shell:
       "../workflow/software/RabbitQC/rabbit_qc {params} -i {input} -h {output.stats}"
