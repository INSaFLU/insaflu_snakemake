def get_raw_input_fastq_se(wildcards):
    return f"user_data/{config_user['samples'][wildcards.sample]['fastq1']}.fastq.gz"

rule raw_fastqc_se:
    input:
        get_raw_input_fastq_se
    output:
        o = "samples/{sample}/raw_fastqc/{sample}_fastqc.html",
        dir=directory("samples/{sample}/raw_fastqc")
    conda:
        "../envs/fastqc.yaml"
    params:
        "--nogroup"
    shell:
        "fastqc {input} -o {output.dir} {params} && python3 utils/fix_fastq_output.py {wildcards.sample}"


def get_raw_input_fastq(wildcards):
    return [f"user_data/{config_user['samples'][wildcards.sample]['fastq1']}.fastq.gz",
            f"user_data/{config_user['samples'][wildcards.sample]['fastq2']}.fastq.gz"]
rule raw_fastqc_pe:
    input:
        get_raw_input_fastq
    output:
        o1 = "samples/{sample}/raw_fastqc/{sample}_1_fastqc.html",
        o2 = "samples/{sample}/raw_fastqc/{sample}_2_fastqc.html",
        dir=directory("samples/{sample}/raw_fastqc")
    conda:
        "../envs/fastqc.yaml"
    params:
        "--nogroup"
    shell:
        "fastqc {input} -o {output.dir} {params} && python3 utils/fix_fastq_output.py {wildcards.sample}"

rule trimmed_fastqc_pe:
    input:
        i1 = "samples/{sample}/trimmed_reads/{sample}_1.trimmed.fastq.gz",
        i2 = "samples/{sample}/trimmed_reads/{sample}_2.trimmed.fastq.gz"
    output:
        o1 = "samples/{sample}/trimmed_fastqc/{sample}_1.trimmed_fastqc.html",
        o2 = "samples/{sample}/trimmed_fastqc/{sample}_2.trimmed_fastqc.html",
        dir=directory("samples/{sample}/trimmed_fastqc")
    conda:
        "../envs/fastqc.yaml"
    params:
        "--nogroup"
    shell:
        "fastqc {input.i1} {input.i2} -o {output.dir} {params}"       

rule trimmed_fastqc_se:
    input:
        "samples/{sample}/trimmed_reads/{sample}.trimmed.fastq.gz"
    output:
        o="samples/{sample}/trimmed_fastqc/{sample}.trimmed_fastqc.html",
        dir=directory("samples/{sample}/trimmed_fastqc/")
    conda:
        "../envs/fastqc.yaml"
    params:
        "--nogroup"
    shell:
        "fastqc {input} -o {output.dir} {params}"

