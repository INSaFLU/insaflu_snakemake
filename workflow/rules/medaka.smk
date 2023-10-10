localrules:medaka_split_depth,mask_between_top_and_50,create_align_file,msa_masker_medaka,get_masked_consensus_medaka,mask_regions_consensus_medaka,
def get_model():
    params =  f"-m {software_parameters['Medaka_model']}"
    params +=  f" -p ../../../{PRIMER_FASTA}" if PRIMER_FASTA else ""
    return params

rule medaka_consensus:
    input:
        i="samples/{sample}/trimmed_reads/nano_{sample}.trimmed.fastq.gz",
        ref=REFERENCE_FASTA,
    output:
        out="align_samples/{sample}/medaka/consensus.fasta",
        i="align_samples/{sample}/medaka/calls_to_draft.bam",
        i2="align_samples/{sample}/medaka/snps.bam",
        hdf="align_samples/{sample}/medaka/consensus_probs.hdf",
        ref_out="align_samples/{sample}/medaka/ref.fasta",
    conda:
        "../envs/medaka_1_4_4.yaml"
    threads: 2
    resources:
        mem_mb=memory["medaka_consensus"],
    log:
        "logs/align_samples/{sample}/medaka/medaka_consensus.log",
    benchmark:
        "benchmark/align_samples/{sample}/medaka/medaka_consensus.tsv"
    params:
        get_model()
    shell:
        "rm -r align_samples/{wildcards.sample}/medaka/ && "
        "mkdir align_samples/{wildcards.sample}/medaka/ && "
        "cp {input.ref} {output.ref_out} && "
        "medaka_consensus -i {input.i} -d {output.ref_out} "
        "-o align_samples/{wildcards.sample}/medaka -t {threads} {params}"
        " && cp {output.i} {output.i2}"


rule medaka_get_depth:
    input:
        i="align_samples/{sample}/medaka/calls_to_draft.bam",
    output:
        only_depth="align_samples/{sample}/medaka/snps.depth.gz",
        depth="align_samples/{sample}/medaka/snps.depth",
        ind="align_samples/{sample}/medaka/snps.depth.gz.tbi",
    conda:
        "../envs/medaka_1_4_4.yaml"
    params:
        "-aa ",
    resources:
        mem_mb=memory["medaka_get_depth"],
    log:
        "logs/align_samples/{sample}/medaka/medaka_get_depth.log",
    benchmark:
        "benchmark/align_samples/{sample}/medaka/medaka_get_depth.tsv"
    shell:
        "samtools depth {params} {input.i} | bgzip -c > {output.only_depth} "
        "&& tabix -p vcf {output.only_depth} && gunzip -c {output.only_depth}  > {output.depth} "


rule medaka_split_depth:
    input:
        "align_samples/{sample}/medaka/snps.depth",
    output:
        "align_samples/{sample}/medaka/{seg}.depth",
    conda:
        "../envs/base.yaml"
    resources:
        mem_mb=memory["medaka_split_depth"],
    log:
        "logs/align_samples/{sample}/medaka/medaka_split_depth/{seg}.log",
    benchmark:
        "benchmark/align_samples/{sample}/medaka/medaka_split_depth/{seg}.tsv"
    shell:
        "python3 {scripts_directory}split_depth_file.py {input} {REFERENCE_GB}"


rule medaka_call_variants:
    input:
        hdf="align_samples/{sample}/medaka/consensus_probs.hdf",
        ref="align_samples/{sample}/medaka/ref.fasta",
        consensus="align_samples/{sample}/medaka/consensus.fasta",
    output:
        vcf="align_samples/{sample}/medaka/round_1.vcf",
    conda:
        "../envs/medaka_1_4_4.yaml"
    resources:
        mem_mb=memory["medaka_call_variants"],
    log:
        "logs/align_samples/{sample}/medaka/medaka_call_variants.log",
    benchmark:
        "benchmark/align_samples/{sample}/medaka/medaka_call_variants.tsv"
    shell:
        "medaka variant --verbose {input.ref} {input.hdf} {output.vcf}"


rule medaka_annotate_variants:
    input:
        ref="align_samples/{sample}/medaka/ref.fasta",
        vcf="align_samples/{sample}/medaka/round_1.vcf",
        bam="align_samples/{sample}/medaka/calls_to_draft.bam",
    output:
        snps="align_samples/{sample}/medaka/snps_ann.vcf",
        vcf_zipped="align_samples/{sample}/medaka/snps_ann.vcf.gz",
    conda:
        "../envs/medaka_1_2_1.yaml"
    resources:
        mem_mb=memory["medaka_annotate_variants"],
    log:
        "logs/align_samples/{sample}/medaka/medaka_annotate_variants.log",
    benchmark:
        "benchmark/align_samples/{sample}/medaka/medaka_annotate_variants.tsv"
    shell:
        "medaka tools annotate {input.vcf} {input.ref} {input.bam} {output.snps} && "
        "bgzip {output.snps} -c > {output.vcf_zipped} && tabix {output.vcf_zipped}"


rule mask_between_top_and_50:
    input:
        vcf_file="align_samples/{sample}/medaka/snps_ann.vcf",
        file_coverage="align_samples/{sample}/medaka/snps.depth.gz",
        normal_reference_fasta="align_samples/{sample}/medaka/ref.fasta",
    output:
        vcf_file_out="align_samples/{sample}/medaka/snps.vcf",
        vcf_file_out_removed_by_filter="align_samples/{sample}/medaka/snps_removed_by_filter.vcf",
        temp_file=temp("align_samples/{sample}/medaka/temp_file.txt"),
        normal_reference_fasta_removed="align_samples/{sample}/medaka/ref_masked.fasta",
        vcf_file_out_compr="align_samples/{sample}/medaka/snps.vcf.gz",
    conda:
        "../envs/filter_medaka.yaml"
    params:
        get_add_freq_medaka(software_parameters),
    resources:
        mem_mb=memory["mask_between_top_and_50"],
    log:
        "logs/align_samples/{sample}/medaka/filter_medaka.log",
    benchmark:
        "benchmark/align_samples/{sample}/medaka/filter_medaka.tsv"
    shell:
        "touch {output.temp_file} && "
        "python {scripts_directory}add_freq_medaka_consensus.py {input.normal_reference_fasta} {input.vcf_file} {input.file_coverage} "
        "{output.vcf_file_out} {output.vcf_file_out_removed_by_filter} {output.temp_file} {output.normal_reference_fasta_removed} {output.temp_file} {params} &&"
        "bgzip -c -f {output.vcf_file_out} > {output.vcf_file_out_compr} && tabix {output.vcf_file_out_compr}"


rule bcf_consensus:
    input:
        vcf_ziped="align_samples/{sample}/medaka/snps.vcf.gz",
        ref="align_samples/{sample}/medaka/ref_masked.fasta",
    output:
        "align_samples/{sample}/medaka/first_consensus.fasta",
    conda:
        "../envs/bcftools.yaml"
    resources:
        mem_mb=memory["bcf_consensus"],
    log:
        "logs/align_samples/{sample}/medaka/bcf_consensus.log",
    benchmark:
        "benchmark/align_samples/{sample}/medaka/bcf_consensus.tsv"
    shell:
        "bcftools consensus -s SAMPLE -f {input.ref} {input.vcf_ziped} -o {output}"


rule create_align_file:
    input:
        first_consensus="align_samples/{sample}/medaka/first_consensus.fasta",
    output:
        align_file="align_samples/{sample}/medaka/medaka_align_{seg}.fasta",
    conda:
        "../envs/base.yaml"
    resources:
        mem_mb=memory["create_align_file"],
    log:
        "logs/align_samples/{sample}/medaka/create_align_file_{seg}.log",
    benchmark:
        "benchmark/align_samples/{sample}/medaka/create_align_file_{seg}.tsv"
    shell:
        "python {scripts_directory}mask_consensus_by_deep.py {REFERENCE_FASTA} {input.first_consensus} {output.align_file} {wildcards.seg}"


rule msa_masker_medaka:
    input:
        align_file="align_samples/{sample}/medaka/medaka_aligned_{seg}.fasta",
        depth="align_samples/{sample}/medaka/{seg}.depth",
    output:
        temp("align_samples/{sample}/medaka/consensus_aligned_{seg}.fasta"),
    conda:
        "../envs/msa_masker.yaml"
    params:
        "--c " + str(software_parameters["mincov_medaka"] - 1),
    resources:
        mem_mb=memory["msa_masker_medaka"],
    log:
        "logs/align_samples/{sample}/medaka/msa_masker_medaka_{seg}.log",
    benchmark:
        "benchmark/align_samples/{sample}/medaka/msa_masker_medaka_{seg}.tsv"
    shell:
        "python {scripts_directory}msa_masker.py -i {input.align_file} -df {input.depth} -o {output} {params}"


rule get_masked_consensus_medaka:
    input:
        lambda wildcards: expand(
            "align_samples/{sample}/medaka/consensus_aligned_{seg}.fasta",
            sample=wildcards.sample,
            seg=get_locus(REFERENCE_GB),
        ),
    output:
        final_consensus="align_samples/{sample}/medaka/pre_{sample}_consensus.fasta",
    conda:
        "../envs/base.yaml"
    resources:
        mem_mb=memory["get_masked_consensus_medaka"],
    log:
        "logs/align_samples/{sample}/medaka/get_masked_consensus_medaka.log",
    benchmark:
        "benchmark/align_samples/{sample}/medaka/get_masked_consensus_medaka.tsv"
    shell:
        "python {scripts_directory}get_consensus_medaka.py '{input}' {output}"


rule mask_regions_consensus_medaka:
    input:
        consensus="align_samples/{sample}/medaka/pre_{sample}_consensus.fasta",
    output:
        final_consensus="align_samples/{sample}/medaka/{sample}_consensus.fasta",
    conda:
        "../envs/base.yaml"
    params:
        mask_regions_parameters(software_parameters),
    resources:
        mem_mb=memory["mask_regions_consensus_medaka"],
    log:
        "logs/align_samples/{sample}/medaka/mask_regions_consensus_medaka.log",
    benchmark:
        "benchmark/align_samples/{sample}/medaka/mask_regions_consensus_medaka.tsv"
    shell:
        "python {scripts_directory}mask_regions.py {input} {output} {params}"
