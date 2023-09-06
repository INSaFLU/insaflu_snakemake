localrules:initiate_folder,makeproject,assemble_consensus,create_segments,
rule initiate_folder:
    conda:
        "../envs/base.yaml"
    output:
        dir("projects/{project}/main_result/"),
    resources:
        mem_mb=memory["initiate_folder"],
    log:
        expand("logs/{project}/main_result/initiate_folder.log", project=PROJECT_NAME),
    benchmark:
        f"benchmark/{PROJECT_NAME}/main_result/initiate_folder.tsv"
    shell:
        "mkdir {output} && cp {software_parameters_path} projects/{wildcards.project}/"


rule makeproject:
    input:
        get_consensus_project,
    output:
        depth="projects/{project}/sample_{sample}/snps.depth.gz",
        vcf="projects/{project}/sample_{sample}/snps.vcf",
        consensus="projects/{project}/sample_{sample}/{sample}_consensus.fasta",
    conda:
        "../envs/base.yaml"
    params:
        get_directory,
    resources:
        mem_mb=memory["makeproject"],
    log:
        "logs/projects/{project}/makeproject/{sample}.log",
    benchmark:
        "benchmark/projects/{project}/makeproject/{sample}.tsv"
    shell:
        "mkdir projects/{wildcards.project}/sample_{wildcards.sample}/ -p && "
        " cp -r {params} projects/{wildcards.project}/sample_{wildcards.sample}/ "


rule assemble_consensus:
    input:
        every_consensus=expand(
            "projects/{project}/sample_{sample}/{sample}_consensus.fasta",
            project=PROJECT_NAME,
            sample=config_user["samples"],
        ),
        coverage="projects/{project}/main_result/coverage_translate.csv",
    output:
        AllConsensus="projects/{project}/main_result/AllConsensus.fasta",
        all_consensus_no_ref="projects/{project}/main_result/AllConsensus_no_ref.fasta",
        All_nt="projects/{project}/main_result/All_nt.fasta",
        All_nt_only_90plus="projects/{project}/main_result/All_nt_only_90plus.fasta",
    conda:
        "../envs/base.yaml"
    resources:
        mem_mb=memory["assemble_consensus"],
    log:
        "logs/projects/{project}/main_result/assemble_consensus.log",
    benchmark:
        "benchmark/projects/{project}/main_result/assemble_consensus.tsv"
    shell:
        "python {scripts_directory}generate_AllConsensus.py {input.coverage} {REFERENCE_GB} '{input.every_consensus}' {REFERENCE_FASTA} {output.AllConsensus} {output.all_consensus_no_ref} "
        "&& python {scripts_directory}concat_segments.py '{input.every_consensus}' {REFERENCE_GB} {output.All_nt} {input.coverage} {REFERENCE_FASTA} {output.All_nt_only_90plus}"


rule create_segments:
    input:
        actual_input=expand(
            "projects/{project}/main_result/Alignment_nt_All_sep.fasta",
            project=PROJECT_NAME,
        ),
    output:
        expand(
            "projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",
            project=PROJECT_NAME,
            seg=SEGMENTS,
        ),
    resources:
        mem_mb=memory["create_segments"],
    log:
        expand(
            "logs/projects/{project}/main_result/create_segments/{seg}.log",
            project=PROJECT_NAME,
            seg=SEGMENTS,
        ),
    benchmark:
        f"benchmark/projects/{PROJECT_NAME}/main_result/create_segments/{SEGMENTS}.tsv"
    conda:
        "../envs/base.yaml"
    shell:
        "python {scripts_directory}split_files_by_locus.py {input} projects/{PROJECT_NAME}/main_result {REFERENCE_GB}"
