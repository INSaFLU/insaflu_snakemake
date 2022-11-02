rule snpeff_concat:
    input:
        expand(
            "projects/{project}/sample_{sample}/freebayes/{sample}_snpeff.vcf",
            sample=config_user["samples"],
            project=config_user["project"],
        ),
    output:
        "projects/{project}/main_result/validated_minor_iSNVs.csv",
    shell:
        "python utils/validated_minor_iSNVs.py '{input}' {output}"


rule snpeff_concat_indels:
    input:
        expand(
            "projects/{project}/sample_{sample}/freebayes/{sample}_snpeff.vcf",
            sample=config_user["samples"],
            project=config_user["project"],
        ),
    output:
        "projects/{project}/main_result/validated_minor_iSNVs_inc_indels.csv",
    shell:
        "python utils/minor_iSNVs_inc_indels.py '{input}' {output}"


rule proportions_iSNVs_graph:
    input:
        expand(
            "projects/{project}/sample_{sample}/freebayes/{sample}_snpeff.vcf",
            sample=config_user["samples"],
            project=config_user["project"],
        ),
    output:
        out_file="projects/{project}/main_result/proportions_iSNVs_graph.csv",
    shell:
        "python utils/proportions_iSNVs_graph.py '{input}' {output.out_file}"
