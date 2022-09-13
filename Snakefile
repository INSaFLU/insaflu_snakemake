from utils.get_gene_bank import get_genes
from utils.import_user_data import Data
import yaml
import os 

user_data = Data("./config_user/sample_info_1.csv")
REFERENCE = "reference/SARS_CoV_2_Wuhan_Hu_1_MN908947.fasta"
REFERENCE_GB = "reference/SARS_CoV_2_Wuhan_Hu_1_MN908947.gb"
REFERENCE_GFF3 = "reference/SARS_CoV_2_Wuhan_Hu_1_MN908947.gff3"



def get_output_files_se(SAMPLES, PROJECT,REFERENCE):
    SAMPLES = SAMPLES
    PROJECT = PROJECT
    config_user = {'samples':SAMPLES, 'project':PROJECT, 'ref':'SARS_CoV_2', 'proteins':get_genes(REFERENCE_GB)}
    with open('config/config_run.yaml', 'w') as file:
        documents = yaml.dump(config_user, file)
    return(
        expand("samples/{sample}/raw_fastqc/{sample}_fastqc.html", sample=SAMPLES),
        expand("samples/{sample}/trimmed_fastqc/{sample}.trimmed_fastqc.html", sample=SAMPLES),
        #expand("samples/{sample}/spades/contigs.fasta", sample=SAMPLES),
        expand("samples/{sample}/abricate/abricate_{sample}.csv", sample=SAMPLES),
        #expand("align_samples/{sample}/snippy/snps.consensus.fa",project=PROJECT, sample=SAMPLES),
        #expand("projects/{project}/main_result/consensus/{sample}__SARS_COV_2_consensus.fasta", sample=SAMPLES, project=PROJECT),
        #expand("projects/{project}/main_result/AllConsensus.fasta", project=PROJECT),
        #expand("projects/{project}/main_result/coverage/{sample}_coverage.tab", sample=SAMPLES, project=PROJECT),
        expand("projects/{project}/main_result/coverage.csv",project=PROJECT),
        #expand("projects/{project}/main_result/freebayes/{sample}_var.vcf", sample=SAMPLES, project=PROJECT),
        expand("projects/{project}/main_result/snpeff/{sample}_snpeff.vcf",sample=SAMPLES, project=PROJECT),
        expand("projects/{project}/main_result/mafft/mafft.fasta", sample=SAMPLES, project=PROJECT),
        expand("projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}.fasta",project=PROJECT,ref=config_user["ref"],protein=config_user["proteins"]),
        expand("projects/{project}/main_result/fasttre/tree", sample=SAMPLES, project=PROJECT), 
        expand("projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}_tree.tree", project=PROJECT,ref=config_user["ref"],protein=config_user["proteins"]),

    )

def get_output_files_pe(SAMPLES, PROJECT,REFERENCE):
    SAMPLES = SAMPLES
    PROJECT = PROJECT
    config_user = {'samples':SAMPLES, 'project':PROJECT, 'ref':'SARS_CoV_2', 'proteins':get_genes(REFERENCE_GB)}
    with open('config/config_run.yaml', 'w') as file:
        documents = yaml.dump(config_user, file)
    return(
        expand("samples/{sample}/raw_fastqc/{sample}_{direction}_fastqc.html", sample=SAMPLES,direction=["1","2"]), #generalizar
        expand("samples/{sample}/trimmed_fastqc/{sample}_{direction}.trimmed_fastqc.html", sample=SAMPLES,direction=["1","2"]),
        #expand("samples/{sample}/spades/contigs.fasta", sample=SAMPLES),
        expand("samples/{sample}/abricate/abricate_{sample}.csv", sample=SAMPLES),
        #expand("align_samples/{sample}/snippy/snps.consensus.fa",project=PROJECT, sample=SAMPLES),
        #expand("projects/{project}/main_result/consensus/{sample}__SARS_COV_2_consensus.fasta", sample=SAMPLES, project=PROJECT),
        #expand("projects/{project}/main_result/AllConsensus.fasta", project=PROJECT),
        #expand("projects/{project}/main_result/coverage/{sample}_coverage.csv", sample=SAMPLES, project=PROJECT),
        expand("projects/{project}/main_result/coverage.csv",project=PROJECT),
        #expand("projects/{project}/main_result/freebayes/{sample}_var.vcf", sample=SAMPLES, project=PROJECT),
        expand("projects/{project}/main_result/snpeff/{sample}_snpeff.vcf",sample=SAMPLES, project=PROJECT),
        expand("projects/{project}/main_result/mafft/mafft.fasta", sample=SAMPLES, project=PROJECT),
        expand("projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}.fasta",project=PROJECT,ref=config_user["ref"],protein=config_user["proteins"]),
        expand("projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}_mafft.fasta",project=PROJECT,ref=config_user["ref"],protein=config_user["proteins"]),
        expand("projects/{project}/main_result/fasttre/tree", sample=SAMPLES, project=PROJECT), 
        expand("projects/{project}/main_result/{ref}/Alignment_aa_{ref}_{protein}_tree.tree", project=PROJECT,ref=config_user["ref"],protein=config_user["proteins"]),
        
    )    




if user_data.get_sample_2() == []:
    get_output = get_output_files_se
else:
    get_output = get_output_files_pe

include: "rules/fastqc.smk"
include: "rules/trimmomatic.smk"
include: "rules/spades.smk"
include: "rules/abricate.smk"
include: "rules/snippy.smk"
include: "rules/makeproject.smk"
include: "rules/getCoverage.smk"
include: "rules/mergeCoverage.smk"
include: "rules/freebayes.smk"
include: "rules/snpeff.smk"
include: "rules/concat.smk"
include: "rules/mafft.smk"
include: "rules/translate.smk"
include: "rules/move_depth.smk"
include: "rules/msa_masker.smk"
include: "rules/fasttree.smk"
include: "rules/seqret.smk"
include: "rules/mafft_proteins.smk"
include: "rules/fastree_proteins.smk"



rule all:
    input:
        get_output(user_data.get_sample_names(),"insaflu_bigger",REFERENCE)

