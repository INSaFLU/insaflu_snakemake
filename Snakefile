from utils.get_gene_bank import get_genes
from utils.import_user_data import Data
from utils.import_user_data import read_yaml
from utils.get_locus import get_locus_protein, get_locus
import yaml
import csv
from Bio import SeqIO

class Checkpoint_MakePattern:
    def __init__(self,prefix,sufix,genbank_file,possible_name, coverage_file):
        self.prefix = prefix
        self.sufix = sufix
        self.genbank_file = genbank_file
        self.possible_name = possible_name
        self.coverage_file = coverage_file

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'check_csv'; this will trigger an
        # exception until that rule has been run.
        checkpoints.mergeCoverage.get(**w)

        # the magic, such as it is, happens here: we create the
        # information used to expand the pattern, using arbitrary
        # Python code.

        pattern = self.get_output(self.genbank_file,self.possible_name,self.coverage_file)
        return pattern
    
    def get_output(self,genbank_file,possible_name, coverage_file):
        locus_protein = []
        valide_locus = self.get_locus_w_coverage(coverage_file)
        handle_gb = open(genbank_file)
        for record in SeqIO.parse(handle_gb, "genbank"):
            for features in record.features:
                if (features.type == 'CDS'):
                    try:
                        a = features.qualifiers["locus_tag"]
                        if int(features.qualifiers['locus_tag'][0][-5:]) in valide_locus:
                            locus_protein.append(f"{self.prefix}{int(features.qualifiers['locus_tag'][0][-5:])}/Alignment_aa_{int(features.qualifiers['locus_tag'][0][-5:])}_{features.qualifiers['gene'][0]}{self.sufix}")
                    except:
                        locus_protein.append(f"{self.prefix}{possible_name}/Alignment_aa_{possible_name}_{features.qualifiers['gene'][0]}{self.sufix}")
        handle_gb.close()
        return locus_protein


    def get_locus_w_coverage(self,coverage_file):
        with open(coverage_file, newline='') as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=',')
            coverage_list = list(csv_reader)
            size = len(coverage_list[0])
            final_output = []
            for times in range(1,size):
                for sample in coverage_list:
                    print(float(sample[times]))
                    if float(sample[times]) >= 90.0:
                        final_output.append(times)
                        break
        return final_output



sample_data = Data("./config_user/sample_info_1.csv")
# sample_data = Data("./config_user/flu.csv")

run_config = read_yaml('./config_user/config_user1.yaml')
# run_config = read_yaml('./config_user/config_user2.yaml')

locus_protein_alignment = get_locus_protein(run_config["gb_reference"],run_config["locus"])

def get_output_sample_se():
    return(
        expand("samples/{sample}/raw_fastqc/{sample}_fastqc.html", sample=config_user['samples']),
        expand("samples/{sample}/trimmed_fastqc/{sample}.trimmed_fastqc.html", sample=config_user['samples']),
        expand("samples/{sample}/abricate/abricate_{sample}.csv", sample=config_user['samples']),
    )

def get_output_sample_pe():
    return(
        expand("samples/{sample}/raw_fastqc/{sample}_{direction}_fastqc.html", sample=config_user['samples'],direction=["1","2"]), #generalizar
        expand("samples/{sample}/trimmed_fastqc/{sample}_{direction}.trimmed_fastqc.html", sample=config_user['samples'],direction=["1","2"]),
        expand("samples/{sample}/abricate/abricate_{sample}.csv", sample=config_user['samples']),
        )    

def get_output_files_se():
    return(
        expand("samples/{sample}/raw_fastqc/{sample}_fastqc.html", sample=config_user['samples']),
        expand("samples/{sample}/trimmed_fastqc/{sample}.trimmed_fastqc.html", sample=config_user['samples']),
        expand("samples/{sample}/spades/contigs.fasta", sample=config_user['samples']),
        expand("samples/{sample}/abricate/abricate_{sample}.csv", sample=config_user['samples']),
        #expand("align_samples/{sample}/snippy/snps.consensus.fa",project=config_user['project'], sample=config_user['samples']),
        #expand("projects/{project}/main_result/consensus/{sample}__SARS_COV_2_consensus.fasta", sample=config_user['samples'], project=config_user['project']),
        #expand("projects/{project}/main_result/AllConsensus.fasta", project=config_user['project']),
        #expand("projects/{project}/main_result/coverage/{sample}_coverage.csv", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/coverage.csv",project=config_user['project']),
        #expand("projects/{project}/main_result/freebayes/{sample}_var.vcf", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",project=config_user['project'], seg=get_locus(run_config["gb_reference"],run_config["locus"])),
        #expand("projects/{project}/main_result/{seg}/Alignment_nt_{seg}_masked.fasta",project=config_user['project'], seg=get_locus(run_config["gb_reference"],run_config["locus"])),
        # expand("projects/{project}/main_result/depth/{sample}__{ref}.depth",sample=config_user['samples'], project=config_user['project'], ref=get_locus(run_config["gb_reference"],run_config["locus"])),        
        # expand("projects/{project}/main_result/depth/{sample}.depth",sample=config_user['samples'], project=config_user['project']),        
        expand("projects/{project}/main_result/snpeff/{sample}_snpeff.vcf",sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/mafft/Alignment_nt_All_concat.fasta", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/{locus_protein_alignment_file}_trans.fasta", project=config_user['project'],locus_protein_alignment_file = locus_protein_alignment),
        expand("projects/{project}/main_result/{locus_protein_alignment_file}_mafft.fasta", project=config_user['project'],locus_protein_alignment_file = locus_protein_alignment),
        expand("projects/{project}/main_result/fasttre/tree", sample=config_user['samples'], project=config_user['project']), 
        expand("projects/{project}/main_result/{locus_protein_alignment_file}_tree.tree", project=config_user['project'],locus_protein_alignment_file = locus_protein_alignment),
    )

def get_output_files_pe():
    return(
        expand("samples/{sample}/raw_fastqc/{sample}_{direction}_fastqc.html", sample=config_user['samples'],direction=["1","2"]), #generalizar
        expand("samples/{sample}/trimmed_fastqc/{sample}_{direction}.trimmed_fastqc.html", sample=config_user['samples'],direction=["1","2"]),
        # expand("samples/{sample}/spades/contigs.fasta", sample=config_user['samples']),
        # expand("samples/{sample}/abricate/abricate_{sample}.csv", sample=config_user['samples']),
        #expand("align_samples/{sample}/snippy/snps.consensus.fa",project=config_user['project'], sample=config_user['samples']),
        #expand("projects/{project}/main_result/consensus/{sample}__SARS_COV_2_consensus.fasta", sample=config_user['samples'], project=config_user['project']),
        #expand("projects/{project}/main_result/AllConsensus.fasta", project=config_user['project']),
        #expand("projects/{project}/main_result/coverage/{sample}_coverage.csv", sample=config_user['samples'], project=config_user['project']),
        expand("projects/{project}/main_result/coverage.csv",project=config_user['project']),
        expand("projects/{project}/main_result/coverage_translate.csv",project=config_user['project']),
        # expand("projects/{project}/main_result/snpeff_samples/{sample}_snpeff.vcf",sample=config_user['samples'], project=config_user['project'], ref=get_locus(run_config["gb_reference"],run_config["locus"])),        

        #expand("projects/{project}/main_result/freebayes/{sample}_var.vcf", sample=config_user['samples'], project=config_user['project']),
        #expand("projects/{project}/main_result/{seg}/Alignment_nt_{seg}.fasta",project=config_user['project'], seg=get_locus(run_config["gb_reference"],run_config["locus"])),
        #expand("projects/{project}/main_result/{seg}/Alignment_nt_{seg}_masked.fasta",project=config_user['project'], seg=get_locus(run_config["gb_reference"],run_config["locus"])),
        expand("projects/{project}/main_result/depth/{sample}__{ref}.depth",sample=config_user['samples'], project=config_user['project'], ref=get_locus(run_config["gb_reference"],run_config["locus"])),        
        expand("projects/{project}/main_result/depth/{sample}.depth",sample=config_user['samples'], project=config_user['project']),        
        # expand("projects/{project}/main_result/snpeff_samples/{sample}_snpeff.vcf",sample=config_user['samples'], project=config_user['project']),
        # expand("projects/{project}/main_result/validated_minor_iSNVs.csv",project=config_user['project']),
        # expand("projects/{project}/main_result/validated_variants.csv",project=config_user['project']),
        # expand("projects/{project}/main_result/validated_minor_iSNVs_inc_indels.csv",project=config_user['project']),
        # expand("projects/{project}/main_result/proportions_iSNVs_graph.csv",project=config_user['project']),
        expand("projects/{project}/main_result/mafft/Alignment_nt_All_concat.fasta", sample=config_user['samples'], project=config_user['project']),
        Checkpoint_MakePattern(f'projects/{run_config["project_name"]}/main_result/',"_trans.fasta",run_config['gb_reference'],run_config["locus"],f"projects/{config_user['project']}/main_result/coverage_translate.csv"),
        Checkpoint_MakePattern(f'projects/{run_config["project_name"]}/main_result/',"_mafft.fasta",run_config['gb_reference'],run_config["locus"],f"projects/{config_user['project']}/main_result/coverage_translate.csv"),
        Checkpoint_MakePattern(f'projects/{run_config["project_name"]}/main_result/',"_tree.tree",run_config['gb_reference'],run_config["locus"], f"projects/{config_user['project']}/main_result/coverage_translate.csv"),
        #expand("projects/{project}/main_result/snpeff/ready.txt",project=config_user['project']),

        # expand("projects/{project}/main_result/{locus_protein_alignment_file}_trans.fasta", project=config_user['project'],locus_protein_alignment_file = locus_protein_alignment),
        # expand("projects/{project}/main_result/{locus_protein_alignment_file}_mafft.fasta", project=config_user['project'],locus_protein_alignment_file = locus_protein_alignment),
        # expand("projects/{project}/main_result/fasttre/tree", sample=config_user['samples'], project=config_user['project']), 
        # expand("projects/{project}/main_result/{locus_protein_alignment_file}_tree.tree", project=config_user['project'],locus_protein_alignment_file = locus_protein_alignment),

    ) 


def prepare_run(settings):
    if settings['only_samples'] == True:
        if sample_data.get_sample_2() == []:
            return get_output_sample_se
        else:
            return get_output_sample_pe
    else:
        if sample_data.get_sample_2() == []:
            return get_output_files_se
        else:
            return get_output_files_pe


REFERENCE_GB =run_config['gb_reference'] 
REFERENCE  =run_config['fasta_reference']
x = re.findall("(?<=/)(.*?)(?=.fasta)",REFERENCE)
REFERENCE_NAME = x[0]


get_output = prepare_run(run_config)
config_user = {'samples':sample_data.get_sample_names(), 
               'project':run_config['project_name'], 
               'locus': run_config['locus'], 
               'proteins':get_genes(run_config['gb_reference'])}


with open('config/config_run.yaml', 'w') as file:
    documents = yaml.dump(config_user, file)
    file.close()

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
include: "rules/freebays_concat.smk"
include: "rules/snpeff_sample.smk"
include: "rules/snp_variant_validated.smk"

rule all:
    input:
        get_output()