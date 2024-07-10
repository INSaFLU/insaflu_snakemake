"""Prepare config file in main_directory."""
import subprocess
import sys
from extract_gb_info import get_id_version, get_locus

CONFIG_FILE = "../workflow/db/snpeff.config"
DATA_DIRECTORY = "../workflow/db/data/"


def prepare_snpeff_run(ref_path_gb, locus, ref_path_fa, reference_name, output):
    """
    The prepare_snpeff_run function creates a snpEFF database for the reference genome.
    It takes as input:
        - The path to the reference genome in GenBank format, and
        - The name of the output file that will be created by this function.

       This function creates a directory called 'config/data' if it does not already exist,
       and then creates three files inside that directory: genes.gbk, sequences.fa,
       and snpeff_genes_prod_snpsiftdb-passed-only-geneidmap
       (the last one is created by running SnpSift).

    :param snpeff_path: Specify the path to snpeff
    :param ref_path_gb: Specify the path to the genbank file of a reference genome
    :param locus: Determine the reference genome to use for snpeff
    :param ref_path_fa: Specify the path to the reference genome fasta file
    :param reference_name: Name the reference files in the snpeff config file
    :param output: Create a file that indicates the database is ready
    :return: A text file with a message &quot;database ready!&quot;
    :doc-author: Trelent
    """
    if len(locus) > 1:
        text = [
            f"{reference_name}.chromosome : {' '.join(str(locus))}\n",
            f"{reference_name}.genome : flu\n",
        ]
        for i in locus:
            text.append(
                f"{reference_name}.{i}.codonTable : Bacterial_and_Plant_Plastid\n"
            )
    else:
        text = [
            f"\n{reference_name}.genome: {reference_name}\n",
            f"{reference_name}.{get_id_version(ref_path_gb)}.codonTable : Bacterial_and_Plant_Plastid\n",
        ]
    with open(CONFIG_FILE, "r", encoding="UTF8") as snpeff_config:
        lines = snpeff_config.readlines()
        joined = "".join(lines)
        find = joined.find("".join(text))

    if find == -1:
        subprocess.run(
            f"mkdir  {DATA_DIRECTORY}{reference_name} -p", shell=True, check=False
        )
        subprocess.run(
            f"cat {ref_path_gb} > {DATA_DIRECTORY}{reference_name}/genes.gbk",
            shell=True,
            check=False,
        )
        subprocess.run(
            f"cat {ref_path_fa} > {DATA_DIRECTORY}{reference_name}/sequences.fa",
            shell=True,
            check=False,
        )
        with open(CONFIG_FILE, "a", encoding="UTF8") as snpeff:
            snpeff.writelines(text)

        subprocess.run(
            f"snpEff build -noCheckProtein -genbank {reference_name} -c {CONFIG_FILE} > {output}",
            shell=True,
            check=False,
        )
    else:
        with open(output, "w", encoding="UTF8") as out:
            out.write("Database Ready!")


if __name__ == "__main__":
    SNPEFF_PATH = sys.argv[1]
    REFERENCE_GB = sys.argv[2]
    REFERENCE_FASTA = sys.argv[3]
    REFERENCE_NAME = sys.argv[4]
    OUTPUT = sys.argv[5]
    SEGMENTS_LIST = get_locus(REFERENCE_GB)
    prepare_snpeff_run(
        REFERENCE_GB, SEGMENTS_LIST, REFERENCE_FASTA, REFERENCE_NAME, OUTPUT
    )
