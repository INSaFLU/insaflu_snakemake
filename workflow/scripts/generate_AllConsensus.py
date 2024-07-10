"""Generate a fasta file with all consensus sequences from samples."""
import sys
import csv
from yaml_io import read_yaml
from Bio import SeqIO
from extract_gb_info import get_locus


def get_locus_out_id(string_id: str):
    return string_id[len(string_id) - string_id[::-1].index("__"):]


def convert_index_to_locus(coverage: dict[str, list[int]], genbank: str) -> dict[str, list[str]]:
    locus_names: list[str] = get_locus(genbank)

    for sample in coverage:
        for idx, value in enumerate(coverage[sample]):
            coverage[sample][idx] = locus_names[value-1]
    return coverage


def generate_AllConsensus(
    coverage_file,
    reference_gb,
    input_files,
    reference_fasta,
    output_w_ref,
    output_no_ref,
):
    """
    The create_consensus_file_for_alignment function takes a list of files, the reference genome in genbank format,
    the species name and an output file name as input. It then creates a fasta file with all the sequences from the
    consensus files that have coverage above 90%. The function also takes into account if it is flu or not. If it is flu
    it will only include samples that have 8 loci covered at 90% or more in one file and in the other it will concatonate all
    segments that have 90% coverage or more.

    :param consensus: Specify the path to the directory containing all of your consensus sequences
    :param reference_gb: Get the name of the reference sequence in a genbank file
    :param species: Determine which reference file to use
    :param output: Name the output file
    :param coverage_file: Specify the file containing all the coverage values for each sample
    :param fasta_file: Get the reference sequence from the fasta file
    :param output_only_90_plus: Create a fasta file with only the sequences that have at least 90% coverage
    :return: A fasta file with the consensus sequences of all samples that have a coverage over 90%
    :doc-author: Trelent
    """
    with open(coverage_file, newline="") as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=",")
        coverage_list = list(csv_reader)

    dic_directory = read_yaml("../config/constants.yaml")

    software_parameters = read_yaml(dic_directory["software_parameters"])

    coverage_value = software_parameters["min_coverage_consensus"]

    for row in coverage_list:
        for value in range(len(row) - 1, 0, -1):
            if float(row[value]) >= coverage_value:
                row[value] = value
            else:
                row.pop(value)
    coverage_dic = {}
    for i in coverage_list:
        coverage_dic[i[0]] = i[1:]

    coverage_dic = convert_index_to_locus(coverage_dic, reference_gb)

    final_consensus_w_ref = []
    final_consensus_no_ref = []
    for record in SeqIO.parse(reference_fasta, "fasta"):
        final_consensus_w_ref.append(record)
    for sample in coverage_dic:
        for file in input_files:
            if sample in file:
                record: SeqIO.SeqRecord
                for record in SeqIO.parse(file, "fasta"):
                    if get_locus_out_id(record.id) in coverage_dic[sample]:
                        final_consensus_w_ref.append(record)
                        final_consensus_no_ref.append(record)

    SeqIO.write(final_consensus_w_ref, output_w_ref, "fasta")
    SeqIO.write(final_consensus_no_ref, output_no_ref, "fasta")


if __name__ == "__main__":
    coverage_file = sys.argv[1]
    reference_gb = sys.argv[2]
    input_files = sys.argv[3].split(" ")
    reference_fasta = sys.argv[4]
    output_w_ref = sys.argv[5]
    output_no_ref = sys.argv[6]
    generate_AllConsensus(
        coverage_file,
        reference_gb,
        input_files,
        reference_fasta,
        output_w_ref,
        output_no_ref,
    )
