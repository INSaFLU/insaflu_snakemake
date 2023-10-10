"""Concatenates multiple fasta files into one"""
import sys
import csv
import re
from Bio import SeqIO
from extract_gb_info import get_locus
from yaml_io import read_yaml


def get_fasta_reference_concat(fasta_ref):
    """
    The get_fasta_reference_concat function takes a fasta file as inputand returns the concatenated
    reference sequence.
    The function is used to create a single, concatenated reference sequence from multiple contigs
    in the assembly.

    :param fasta_ref: Specify the reference genome to be used for mapping
    :return: The concatenated reference sequence in the form of a seqrecord
    """
    concat_reference = []
    count = 0
    for record in SeqIO.parse(fasta_ref, "fasta"):
        if count == 0:
            concat_reference.append(record)
            start = fasta_ref.index("reference")
            end = fasta_ref.index(".fasta")
            # why + 11 => 16/05
            concat_reference[0].id = fasta_ref[start + 11 : end]
        else:
            concat_reference[0].seq += record.seq
        count += 1
    return concat_reference[0]


def get_id_sequence_from_consensus(consensus, sample_name, value_list):
    """
    The get_id_sequence_from_consensus function takes a fasta file of consensus sequences and
    returns the sequence of the locus that is specified by the value_list. The first record in this
    fasta file will be returned with its id changed to sample_name.
    This function assumes that there are no duplicate records in consensus, which may not be true.

    :param consensus: Define the consensus sequence file
    :param sample_name: Name the sequence in the seqrecord object
    :param value_list: Specify which segments to concatenate
    :return: A seqrecord object that contains the consensus sequence for a given sample
    :doc-author: Trelent
    """

    concat_consensus = []
    count = 1

    starting_segment = value_list[0]
    for record in SeqIO.parse(consensus, "fasta"):
        if count == starting_segment:
            concat_consensus.append(record)
            concat_consensus[0].id = sample_name
        elif count in value_list:
            concat_consensus[0].seq += record.seq
        count += 1
    return concat_consensus[0]


def get_id_sequence_from_consensus_strict(consensus, sample_name, value_list):
    """
    The get_id_sequence_from_consensus_strict function takes a consensus sequence file, the sample name, and a list of integers.
    It then parses through the consensus sequence file to find all sequences that match the starting segment in the list.
    The function concatenates these sequences together into one record and assigns it an ID based on its sample name.  It returns this record.

    :param consensus: Specify the path to a fasta file containing the consensus sequence for each sample
    :param sample_name: Name the sequence
    :param value_list: Specify the segments that will be concatenated
    :return: A seqrecord object from a fasta file
    :doc-author: Trelent
    """
    concat_consensus = []
    count = 1

    starting_segment = value_list[0]
    for record in SeqIO.parse(consensus, "fasta"):
        if count == starting_segment:
            concat_consensus.append(record)
            concat_consensus[0].id = ""
            concat_consensus[0].id = sample_name
            concat_consensus[0].description = ""
        elif count in value_list:
            concat_consensus[0].seq += record.seq
        count += 1
    return concat_consensus[0]


def create_consensus_file_for_alignment(
    consensus_list: list,
    reference_gb: str,
    output: str,
    coverage_file: str,
    reference_fasta: str,
    output_only_90_plus: str,
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
    with open(coverage_file, encoding="UTF8", newline="") as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=",")
        coverage_list = list(csv_reader)
    dic_directory = read_yaml("../config/constants.yaml")

    software_parameters = read_yaml(dic_directory["software_parameters"])
    coverage_value = software_parameters["min_coverage_consensus"]
    for row in coverage_list:
        for value in range(len(row) - 1, 0, -1):
            if float(row[value]) > coverage_value:
                row[value] = value
            else:
                row.pop(value)

    coverage_dic = {}
    for i in coverage_list:
        coverage_dic[i[0]] = i[1:]

    fasta_reference_concat = get_fasta_reference_concat(reference_fasta)
    new_concat_file = [fasta_reference_concat]
    new_concat_file_strict = [fasta_reference_concat]
    if len(get_locus(reference_gb)) > 1:
        for file in consensus_list:
            sample_name = re.findall("(?<=sample_)(.*?)(?=/)", file)[0]
            if len(coverage_dic[sample_name]) > 0:
                new_sequence = get_id_sequence_from_consensus(
                    file, sample_name, coverage_dic[sample_name]
                )
                new_concat_file.append(new_sequence)
            if len(coverage_dic[sample_name]) == 8:
                new_sequence = get_id_sequence_from_consensus_strict(
                    file, sample_name, coverage_dic[sample_name]
                )
                new_concat_file_strict.append(new_sequence)
        SeqIO.write(new_concat_file, output_only_90_plus, "fasta")
        SeqIO.write(new_concat_file_strict, output, "fasta")

    else:
        for file in consensus_list:
            sample_name = re.findall("(?<=sample_)(.*?)(?=/)", file)[0]
            if len(coverage_dic[sample_name]) == 1:
                new_sequence = get_id_sequence_from_consensus_strict(
                    file, sample_name, coverage_dic[sample_name]
                )
                new_concat_file_strict.append(new_sequence)
        SeqIO.write(new_concat_file_strict, output_only_90_plus, "fasta")
        SeqIO.write(new_concat_file_strict, output, "fasta")


if __name__ == "__main__":
    CONSENSUS_LIST = sys.argv[1].split(" ")
    REFERENCE_GB = sys.argv[2]
    OUTPUT_FILE = sys.argv[3]
    COVERAGE_FILE = sys.argv[4]
    REFERENCE_FASTA = sys.argv[5]
    OUTPUT_FILE_WITH_90_PLUS_COVERAGE = sys.argv[6]

    create_consensus_file_for_alignment(
        CONSENSUS_LIST,
        REFERENCE_GB,
        OUTPUT_FILE,
        COVERAGE_FILE,
        REFERENCE_FASTA,
        OUTPUT_FILE_WITH_90_PLUS_COVERAGE,
    )
