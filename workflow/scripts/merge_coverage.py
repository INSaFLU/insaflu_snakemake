"""Create coverage report"""
import re
import sys
import csv
import yaml

from Bio import SeqIO


def get_locus(genbank_file):
    """
    The get_locus function takes a genbank file and returns the locus number of that record.
    If there is only one record in the genbank file, it will return the possible_name.
    If there are multiple records, it will return a list of all locus numbers.

    :param genbank_file: Open the genbank file, and then parse it using seqio
    :param possible_name: Determine if the file is a single genbank file or not
    :return: A list of locus numbers
    :doc-author: Trelent
    """
    locus = []
    with open(genbank_file, encoding="utf-8") as handle_gb:
        for record in SeqIO.parse(handle_gb, "genbank"):
            locus.append(record.name)
    return locus


def get_locus_len(genbank_file):
    """
    The get_locus function takes a genbank file and returns the locus number of that record.
    If there is only one record in the genbank file, it will return the possible_name.
    If there are multiple records, it will return a list of all locus numbers.

    :param genbank_file: Open the genbank file, and then parse it using seqio
    :param possible_name: Determine if the file is a single genbank file or not
    :return: A list of locus numbers
    :doc-author: Trelents
    """
    locus = []
    with open(genbank_file, encoding="utf-8") as handle_gb:
        for record in SeqIO.parse(handle_gb, "genbank"):
            locus.append(len(record.seq))
    return locus


def read_yaml(yaml_file_path: str) -> dict:
    """
    The read_yaml function reads a yaml file and returns the contents as a dictionary.
    The read_yaml function accepts one argument, which is the path to the yaml file.

    :param yaml_file_path: Specify the path to the yaml file that is going to be read
    :return: A dictionary
    :doc-author: Trelent
    """
    with open(yaml_file_path, encoding="utf-8") as yaml_file:
        return yaml.load(yaml_file, Loader=yaml.FullLoader)


def merge_coverage(file_list, output_file):
    """
    The merge_coverage function takes a list of coverage files and merges them into one file.
    The first argument is the list of files, the second argument is the name for the output file.

    :param file_list: Store the list of files that will be merged
    :return: A list of lists with the coverage information for each sample
    :doc-author: Trelent
    """

    with open(file_list[0], "r", encoding="utf8") as file_handler:
        length = file_handler.readlines()[3].split()[1]

    species = get_locus(REFERENCE_GB)
    length = get_locus_len(REFERENCE_GB)
    super_header = ["Name", *species, "", ""]
    super_header_1 = ["Length", *length, "", ""]

    software_parameters = read_yaml(dic_directory["software_parameters"])

    header = [
        "SAMPLES",
        "Mean depth of coverage",
        r"% of size covered by at least 1-fold",
        r"% of size covered by at least X-fold",
    ]
    sub_header = ["", *species, "", *species, "X-Fold", *species]
    final_output = [super_header, super_header_1, header, sub_header]

    for file in file_list:
        with open(file, "r", encoding="utf8") as f:
            sample_name = re.findall("(?<=/)(.*?)(?=_coverage.csv)", file)[0].split(
                "/"
            )[0]
            info = f.readlines()[6].split()
            sample_t = sample_type[sample_name]
            if sample_t in ("snippy", "iVar"):
                sample_cov = software_parameters["mincov"]
            elif sample_t == "medaka":
                sample_cov = software_parameters["mincov_medaka"]
            else:
                raise Exception(
                    "There is no implementation for other tools than snippy, iVar and medaka."
                )
            info[0] = sample_name
            info.insert(1 + len(species), "")
            info.insert(2 + len(species) * 2, sample_cov)
            final_output.append(info)

    with open(output_file, mode="w", encoding="utf8") as f:
        f_writer = csv.writer(
            f, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL
        )
        f_writer.writerows(final_output)


def get_coverage(filename: str, n_locus: int):
    """
    The get_coverage function takes a filename as an argument and returns the coverage of each
    locus in that file.
    The function also takes n_locus as an argument, which is the number of loci in the file.

    :param filename: Specify the file to be read
    :param n_locus: Specify the number of loci in the file
    :return: A list of the coverage for each locus in a given sample
    :doc-author: Trelent
    """
    with open(filename, newline="", encoding="UTF8") as csvfile:
        coverage_list = list(csv.reader(csvfile, delimiter="\t"))
        sample_name = re.findall("(?<=/)(.*?)(?=_coverage.csv)", filename)[0].split(
            "/"
        )[0]
        coverage_list = coverage_list[-1][-n_locus:]
        coverage_list.insert(0, sample_name)
        csvfile.close()
    return coverage_list


def create_script_readable_coverage_file(coverage_files, output, n_locus):
    """
    The create_script_readable_coverage_file function takes a list of coverage files and outputs
    a script readable coverage file. The output is the name of the file you want to create.
    The last argument is n_locus which specifies how many loci are in each coverage file.

    :param coverage_files: Specify the files that contain the coverage data
    :param output: Specify the name of the output file
    :param n_locus: Indicate the number of locus in the genome
    :return: A list of lists, where each sublist contains the coverage each locus
    :doc-author: Trelent
    """
    coverage_list = coverage_files

    n_locus = int(n_locus)
    coverage_translate = []
    for filename in coverage_list:
        coverage_translate.append(get_coverage(filename, n_locus))

    with open(output, mode="w", encoding="UTF8") as output_file:
        f_writer = csv.writer(output_file, delimiter=",")
        f_writer.writerows(coverage_translate)
        output_file.close()


if __name__ == "__main__":
    COVERAGE_FILES = sys.argv[1].split()
    OUTPUT_FILE_1 = sys.argv[2]
    OUTPUT_FILE_2 = sys.argv[3]
    REFERENCE_GB = sys.argv[4]

    dic_directory = read_yaml("../config/constants.yaml")
    sample_type = read_yaml(dic_directory["config_file"])["sample_type"]

    NUMBER_OF_LOCUS = len(get_locus(REFERENCE_GB))
    merge_coverage(COVERAGE_FILES, OUTPUT_FILE_1)
    create_script_readable_coverage_file(
        COVERAGE_FILES, OUTPUT_FILE_2, NUMBER_OF_LOCUS)
