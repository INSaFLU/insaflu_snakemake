import sys
import csv
from Bio import SeqIO
import extract_gb_info as ggb
from yaml_io import read_yaml
import textwrap


def write_fasta(dictionary: dict[str, str], filename: str):
    with open(filename, "w") as fasta:
        for key, value in dictionary.items():
            fasta.write(f">{key}\n")
            fasta.write("\n".join(textwrap.wrap(str(value), 70)))
            fasta.write("\n")


def get_reference(file_name: str) -> SeqIO.SeqRecord:
    with open(file_name, "r") as handle:
        return next(SeqIO.parse(handle, "fasta"))


def get_ref_adjusted_positions(
    alignment: str, positions: list[tuple[str, list[list[int]]]], gene
) -> list[tuple[str, list[list[int]]]]:
    new_positions = []
    ref = get_reference(alignment)

    index_list = [idx for idx, nucleotide in enumerate(
        ref) if nucleotide != "-"]
    for gene_group in positions:
        if gene_group[0] == gene:
            for group in gene_group[1]:
                group[0] = index_list[group[0]]
                if group[1] >= len(index_list):
                    group[1] = index_list[-1]
                else:
                    group[1] = index_list[group[1]]
            new_positions.append(gene_group)
    return new_positions


def get_coverage_to_translate_matrix(filename: str) -> dict[str, list[str]]:
    with open(filename) as csvfile:
        csv_reader = csv.reader(csvfile)
        coverage_dic: dict[str, list[str]] = {
            row[0]: row[1:] for row in csv_reader}
    return coverage_dic


def get_position_in_list(reference: str, locus: str) -> int:
    list_of_locus: list[str] = ggb.get_locus(reference)
    return list_of_locus.index(locus)


def get_best_translation(record, pos: list[int]) -> str:
    options = ["", "", ""]
    options[0] = record.seq[pos[0]: pos[1]].ungap(
    ).translate(table=11, to_stop=False)
    options[1] = record.seq[pos[0] + 1: pos[1]].ungap(
    ).translate(table=11, to_stop=False)
    options[2] = record.seq[pos[0] + 2: pos[1]].ungap(
    ).translate(table=11, to_stop=False)

    idx = 0
    min_idx = 0
    min_value = options[0][:-1].count("*")
    for entry in options[1:]:
        idx += 1
        new_value = entry[:-1].count("*")
        if new_value < min_value:
            min_value = new_value
            min_idx = idx
    return options[min_idx]


def get_new_consensus(positions:  list[tuple[str, list[list[int]]]],
                      coverage_dic: dict[str, list[str]],
                      coverage_value: float,
                      gene_idx: int) -> dict[str, dict[str, str]]:
    new_consensus: dict[str, dict[str, str]] = {}
    for gene_structure in positions:
        gene_name = gene_structure[0]
        gene_positions = gene_structure[1]
        new_consensus[gene_name] = {}
        for pos in gene_positions:
            for record in SeqIO.parse(alignment, "fasta"):
                if record.id == reference_id:
                    if new_consensus[gene_name].get(record.id, False):
                        new_consensus[gene_name][record.id] += get_best_translation(
                            record, pos)
                    else:
                        new_consensus[gene_name][record.id] = get_best_translation(
                            record, pos)
                    continue

                seq = "".join([str(record.seq)[pos[0]: pos[1]].replace("-", "")
                              for pos in gene_positions])
                gene_length = sum([pos[1] - pos[0] for pos in gene_positions])
                if seq.count("N") / gene_length > 0.1:
                    continue
                identifier = record.id
                record_coverage = float(
                    coverage_dic[identifier[: identifier.index(
                        f"__{reference_id}")]][gene_idx]
                )
                if record_coverage >= coverage_value:
                    if new_consensus[gene_name].get(identifier, False):
                        new_consensus[gene_name][record.id] += get_best_translation(
                            record, pos)
                    else:
                        new_consensus[gene_name][identifier] = get_best_translation(
                            record, pos)
    return new_consensus


def write_fast_aa(
    reference: str,
    alignment: str,
    output: str,
    reference_id: str,
    gene: str,
    coverage: str,
):
    coverage_dic = get_coverage_to_translate_matrix(coverage)

    constants: dict[str, str] = read_yaml("../config/constants.yaml")

    coverage_value: int = read_yaml(constants["software_parameters"])[
        "min_coverage_consensus"
    ]

    gene_idx = get_position_in_list(reference, reference_id)

    positions: list[tuple[str, list[list[int]]]
                    ] = ggb.get_positions_gb(reference)
    positions = get_ref_adjusted_positions(alignment, positions, gene)

    new_consensus = get_new_consensus(positions,
                                      coverage_dic,
                                      coverage_value,
                                      gene_idx)

    for value in new_consensus.values():
        write_fasta(value, output)


if __name__ == "__main__":
    reference = sys.argv[1]
    alignment = sys.argv[2]
    output = sys.argv[3]
    reference_id = sys.argv[4]
    gene = sys.argv[5]
    coverage = sys.argv[6]

    write_fast_aa(reference, alignment, output, reference_id, gene, coverage)
