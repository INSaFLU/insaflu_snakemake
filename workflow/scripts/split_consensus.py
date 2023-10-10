"""import Statment"""
import sys
import csv
from Bio import SeqIO
from extract_gb_info import get_locus


def parse_depth(depth_file):
    with open(depth_file, "r", encoding="UTF8") as depth:
        reader = csv.reader(depth)
        out = []
        for i in reader:
            out.append(i[0].split("\t"))
    return out


def get_consensus(consensus):
    with open(consensus, "r", encoding="UTF8") as handle_fasta:
        return list(SeqIO.parse(handle_fasta, "fasta"))[0].seq


def split_consensus(consensus, depth, locus, output):
    depth = parse_depth(depth)
    locus = get_locus(locus)
    dic = {}
    for i in locus:
        dic[i] = ""
    for entry in depth:
        dic[entry[0]] = entry[1]
    numbers_list = []
    for i in locus:
        numbers_list.append(int(dic[i]))
    consensus_seq = get_consensus(consensus)
    # consensus_seq = "".join(['A' for x in range(0,13136)])

    last_list = []
    # count_nt = 0
    seq = " "

    new_dic = {}
    count = 0
    for loc in locus:
        new_dic[loc] = [count, count + int(dic[loc])]
        count += int(dic[loc])

    for loc in new_dic:
        small_part = consensus_seq[new_dic[loc][0] : new_dic[loc][1]]
        seq = small_part
        last_list.append(seq)

    with open(output, "w", encoding="UTF8") as out:
        for idx, seq in enumerate(last_list):
            out.writelines(f">{locus[idx]}\n")
            out.writelines(seq)
            out.writelines("\n")


if __name__ == "__main__":
    CONSENSUS = sys.argv[1]
    DEPTH = sys.argv[2]
    LOCUS = sys.argv[3]
    OUTPUT = sys.argv[4]
    split_consensus(consensus=CONSENSUS, depth=DEPTH, locus=LOCUS, output=OUTPUT)
