import argparse

from Bio import SeqIO
from Bio.Seq import MutableSeq


def compute_masking_sites(
    sequence,
    ranges=None,
    singles_positions=None,
    from_beggining=None,
    from_end=None,
):
    length = len(sequence)
    masking_sites = []
    if ranges is not None:
        for pos in ranges:
            pos[0] = int(pos[0])
            pos[1] = int(pos[1])
            if pos[1] < length and pos[1] > pos[0]:
                masking_sites.extend(iter(range(pos[0] - 1, pos[1] - 1)))
    if single_positions is not None:
        for pos in single_positions:
            pos = int(pos)
            if pos < length:
                masking_sites.append(pos - 1)
    if from_beggining is not None:
        masking_sites.extend(iter(range(int(from_beggining))))
    if from_end is not None:
        masking_sites.extend(iter(range(length - int(from_end), length)))
    return masking_sites


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("consensus", type=str, help="Consensus File Path")
    parser.add_argument("output", type=str, help="Output File Path")
    parser.add_argument(
        "-r",
        "--regions",
        type=str,
        help="Pass a string with the format x-y,n-z,...",
    )
    parser.add_argument(
        "-s",
        "--single",
        type=str,
        help="Pass a string with the format x,y,n,z,...",
    )
    parser.add_argument(
        "-b", "--beggining", type=int, help="Pass a integer with the format x"
    )
    parser.add_argument(
        "-e", "--end", type=int, help="Pass a integer with the format x"
    )

    args = parser.parse_args()

    ranges = None
    consensus = args.consensus or None
    output = args.output or None
    # Need to Error if output is missing
    single_positions = args.single.split(",") if args.single else None
    if args.regions:
        ranges = [x.split("-") for x in args.regions.split(",")]
    from_beggining = args.beggining or None
    from_end = args.end or None
    final_mask_consensus = []

    for record in SeqIO.parse(consensus, "fasta"):
        sequence = MutableSeq(record.seq)
        masking_sites = compute_masking_sites(
            sequence, ranges, single_positions, from_beggining, from_end
        )
        # Taken from insaflu
        ref_pos = 0
        ref_insertions = 0
        # for _ in range(len(sequence)):
        #     if sequence[_] == "-":
        #         ref_insertions += 1
        #         continue
        #     if ref_pos in masking_sites:
        #         sequence[ref_pos + ref_insertions] = "N"
        #     ref_pos += 1
        #     if (ref_pos + ref_insertions) >= len(record.seq):
        #         break
        for position in masking_sites:
            sequence[position] = "N"
        # End of insaflu code
        record.seq = sequence
        final_mask_consensus.append(record)

    with open(output, "w") as handle_fasta_out:
        SeqIO.write(final_mask_consensus, handle_fasta_out, "fasta")
