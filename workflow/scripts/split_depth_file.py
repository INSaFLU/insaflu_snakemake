"""tbd"""
import sys
from extract_gb_info import get_locus


def split_depth_file(file_path, reference_gb):
    """
    The split_depth_file function takes a file path as an argument and splits the depth file into separate files for each reference sequence.
    The function takes the last position of a string containing the full path to the depth file as its starting point, then finds where in that string
    the first instance of "/" occurs. It then uses this index to slice out everything after it and stores it in a variable called filename. The function
    then creates another variable called new_file; which is empty at first but will eventually contain all lines from our original depth file except for
    those corresponding to our reference sequence (which we know because they start with the number)

    :param file_path: Specify the path to the file that is being split
    :return: A list of lines from the depth file that correspond to the reference sequence
    :doc-author: Trelent
    """
    last_pos = len(file_path) - file_path[::-1].index("/")
    path = file_path[:last_pos]
    reference_list = get_locus(reference_gb)
    for name in reference_list:
        with open(file_path, "r", encoding="utf-8") as handle_depth:
            new_file = []
            for line in handle_depth.readlines():
                if line.split()[0] == str(name):
                    new_file.append(line)
            if new_file:
                with open(f"{path}{name}.depth", "w", encoding="utf8") as output_file:
                    output_file.writelines(new_file)


if __name__ == "__main__":
    DEPTH_FILE = sys.argv[1]
    REFERENCE_GB = sys.argv[2]
    split_depth_file(DEPTH_FILE, REFERENCE_GB)
