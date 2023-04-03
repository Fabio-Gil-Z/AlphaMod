import argparse
from textwrap import wrap

parser = argparse.ArgumentParser(description="Opens a '.fasta file and transforms it to '.ali")

parser.add_argument('-d', '--fastaFilePath', metavar='', required=True, help="Path to .fasta file")

args = parser.parse_args()

if __name__ == '__main__':
    with open(args.fastaFilePath) as fasta:
        tempList = fasta.readlines()
    with open("inputSequence.ali", 'w') as inputSequence:
        inputSequence.write(">P1;inputSequence\n")
        inputSequence.write("sequence: inputSequence:::::::0.00: 0.00\n")
        sequence = wrap(tempList[1], 75)  # cuts the sequence in slices of 76 amiacids
        sequence[-1] = sequence[-1] + "*"  # adds the asterisk to the last slice of the proteinSequence to indicate the end of the sequence for Modeller
        for slice in sequence:
            inputSequence.write(slice + "\n")




