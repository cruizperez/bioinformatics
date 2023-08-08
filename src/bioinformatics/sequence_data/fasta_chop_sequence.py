#!/usr/bin/env python

"""
Chops sequences in a FASTA file into smaller non-overlapping sequences of a given length.
"""

import argparse
from pathlib import Path
from sys import argv
from typing import List, Optional
from Bio import SeqIO
from bioinformatics.sequence_data.helpers import wrap_sequence_string


def chop_sequences(
    input_file: Path, output_file: Path, config_file: Path, seq_len: int, wrap_len: Optional[int]
) -> None:
    """
    Reads a FASTA file and chops the sequences into smaller non-overlapping sequences of a given
    length.

    Args:
        input_file (Path): Input FASTA file.
        output_file (Path): Output FASTA file.
        config_file (Path): File with the list of sequences to chop.
        seq_len (int): Length of the chopped sequences.
        wrap_len (Optional[int], optional): Line length to which sequence should be wrapped. Defaults to None.
    """
    # Read the FASTA file
    sequences = SeqIO.parse(input_file, "fasta")

    # Read target sequences
    target_sequences: List[str] = config_file.read_text().splitlines()

    # Open output file
    with open(output_file, "w", encoding="utf-8") as output:
        # Iterate over the sequences in the file
        for record in sequences:
            # Check if the current sequence is in the list of target sequences
            if record.id in target_sequences:
                for position in range(0, len(record.seq), seq_len):
                    new_id = f"{record.id}_{position + 1}_{position + seq_len}"
                    new_seq = record.seq[position : position + seq_len]
                    if wrap_len:
                        wrapped_sequence = wrap_sequence_string(sequence_string=str(new_seq), line_width=wrap_len)
                        # Write new masked sequence to file
                        output.write(f">{new_id}\n{wrapped_sequence}\n")
                    else:
                        # Write original sequence to file
                        output.write(f">{new_id}\n{new_seq}\n")


def parse_args() -> argparse.Namespace:  # pragma: no cover
    """
    Parses the command-line arguments.

    Returns:
        argparse.Namespace: The parsed arguments.
    """
    # Create a new ArgumentParser object
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "This script reads a FASTA file and chops target sequences into non-overlapping "
            "substrings of length N.\n"
            f"Usage: {argv[0]} -i [input fasta] -o [output fasta] -c [configuration file]\n"
        ),
    )
    required_args = parser.add_argument_group("Required arguments:")

    # Define the command-line arguments
    # The input FASTA file path
    required_args.add_argument(
        "-i",
        "--input",
        required=True,
        type=Path,
        help="The path to the input FASTA file.",
    )

    # The output FASTA file path
    required_args.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="The path to the output FASTA file.",
    )

    # The input file describing the sequences to chop

    required_args.add_argument(
        "-c",
        "--config",
        required=True,
        type=Path,
        help="The path to the configuration file. Must contain a single sequence ID per line."
    )

    # The length of the chopped sequences
    required_args.add_argument(
        "-l",
        "--length",
        required=True,
        type=int,
        help="The length of the chopped sequences.",
    )

    optional_args = parser.add_argument_group("Optional arguments:")

    # The length of the wrapped sequences
    optional_args.add_argument(
        "-w",
        "--wrap",
        required=False,
        type=int,
        help=(
            "The length to which sequences should be wrapped.\n"
            "By default no wrapping is applied, i.e. the whole sequence\n"
            "is stored in a single line."
        ),
    )

    # If no arguments are provided
    if len(argv) == 1:
        parser.print_help()
        exit(0)

    # Parse the command-line arguments
    return parser.parse_args()


def main() -> None:  # pragma: no cover
    """
    The main function.
    """
    # Parse the command-line arguments
    args = parse_args()
    input_file = args.input
    output_file = args.output
    config = args.config
    length = args.length
    wrap_len = args.wrap

    # Read configuration file
    chop_sequences(input_file, output_file, config, length, wrap_len)


# Call the main() function if the script is run directly
if __name__ == "__main__":  # pragma: no cover
    main()
