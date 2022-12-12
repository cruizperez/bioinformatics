#!/usr/bin/env python

import argparse
from pathlib import Path
from sys import argv
from typing import Any, Dict, Optional, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO

from bioinformatics.sequence_data.helpers import wrap_sequence_string


def get_targets_from_config(input_config: Path) -> Dict[str, Tuple[Any, Any]]:
    """
    Reads a tab-delimited configuration file with the following columns:
    sequence_id   start_position_to_mask   end_position_to_mask
    Returns a dictionare whose keys are the sequence ids and the values a tuple containing the start and end positions.

    Args:
        input_config (Path): Path to input configuration file.

    Returns:
        Dict[str, Tuple[Any, Any]]: Dictionary with sequence ids and coordinates to mask.
    """
    # Initialize empty dictionary to store sequence ids and positions
    target_sequences = {}

    # Read configuration file
    config = pd.read_csv(input_config, sep="\t", header=None)
    config.fillna(-1, inplace=True)

    # Iterate over configuration file and get targets
    for info_row in config.itertuples():
        target_sequences[info_row._1] = (info_row._2, info_row._3)

    # Return target sequence dictionary
    return target_sequences


def modify_fasta(
    input_file: Path,
    output_file: Path,
    target_sequences: Dict[str, Tuple[Any, Any]],
    masking_char: Optional[str] = "N",
    wrapping_len: Optional[int] = None,
) -> None:
    """
    Reads a FASTA file, and identifies sequences that should be masked with a given masking character (N) by default.
    Then saves the entire original file and the masked sequences into another output FASTA file.
    Lines can be wrapped to a desired length. 

    Args:
        input_file (Path): Input FASTA file.
        output_file (Path): Output FASTA file.
        target_sequences (Dict[str, Tuple[Any, Any]]): Dictionary with sequences to mask in the keys and a tuple
                                                       containing the coordinates to mask as values.
        masking_char (Optional[str], optional): Character to mask original sequence. Defaults to "N".
        wrapping_len (Optional[int], optional): Line length to which sequence should be wrapped. Defaults to None.
    """
    # Read the FASTA file
    sequences = SeqIO.parse(input_file, "fasta")

    # Open output file
    with open(output_file, "w") as output:
        # Iterate over the sequences in the file
        for record in sequences:
            # Check if the current sequence is the one we want to modify
            if record.id in target_sequences:
                # Get target positions
                start_pos, end_pos = target_sequences[record.id]
                start_index = 0 if np.isnan(start_pos) else start_pos - 1
                end_index = len(record.seq) if np.isnan(end_pos) else end_pos

                # Replace the section of the sequence with 'N's
                record.seq = (
                    record.seq[:start_index] + masking_char * (end_index - start_index) + record.seq[end_index:]
                )
            if wrapping_len:
                wrapped_sequence = wrap_sequence_string(sequence_string=str(record.seq), line_width=wrapping_len)
                # Write new masked sequence to file
                output.write(f">{record.id}\n{wrapped_sequence}\n")
            else:
                # Write original sequence to file
                output.write(f">{record.id}\n{record.seq}\n")


def parse_args() -> argparse.Namespace:  # pragma: no cover
    # Create a new ArgumentParser object
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description=(
            "This script reads a FASTA file and masks a sequence or part of it with Ns.\n"
            f"Usage: {argv[0]} -i [input fasta] -o [output fasta] -c [configuration file]\n"
        ),
    )
    required_args = parser.add_argument_group("Required arguments:")

    # Define the command-line arguments
    # The input FASTA file path
    required_args.add_argument("-i", "--input", required=True, type=Path, help="The path to the input FASTA file.")

    # The output FASTA file path
    required_args.add_argument("-o", "--output", required=True, type=Path, help="The path to the output FASTA file.")

    # The input file describing the sequences and positions to mask
    required_args.add_argument(
        "-c",
        "--config",
        required=True,
        type=Path,
        help=(
            "The path to the configuration file.\n"
            "This tab-separated file must contain three columns:\n"
            "sequence ID, start and stop positions to mask.\n"
            "If the second or third column are empty,\n"
            "they will be assumed to be the first\n"
            "and last positions in the sequence."
        ),
    )

    optional_args = parser.add_argument_group("Optional arguments:")
    # The masking character
    optional_args.add_argument("-m", "--mask", default="N", help="The character to use for masking. Default: 'N'.")

    # The length of the wrapped sequences
    optional_args.add_argument(
        "-w",
        "--wrap",
        required=False,
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
    # Parse the command-line arguments
    args = parse_args()
    input_file = args.input
    output_file = args.output
    config_file = args.config
    masking_char = args.mask
    wrapping_len = args.wrap

    # Read configuration file
    target_sequences = get_targets_from_config(config_file)

    # Modify the FASTA file
    modify_fasta(input_file, output_file, target_sequences, masking_char, wrapping_len)


# Call the main() function if the script is run directly
if __name__ == "__main__":  # pragma: no cover
    main()
