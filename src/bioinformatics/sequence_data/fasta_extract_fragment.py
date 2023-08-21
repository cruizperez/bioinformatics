#!/usr/bin/env python

import argparse
from pathlib import Path
from sys import argv
from typing import Optional

import pandas as pd
from Bio import SeqIO

from bioinformatics.sequence_data.helpers import wrap_sequence_string


def get_targets_from_config(input_config: Path) -> pd.DataFrame:
    """
    Reads a tab-delimited configuration file with the following columns:
    sequence_id   start_position_to_extract   end_position_to_extract
    Converts to a pandas dataframe, fills NAs with -1 and returns it.

    Args:
        input_config (Path): Path to input configuration file.

    Returns:
        Dict[str, Tuple[Any, Any]]: Dictionary with sequence ids and coordinates to mask.
    """
    # Read configuration file
    config = pd.read_csv(input_config, sep="\t", header=None)
    config.fillna(-1, inplace=True)
    config.columns = pd.core.indexes.base.Index(data=["seq_id", "start", "end"])
    config["start"] = config["start"].astype(int)
    config["end"] = config["end"].astype(int)

    # Return dataframe
    return config


def extract_from_fasta(
    input_file: Path,
    output_file: Path,
    target_sequences: pd.DataFrame,
    wrapping_len: Optional[int] = None,
) -> None:
    """
    Reads a FASTA file and extracts fragments 

    Args:
        input_file (Path): Input FASTA file.
        output_file (Path): Output FASTA file.
        target_sequences (pd.DataFrame): Dataframe with name of sequences to extract from and the start and end
            positions to extract.
        wrapping_len (Optional[int], optional): Line length to which sequence should be wrapped. Defaults to None.
    """
    # Read the FASTA file
    sequences = SeqIO.parse(input_file, "fasta")

    # Get the target sequence IDs
    seq_ids = set(target_sequences["seq_id"])

    # Open output file
    with open(output_file, "w") as output:
        # Iterate over the sequences in the file
        for record in sequences:
            # Check if the current sequence is the one we want to modify
            if record.id in seq_ids:
                tmp_data = target_sequences[target_sequences["seq_id"] == record.id]
                for section in tmp_data.itertuples():
                    # Get target positions
                    start_idx = 0 if section.start == -1 else section.start - 1
                    end_idx = len(record.seq) if section.end == -1 else section.end
                    # Get the new sequence id to extract and the sequence
                    new_id = f"{record.id}_{start_idx + 1}_{end_idx}"
                    sub_sequence = record.seq[start_idx:end_idx]
            if wrapping_len:
                sub_sequence = wrap_sequence_string(sequence_string=str(sub_sequence), line_width=wrapping_len)
            # Write new sequences to file
            output.write(f">{new_id}\n{sub_sequence}\n")


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

    # The input file describing the sequences and positions to mask
    required_args.add_argument(
        "-c",
        "--config",
        required=True,
        type=Path,
        help=(
            "The path to the configuration file.\n"
            "This tab-separated file must contain three columns:\n"
            "sequence ID, start and stop positions to extract.\n"
            "If the second or third column are empty,\n"
            "they will be assumed to be the first\n"
            "and last positions in the sequence."
        ),
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
    # Parse the command-line arguments
    args = parse_args()
    input_file = args.input
    output_file = args.output
    config_file = args.config
    wrapping_len = args.wrap

    # Read configuration file
    target_sequences = get_targets_from_config(config_file)

    # Extract from input FASTA file
    extract_from_fasta(input_file, output_file, target_sequences, wrapping_len)


# Call the main() function if the script is run directly
if __name__ == "__main__":  # pragma: no cover
    main()
