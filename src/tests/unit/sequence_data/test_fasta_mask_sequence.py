# Import the script and the required modules
from pathlib import Path

import pytest

from bioinformatics.sequence_data.fasta_mask_sequence import (
    get_targets_from_config,
    modify_fasta,
)


# Create a FASTA file for testing
@pytest.fixture
def input_fasta(tmp_path: Path):
    input_fasta = tmp_path / "temp.fasta"
    input_fasta_fh = open(input_fasta, "w")
    input_fasta_fh.write(">seq1\nATCGATCGATCG\n")
    input_fasta_fh.write(">seq2\nATCGATCGATCG\n")
    input_fasta_fh.write(">seq3\nATCGATCGATCG\n")
    input_fasta_fh.close()
    return input_fasta


# Create a configuration file for testing
@pytest.fixture
def config_file(tmp_path):
    config_file = tmp_path / "temp.config"
    config_file_fh = open(config_file, "w")
    config_file_fh.write("seq1\t1\t7\n")
    config_file_fh.write("seq2\t\t\n")
    config_file_fh.write("seq3\t5\t\n")
    config_file_fh.close()
    return config_file


# Test the get_targets_from_config() function
def test_get_targets_from_config(config_file):
    targets = get_targets_from_config(config_file)
    assert targets == {"seq1": (1, 7), "seq2": (-1, -1), "seq3": (5, -1)}


# Test the modify_fasta() function
def test_modify_fasta(tmp_path: Path, input_fasta: Path):
    # Define output file
    outfile = tmp_path / "output.fasta"

    # Define target sequences and positions
    target_sequences = {
        "seq1": (1, 9),
        "seq2": (5, -1),
        "seq3": (-1, -1),
    }

    expected_output = ">seq1\nNNNNNNNNNTCG\n" ">seq2\nATCGNNNNNNNN\n" ">seq3\nNNNNNNNNNNNN\n"
    # Execute function
    modify_fasta(input_fasta, outfile, target_sequences)
    observed_output = outfile.read_text()
    assert observed_output == expected_output


# Test the modify_fasta() function with alternative masking character
def test_modify_fasta_alternative_character(tmp_path: Path, input_fasta: Path):
    # Define output file
    outfile = tmp_path / "output.fasta"

    # Define target sequences and positions
    target_sequences = {
        "seq1": (1, 9),
        "seq2": (5, -1),
        "seq3": (-1, -1),
    }

    expected_output = ">seq1\nXXXXXXXXXTCG\n" ">seq2\nATCGXXXXXXXX\n" ">seq3\nXXXXXXXXXXXX\n"
    # Execute function
    modify_fasta(input_fasta, outfile, target_sequences, masking_char="X")
    observed_output = outfile.read_text()
    assert observed_output == expected_output


# Test the modify_fasta() function with alternative line wrapping
def test_modify_fasta_alternative_wrapping(tmp_path: Path, input_fasta: Path):
    # Define output file
    outfile = tmp_path / "output.fasta"

    # Define target sequences and positions
    target_sequences = {
        "seq1": (1, 9),
        "seq2": (5, -1),
        "seq3": (-1, -1),
    }

    expected_output = ">seq1\nXXX\nXXX\nXXX\nTCG\n" ">seq2\nATC\nGXX\nXXX\nXXX\n" ">seq3\nXXX\nXXX\nXXX\nXXX\n"
    # Execute function
    modify_fasta(input_fasta, outfile, target_sequences, masking_char="X", wrapping_len=3)
    observed_output = outfile.read_text()
    assert observed_output == expected_output
