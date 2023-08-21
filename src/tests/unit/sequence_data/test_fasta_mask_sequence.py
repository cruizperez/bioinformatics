# Import the script and the required modules
from pathlib import Path

import pandas as pd
import pytest

from bioinformatics.sequence_data.fasta_mask_sequence import (
    get_targets_from_config,
    modify_fasta,
)


# Create a FASTA file for testing
@pytest.fixture
def input_fasta(tmp_path: Path) -> Path:
    """Create a FASTA file for testing"""
    input_fasta = tmp_path / "temp.fasta"
    with open(input_fasta, "w") as input_fasta_fh:
        input_fasta_fh.write(">seq1\nATCGATCGATCG\n")
        input_fasta_fh.write(">seq2\nATCGATCGATCG\n")
        input_fasta_fh.write(">seq3\nATCGATCGATCG\n")
    return input_fasta


# Create a configuration file for testing
@pytest.fixture
def config_file(tmp_path: Path) -> Path:
    """Create a configuration file for testing"""
    config_file = tmp_path / "temp.config"
    with open(config_file, "w") as config_file_fh:
        config_file_fh.write("seq1\t1\t7\n")
        config_file_fh.write("seq2\t\t\n")
        config_file_fh.write("seq3\t5\t\n")
    return config_file


# Test the get_targets_from_config() function
def test_get_targets_from_config(config_file: Path) -> None:
    """Test the get_targets_from_config() function"""
    # Arrange
    expected = pd.DataFrame(
        [["seq1", 1, 7], ["seq2", -1, -1], ["seq3", 5, -1]],
        columns=["seq_id", "start", "end"],
    )

    # Act
    targets = get_targets_from_config(config_file)

    # Assert
    pd.testing.assert_frame_equal(expected, targets)


# Test the modify_fasta() function
def test_modify_fasta(tmp_path: Path, input_fasta: Path) -> None:
    """Test the modify_fasta() function"""
    # Arrange
    outfile = tmp_path / "output.fasta"
    target_sequences = pd.DataFrame(
        [["seq1", 1, 9], ["seq2", 5, -1], ["seq3", -1, -1]],
        columns=["seq_id", "start", "end"],
    )
    expected_output = ">seq1\nNNNNNNNNNTCG\n" ">seq2\nATCGNNNNNNNN\n" ">seq3\nNNNNNNNNNNNN\n"

    # Act
    modify_fasta(input_fasta, outfile, target_sequences)
    observed_output = outfile.read_text()

    # Assert
    assert observed_output == expected_output


# Test the modify_fasta() function with alternative masking character
def test_modify_fasta_alternative_character(tmp_path: Path, input_fasta: Path) -> None:
    """Test the modify_fasta() function with alternative masking character"""
    # Define output file
    outfile = tmp_path / "output.fasta"

    # Define target sequences and positions
    target_sequences = pd.DataFrame(
        [["seq1", 1, 9], ["seq2", 5, -1], ["seq3", -1, -1]],
        columns=["seq_id", "start", "end"],
    )

    expected_output = ">seq1\nXXXXXXXXXTCG\n" ">seq2\nATCGXXXXXXXX\n" ">seq3\nXXXXXXXXXXXX\n"
    # Execute function
    modify_fasta(input_fasta, outfile, target_sequences, masking_char="X")
    observed_output = outfile.read_text()
    assert observed_output == expected_output


# Test the modify_fasta() function with alternative line wrapping
def test_modify_fasta_alternative_wrapping(tmp_path: Path, input_fasta: Path) -> None:
    """Test the modify_fasta() function with alternative line wrapping"""
    # Define output file
    outfile = tmp_path / "output.fasta"

    # Define target sequences and positions
    target_sequences = pd.DataFrame(
        [["seq1", 1, 9], ["seq2", 5, -1], ["seq3", -1, -1]],
        columns=["seq_id", "start", "end"],
    )

    expected_output = ">seq1\nXXX\nXXX\nXXX\nTCG\n" ">seq2\nATC\nGXX\nXXX\nXXX\n" ">seq3\nXXX\nXXX\nXXX\nXXX\n"
    # Execute function
    modify_fasta(input_fasta, outfile, target_sequences, masking_char="X", wrapping_len=3)
    observed_output = outfile.read_text()
    assert observed_output == expected_output


def test_modify_fasta_mulitple_entries(tmp_path: Path, input_fasta: Path) -> None:
    """Test the modify_fasta() function with multiple entries"""
    # Arrange
    outfile = tmp_path / "output.fasta"
    target_sequences = pd.DataFrame(
        [
            ["seq1", 1, 4],
            ["seq1", 7, 9],
            ["seq2", 5, -1],
            ["seq2", 7, 9],
            ["seq3", 1, 3],
            ["seq3", 3, -1],
        ],
        columns=["seq_id", "start", "end"],
    )
    expected_output = ">seq1\nNNNNATNNNTCG\n" ">seq2\nATCGNNNNNNNN\n" ">seq3\nNNNNNNNNNNNN\n"

    # Act
    modify_fasta(input_fasta, outfile, target_sequences)
    observed_output = outfile.read_text()

    # Assert
    assert observed_output == expected_output
