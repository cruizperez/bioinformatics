from pathlib import Path

import pytest

from bioinformatics.sequence_data.fasta_chop_sequence import chop_sequences


@pytest.fixture
def mock_fasta_file(tmp_path: Path) -> Path:
    """Create a FASTA file for testing"""
    input_fasta = tmp_path / "temp.fasta"
    with open(input_fasta, "w") as input_fasta_fh:
        input_fasta_fh.write(">seq1\nATCGATCGATCGATACGACTAGCA\n")
        input_fasta_fh.write(">seq2\nCGAATTAACGACAATCAGCTACGG\n")
        input_fasta_fh.write(">seq3\nTACGATACTGACCATACGATCAGC\n")
    return input_fasta


@pytest.fixture
def mock_config(tmp_path: Path) -> Path:
    """Create a configuration file for testing"""
    config_file = tmp_path / "temp.config"
    with open(config_file, "w") as config_file_fh:
        config_file_fh.write("seq1\n")
        config_file_fh.write("seq2\n")
    return config_file


def test_chop_sequences(tmp_path: Path, mock_fasta_file: Path, mock_config: Path) -> None:
    output_file = tmp_path / "output.fasta"

    chop_sequences(mock_fasta_file, output_file, mock_config, seq_len=10)

    expected_output = [
        ">seq1_1_10",
        "ATCGATCGAT",
        ">seq1_11_20",
        "CGATACGACT",
        ">seq1_21_30",
        "AGCA",
        ">seq2_1_10",
        "CGAATTAACG",
        ">seq2_11_20",
        "ACAATCAGCT",
        ">seq2_21_30",
        "ACGG",
    ]

    assert output_file.read_text().splitlines() == expected_output


def test_chop_sequences_wrap(tmp_path: Path, mock_fasta_file: Path, mock_config: Path) -> None:
    output_file = tmp_path / "output.fasta"

    chop_sequences(mock_fasta_file, output_file, mock_config, seq_len=10, wrap_len=5)

    expected_output = [
        ">seq1_1_10",
        "ATCGA",
        "TCGAT",
        ">seq1_11_20",
        "CGATA",
        "CGACT",
        ">seq1_21_30",
        "AGCA",
        ">seq2_1_10",
        "CGAAT",
        "TAACG",
        ">seq2_11_20",
        "ACAAT",
        "CAGCT",
        ">seq2_21_30",
        "ACGG",
    ]

    assert output_file.read_text().splitlines() == expected_output
