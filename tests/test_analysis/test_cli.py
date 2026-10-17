"""CLI smoke tests for analysis commands."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

from click.testing import CliRunner

from pyjaspar.cli.main import cli


def _write_fasta(path: Path, sequences: dict[str, str]) -> None:
    """Write a simple FASTA file."""
    with open(path, "w") as f:
        for name, seq in sequences.items():
            f.write(f">{name}\n{seq}\n")


def test_analysis_commands_in_help():
    result = CliRunner().invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "scan" in result.output
    assert "similarity" in result.output
    assert "align" in result.output
    assert "enrichment" in result.output
    assert "search" in result.output
    assert "infer" in result.output
    assert "logo" in result.output


# --- scan tests ---


def test_scan_literal():
    result = CliRunner().invoke(
        cli,
        [
            "scan",
            "ACGTACGTACGTACGT" * 5,
            "--motif-id",
            "MA0001.1",
            "--threshold",
            "0.3",
        ],
    )
    assert result.exit_code == 0


def test_scan_json():
    result = CliRunner().invoke(
        cli,
        [
            "scan",
            "ACGTACGTACGTACGT" * 5,
            "--motif-id",
            "MA0001.1",
            "--threshold",
            "0.3",
            "--format",
            "json",
        ],
    )
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert "results" in data
    assert "count" in data
    assert data["motif_id"] == "MA0001.1"


def test_scan_fasta_file():
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = Path(tmpdir) / "test.fasta"
        _write_fasta(
            fasta_path,
            {
                "seq1": "ACGTACGTACGTACGT" * 5,
                "seq2": "TGCATGCATGCATGCA" * 5,
            },
        )
        result = CliRunner().invoke(
            cli,
            [
                "scan",
                str(fasta_path),
                "--motif-id",
                "MA0001.1",
                "--threshold",
                "0.3",
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "results" in data


def test_scan_invalid_motif():
    result = CliRunner().invoke(
        cli,
        ["scan", "ACGTACGT", "--motif-id", "MA9999.1"],
    )
    assert result.exit_code == 1


def test_scan_no_reverse():
    result = CliRunner().invoke(
        cli,
        [
            "scan",
            "ACGTACGTACGTACGT" * 5,
            "--motif-id",
            "MA0001.1",
            "--threshold",
            "0.3",
            "--no-reverse",
            "--format",
            "json",
        ],
    )
    assert result.exit_code == 0


# --- similarity tests ---


def test_similarity_pearson():
    result = CliRunner().invoke(cli, ["similarity", "MA0001.1", "MA0002.1"])
    assert result.exit_code == 0


def test_similarity_euclidean():
    result = CliRunner().invoke(
        cli, ["similarity", "MA0001.1", "MA0002.1", "--metric", "euclidean"]
    )
    assert result.exit_code == 0


def test_similarity_kl():
    result = CliRunner().invoke(cli, ["similarity", "MA0001.1", "MA0002.1", "--metric", "kl"])
    assert result.exit_code == 0


def test_similarity_best():
    result = CliRunner().invoke(cli, ["similarity", "MA0001.1", "MA0002.1", "--best"])
    assert result.exit_code == 0


def test_similarity_json():
    result = CliRunner().invoke(cli, ["similarity", "MA0001.1", "MA0002.1", "--format", "json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert "score" in data
    assert data["motif1"] == "MA0001.1"
    assert data["motif2"] == "MA0002.1"


def test_similarity_best_json():
    result = CliRunner().invoke(
        cli, ["similarity", "MA0001.1", "MA0002.1", "--best", "--format", "json"]
    )
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert "best_score" in data
    assert "best_offset" in data
    assert "is_reverse_complement" in data


def test_similarity_invalid_motif():
    result = CliRunner().invoke(cli, ["similarity", "MA9999.1", "MA0002.1"])
    assert result.exit_code == 1


def test_similarity_with_offset():
    result = CliRunner().invoke(
        cli, ["similarity", "MA0001.1", "MA0002.1", "--offset", "2", "--format", "json"]
    )
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert data["offset"] == 2


# --- align tests ---


def test_align_basic():
    result = CliRunner().invoke(cli, ["align", "MA0001.1", "MA0002.1"])
    assert result.exit_code == 0
    assert "target" in result.output
    assert "query" in result.output


def test_align_no_reverse_complement():
    result = CliRunner().invoke(cli, ["align", "MA0001.1", "MA0002.1", "--no-reverse-complement"])
    assert result.exit_code == 0


def test_align_json():
    result = CliRunner().invoke(cli, ["align", "MA0001.1", "MA0002.1", "--format", "json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert data["motif1"] == "MA0001.1"
    assert data["motif2"] == "MA0002.1"
    assert "score" in data
    assert "offset" in data
    assert "is_reverse_complement" in data


def test_align_invalid_motif():
    result = CliRunner().invoke(cli, ["align", "MA9999.1", "MA0002.1"])
    assert result.exit_code == 1


# --- search tests ---


def test_search_tsv_lists_the_top_hits():
    result = CliRunner().invoke(
        cli, ["search", "MA0139.1", "--tax-group", "Vertebrates", "--top", "3"]
    )
    assert result.exit_code == 0
    lines = result.output.strip().splitlines()
    assert lines[0].split("\t") == [
        "rank",
        "matrix_id",
        "name",
        "score",
        "offset",
        "overlap",
        "gaps",
        "is_reverse_complement",
    ]
    assert len(lines) == 4
    assert [line.split("\t")[2] for line in lines[1:]] == ["CTCF", "CTCF", "CTCF"]


def test_search_json_pearson():
    result = CliRunner().invoke(
        cli,
        [
            "search",
            "MA0139.1",
            "--tax-group",
            "Vertebrates",
            "--method",
            "pearson",
            "--top",
            "2",
            "--format",
            "json",
        ],
    )
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert data["method"] == "pearson"
    assert data["count"] == 2
    assert data["results"][0]["percent_score"] is None


def test_search_invalid_motif():
    result = CliRunner().invoke(cli, ["search", "MA9999.1"])
    assert result.exit_code == 1


# --- infer tests ---


def test_infer_tsv(monkeypatch):
    from pyjaspar.analysis import InferenceHit

    def fake_infer(sequence, release):
        return [InferenceHit("MA0162.2", "EGR1", 1.9e-63, 0.942, "logo.svg", release)]

    monkeypatch.setattr("pyjaspar.analysis.infer_profiles", fake_infer)
    result = CliRunner().invoke(cli, ["infer", "MSNEQ"])
    assert result.exit_code == 0
    assert "input\tJASPAR2024\tMA0162.2\tEGR1" in result.output


def test_infer_unavailable_release():
    result = CliRunner().invoke(cli, ["infer", "MSNEQ", "--release", "2026"])
    assert result.exit_code == 1
    assert "not available" in result.output


# --- logo tests ---


def test_logo_writes_an_image(tmp_path):
    out = tmp_path / "ctcf.png"
    result = CliRunner().invoke(cli, ["logo", "MA0139.1", "-o", str(out)])
    assert result.exit_code == 0
    assert out.stat().st_size > 0


def test_logo_invalid_motif(tmp_path):
    result = CliRunner().invoke(cli, ["logo", "MA9999.1", "-o", str(tmp_path / "x.png")])
    assert result.exit_code == 1


# --- enrichment tests ---


def test_enrichment_basic():
    with tempfile.TemporaryDirectory() as tmpdir:
        fg_path = Path(tmpdir) / "fg.fasta"
        bg_path = Path(tmpdir) / "bg.fasta"
        _write_fasta(fg_path, {f"fg{i}": "ACGTACGTACGTACGT" * 5 for i in range(5)})
        _write_fasta(bg_path, {f"bg{i}": "TGCATGCATGCATGCA" * 5 for i in range(5)})

        result = CliRunner().invoke(
            cli,
            [
                "enrichment",
                "--foreground",
                str(fg_path),
                "--background",
                str(bg_path),
                "--motif-ids",
                "MA0001.1",
                "--threshold",
                "0.3",
            ],
        )
        assert result.exit_code == 0


def test_enrichment_json():
    with tempfile.TemporaryDirectory() as tmpdir:
        fg_path = Path(tmpdir) / "fg.fasta"
        bg_path = Path(tmpdir) / "bg.fasta"
        _write_fasta(fg_path, {f"fg{i}": "ACGTACGTACGTACGT" * 5 for i in range(5)})
        _write_fasta(bg_path, {f"bg{i}": "TGCATGCATGCATGCA" * 5 for i in range(5)})

        result = CliRunner().invoke(
            cli,
            [
                "enrichment",
                "--foreground",
                str(fg_path),
                "--background",
                str(bg_path),
                "--motif-ids",
                "MA0001.1",
                "--threshold",
                "0.3",
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert "results" in data
        assert "count" in data


def test_enrichment_multiple_motifs():
    with tempfile.TemporaryDirectory() as tmpdir:
        fg_path = Path(tmpdir) / "fg.fasta"
        bg_path = Path(tmpdir) / "bg.fasta"
        _write_fasta(fg_path, {f"fg{i}": "ACGTACGTACGTACGT" * 5 for i in range(3)})
        _write_fasta(bg_path, {f"bg{i}": "TGCATGCATGCATGCA" * 5 for i in range(3)})

        result = CliRunner().invoke(
            cli,
            [
                "enrichment",
                "--foreground",
                str(fg_path),
                "--background",
                str(bg_path),
                "--motif-ids",
                "MA0001.1,MA0002.1",
                "--threshold",
                "0.3",
                "--format",
                "json",
            ],
        )
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert data["count"] == 2


def test_enrichment_invalid_motif_warns():
    with tempfile.TemporaryDirectory() as tmpdir:
        fg_path = Path(tmpdir) / "fg.fasta"
        bg_path = Path(tmpdir) / "bg.fasta"
        _write_fasta(fg_path, {"fg1": "ACGTACGTACGTACGT" * 5})
        _write_fasta(bg_path, {"bg1": "TGCATGCATGCATGCA" * 5})

        result = CliRunner().invoke(
            cli,
            [
                "enrichment",
                "--foreground",
                str(fg_path),
                "--background",
                str(bg_path),
                "--motif-ids",
                "MA9999.1",
            ],
        )
        assert result.exit_code == 1


# --- FASTA detection ---


def test_is_fasta_file_ignores_long_literal_sequences():
    from pyjaspar.cli.analysis import _is_fasta_file

    assert _is_fasta_file("ACGT" * 200) is False
