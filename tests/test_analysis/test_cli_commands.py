"""Tests of the align and similarity command objects, invoked directly."""

from __future__ import annotations

import json

import pytest
from click.testing import CliRunner

from pyjaspar import JasparDB
from pyjaspar.analysis import align_motifs, format_alignment
from pyjaspar.cli.analysis import align, similarity


@pytest.fixture(scope="module")
def motifs():
    jdb = JasparDB()
    return jdb.fetch_motif_by_id("MA0001.1"), jdb.fetch_motif_by_id("MA0002.1")


def test_similarity_best_json_is_the_pearson_alignment(motifs):
    result = CliRunner().invoke(similarity, ["MA0001.1", "MA0002.1", "--best", "--format", "json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert set(data) == {
        "motif1",
        "motif2",
        "metric",
        "best_score",
        "best_offset",
        "is_reverse_complement",
    }
    expected = align_motifs(*motifs, method="pearson")
    assert data["metric"] == "pearson"
    assert data["best_score"] == round(expected.score, 6)
    assert data["best_offset"] == expected.offset
    assert data["is_reverse_complement"] == expected.is_reverse_complement


def test_similarity_best_tsv_header():
    result = CliRunner().invoke(similarity, ["MA0001.1", "MA0002.1", "--best"])
    assert result.exit_code == 0
    header = result.output.splitlines()[0].split("\t")
    assert header == [
        "motif1",
        "motif2",
        "metric",
        "best_score",
        "best_offset",
        "is_reverse_complement",
    ]


def test_similarity_at_a_fixed_offset_is_unchanged():
    result = CliRunner().invoke(similarity, ["MA0001.1", "MA0002.1", "--offset", "2"])
    assert result.exit_code == 0
    assert result.output.splitlines()[1].split("\t") == [
        "MA0001.1",
        "MA0002.1",
        "pearson",
        "2",
        "0.471837",
    ]


@pytest.mark.parametrize("method", ["matrix_align", "pearson"])
def test_align_text_shows_the_alignment(motifs, method):
    result = CliRunner().invoke(align, ["MA0001.1", "MA0002.1", "--method", method])
    assert result.exit_code == 0
    expected = align_motifs(*motifs, method=method)
    assert format_alignment(*motifs, expected) in result.output
    assert (
        f"method={method} score={expected.score:.4f} offset={expected.offset} "
        f"overlap={expected.overlap} gaps={expected.gaps} "
        f"is_reverse_complement={expected.is_reverse_complement}"
    ) in result.output


def test_align_uses_matrix_align_by_default(motifs):
    result = CliRunner().invoke(align, ["MA0001.1", "MA0002.1", "--format", "json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert set(data) == {
        "motif1",
        "motif2",
        "method",
        "score",
        "offset",
        "overlap",
        "gaps",
        "is_reverse_complement",
    }
    expected = align_motifs(*motifs)
    assert data["method"] == "matrix_align"
    assert data["score"] == round(expected.score, 6)
    assert data["offset"] == expected.offset


def test_align_no_reverse_complement():
    result = CliRunner().invoke(
        align, ["MA0001.1", "MA0002.1", "--no-reverse-complement", "--format", "json"]
    )
    assert result.exit_code == 0
    assert json.loads(result.output)["is_reverse_complement"] is False


def test_align_min_overlap_larger_than_a_motif_fails_cleanly():
    # MA2557.1 has 4 columns
    result = CliRunner().invoke(
        align, ["MA0139.2", "MA2557.1", "--method", "pearson", "--min-overlap", "5"]
    )
    assert result.exit_code == 1
    assert "fewer than min_overlap=5" in result.output


def test_align_min_overlap_with_matrix_align_fails_cleanly():
    result = CliRunner().invoke(align, ["MA0001.1", "MA0002.1", "--min-overlap", "5"])
    assert result.exit_code == 1
    assert "min_overlap" in result.output


def test_align_unknown_motif():
    result = CliRunner().invoke(align, ["MA9999.9", "MA0002.1"])
    assert result.exit_code == 1
    assert "No motif found" in result.output
