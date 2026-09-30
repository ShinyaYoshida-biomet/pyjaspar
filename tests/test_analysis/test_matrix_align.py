"""Tests for searching profiles with a query motif.

Expected scores come from two public sources, recorded next to each value:
the TFBSTools manual example, and the "Score" column of the JASPAR web tool
(https://jaspar.elixir.no/align/, queried on 2026-09-30 with collection CORE,
taxonomic group vertebrates, latest versions).
"""

from __future__ import annotations

import pytest
from Bio.motifs.jaspar import Motif

from pyjaspar import JasparDB
from pyjaspar.analysis.matrix_align import ProfileHit, align_score, search_profiles


@pytest.fixture(scope="module")
def jdb():
    return JasparDB("JASPAR2026")


@pytest.fixture(scope="module")
def egr1(jdb):
    return jdb.fetch_motif_by_id("MA0162.2")


@pytest.fixture(scope="module")
def ctcf(jdb):
    return jdb.fetch_motif_by_id("MA0139.2")


def with_flat_columns(motif, position, how_many):
    """Copy of ``motif`` with uninformative columns inserted at ``position``."""
    counts = {b: list(motif.counts[b]) for b in "ACGT"}
    total = sum(counts[b][position] for b in "ACGT")
    for b in "ACGT":
        counts[b][position:position] = [total / 4] * how_many
    return Motif(matrix_id="", name="flat", counts=counts)


def test_tfbstools_manual_example():
    # TFBSTools manual, PFMSimilarity(MA0003.2, MA0004.1) on JASPAR2014 data: score 7.294736
    db = JasparDB("JASPAR2014")
    a = db.fetch_motif_by_id("MA0003.2")
    b = db.fetch_motif_by_id("MA0004.1")
    assert align_score(a, b).score == pytest.approx(7.294736, abs=1e-4)
    assert align_score(b, a).score == pytest.approx(7.294736, abs=1e-4)


@pytest.mark.parametrize(
    ("matrix_id", "web_score"),
    [
        # JASPAR web tool, query = human EGR1 MA0162.2 (JASPAR2026 counts)
        ("MA0002.3", 13.1023),
        ("MA0003.5", 13.3877),
        ("MA0004.1", 8.64621),
        ("MA0006.2", 8.67515),
        ("MA0007.4", 17.2406),
        ("MA1723.2", 24.6726),
        ("MA0073.2", 23.9691),
        ("MA1929.2", 23.8812),
        ("MA0162.5", 18.9821),
        # the best alignments of these two contain a gap
        ("MA2457.1", 16.9322),
        ("MA1654.2", 20.3021),
    ],
)
def test_matches_web_tool_for_egr1_query(jdb, egr1, matrix_id, web_score):
    candidate = jdb.fetch_motif_by_id(matrix_id)
    assert align_score(egr1, candidate).score == pytest.approx(web_score, abs=1e-3)


def test_gapped_rows_use_a_gap(jdb, egr1):
    assert align_score(egr1, jdb.fetch_motif_by_id("MA2457.1")).gaps == 1
    assert align_score(egr1, jdb.fetch_motif_by_id("MA0002.3")).gaps == 0


@pytest.mark.parametrize(
    ("matrix_id", "web_score"),
    [
        # JASPAR web tool, query = CTCF MA0139.2: its own row is 2 * 15 columns
        ("MA0139.2", 30.0),
        ("MA1930.2", 29.888),
        ("MA1929.2", 27.8202),
        ("MA2510.1", 24.2936),
    ],
)
def test_matches_web_tool_for_ctcf_query(jdb, ctcf, matrix_id, web_score):
    candidate = jdb.fetch_motif_by_id(matrix_id)
    assert align_score(ctcf, candidate).score == pytest.approx(web_score, abs=1e-3)


@pytest.mark.parametrize(
    ("matrix_id", "web_score"),
    [
        # JASPAR web tool, query = CTCF MA0139.2 with two flat columns inserted after column 7;
        # its own row is 30 - (3 + 0.01): one gap of two columns
        ("MA0139.2", 26.99),
        ("MA1930.2", 28.9098),
        ("MA1929.2", 27.8626),
        ("MA1987.2", 27.3149),
    ],
)
def test_matches_web_tool_for_query_with_inserted_columns(jdb, ctcf, matrix_id, web_score):
    query = with_flat_columns(ctcf, position=7, how_many=2)
    candidate = jdb.fetch_motif_by_id(matrix_id)
    assert align_score(query, candidate).score == pytest.approx(web_score, abs=1e-3)


def test_self_alignment_scores_two_per_column(ctcf):
    result = align_score(ctcf, ctcf)
    assert result.score == pytest.approx(2 * ctcf.length)
    assert result.gaps == 0
    assert result.is_reverse_complement is False


def test_reverse_complement_does_not_change_the_score(jdb, egr1):
    candidate = jdb.fetch_motif_by_id("MA0002.3")
    forward = align_score(egr1, candidate)
    reverse = align_score(egr1, candidate.reverse_complement())
    assert forward.score == pytest.approx(reverse.score)
    assert forward.is_reverse_complement != reverse.is_reverse_complement


def test_inserted_columns_cost_the_gap_penalty(ctcf):
    one = align_score(with_flat_columns(ctcf, 7, 1), ctcf)
    two = align_score(with_flat_columns(ctcf, 7, 2), ctcf)
    assert one.gaps == two.gaps == 1
    assert one.score == pytest.approx(30.0 - 3.0)
    assert two.score == pytest.approx(30.0 - 3.0 - 0.01)


def test_gap_penalties_are_parameters(ctcf):
    query = with_flat_columns(ctcf, 7, 2)
    cheap = align_score(query, ctcf, open_penalty=1.0, ext_penalty=0.5)
    assert cheap.score == pytest.approx(30.0 - 1.0 - 0.5)


def test_column_without_counts_is_rejected(ctcf):
    counts = {b: list(ctcf.counts[b]) for b in "ACGT"}
    for b in "ACGT":
        counts[b][3] = 0
    broken = Motif(matrix_id="broken", name="broken", counts=counts)
    with pytest.raises(ValueError, match="no counts"):
        align_score(ctcf, broken)


def test_search_returns_hits_best_first(jdb, ctcf):
    candidates = [jdb.fetch_motif_by_id(i) for i in ("MA1929.2", "MA0139.2", "MA1930.2")]
    hits = search_profiles(ctcf, candidates)
    assert [h.matrix_id for h in hits] == ["MA0139.2", "MA1930.2", "MA1929.2"]
    assert all(isinstance(h, ProfileHit) for h in hits)
    assert hits[0].name == "CTCF"
    assert hits[0].score == pytest.approx(30.0)
    assert hits[0].width == 15


def test_search_relative_score_is_per_pair(jdb, ctcf):
    # CTCF MA1929.2 is 31 columns wide: 100 * score / (2 * min(15, 31))
    (hit,) = search_profiles(ctcf, [jdb.fetch_motif_by_id("MA1929.2")])
    assert hit.width == 31
    assert hit.relative_score == pytest.approx(100 * hit.score / (2 * 15))


def test_search_sort_by_relative_score(jdb, egr1):
    candidates = [jdb.fetch_motif_by_id(i) for i in ("MA0006.2", "MA1723.2", "MA0162.5")]
    by_score = search_profiles(egr1, candidates)
    by_relative = search_profiles(egr1, candidates, sort_by="relative_score")
    assert [h.score for h in by_score] == sorted((h.score for h in by_score), reverse=True)
    assert [h.relative_score for h in by_relative] == sorted(
        (h.relative_score for h in by_relative), reverse=True
    )


def test_search_top_truncates(jdb, ctcf):
    candidates = [jdb.fetch_motif_by_id(i) for i in ("MA1929.2", "MA0139.2", "MA1930.2")]
    assert len(search_profiles(ctcf, candidates, top=2)) == 2


def test_search_with_no_candidates(ctcf):
    assert search_profiles(ctcf, []) == []


def test_search_rejects_unknown_sort_key(ctcf):
    with pytest.raises(ValueError, match="sort_by"):
        search_profiles(ctcf, [], sort_by="pvalue")
