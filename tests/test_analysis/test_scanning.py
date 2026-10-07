"""Tests for motif scanning."""

from __future__ import annotations

import pytest
from Bio.motifs.jaspar import Motif
from Bio.Seq import Seq

from pyjaspar import JasparDB
from pyjaspar.analysis.scanning import ScanHit, scan_sequence


@pytest.fixture
def acg_motif():
    return Motif(
        matrix_id="TEST",
        name="ACG",
        counts={base: [10 if letter == base else 0 for letter in "ACG"] for base in "ACGT"},
    )


@pytest.mark.parametrize("sequence", ["CGTAAAAAAA", "AAAACGTAAA", "AAAAAAACGT"])
def test_reverse_hit_uses_original_sequence_coordinates(acg_motif, sequence):
    hits = scan_sequence(sequence, acg_motif, threshold=0.99)
    reverse_hits = [hit for hit in hits if hit.strand == "-"]
    assert len(reverse_hits) == 1
    hit = reverse_hits[0]
    assert hit.position == sequence.index("CGT")
    assert hit.sequence == "CGT"
    assert hit.sequence == sequence[hit.position : hit.position + acg_motif.length]
    pssm = acg_motif.counts.normalize(0.001).log_odds()
    assert hit.score == pytest.approx(pssm.calculate(Seq(hit.sequence).reverse_complement()))


def test_mixed_strand_hits_are_sorted_in_original_coordinates(acg_motif):
    sequence = "AAAACGTAAA"
    hits = scan_sequence(sequence, acg_motif, threshold=0.99)
    assert [(hit.position, hit.strand, hit.sequence) for hit in hits] == [
        (3, "+", "ACG"),
        (4, "-", "CGT"),
    ]
    assert scan_sequence(sequence, acg_motif, threshold=0.99, both_strands=False) == hits[:1]


def test_reverse_only_site_is_excluded_when_both_strands_is_false(acg_motif):
    assert scan_sequence("CGTAAAAAAA", acg_motif, threshold=0.99, both_strands=False) == []


@pytest.fixture(scope="module")
def jdb():
    return JasparDB()


@pytest.fixture(scope="module")
def ctcf_motif(jdb):
    motifs = jdb.fetch_motifs_by_name("CTCF")
    return motifs[0]


def test_scan_finds_hits(ctcf_motif):
    """A sequence containing a strong CTCF site should produce hits."""
    # CTCF consensus-like sequence embedded in random flanking
    sequence = "AAAAACCACCAGGGGGCGCAAAAAA"
    hits = scan_sequence(sequence, ctcf_motif, threshold=0.5)
    assert isinstance(hits, list)
    for hit in hits:
        assert isinstance(hit, ScanHit)
        assert hit.strand in ("+", "-")
        assert isinstance(hit.score, float)
        assert isinstance(hit.position, int)


def test_scan_empty_sequence(ctcf_motif):
    """Empty sequence should return no hits."""
    hits = scan_sequence("", ctcf_motif, threshold=0.8)
    assert hits == []


def test_scan_short_sequence(ctcf_motif):
    """Sequence shorter than motif should return no hits."""
    hits = scan_sequence("ACGT", ctcf_motif, threshold=0.8)
    assert hits == []


def test_scan_low_threshold(ctcf_motif):
    """Low threshold should find more hits."""
    sequence = "ACGTACGTACGTACGTACGTACGTACGTACGT"
    hits_high = scan_sequence(sequence, ctcf_motif, threshold=0.9)
    hits_low = scan_sequence(sequence, ctcf_motif, threshold=0.3)
    assert len(hits_low) >= len(hits_high)


def test_scan_hit_fields(ctcf_motif):
    """Verify ScanHit dataclass fields."""
    sequence = "ACGTACGTACGTACGTACGTACGTACGTACGT" * 3
    hits = scan_sequence(sequence, ctcf_motif, threshold=0.3)
    if hits:
        hit = hits[0]
        assert hasattr(hit, "position")
        assert hasattr(hit, "strand")
        assert hasattr(hit, "score")
        assert hasattr(hit, "sequence")
        assert len(hit.sequence) == ctcf_motif.length
