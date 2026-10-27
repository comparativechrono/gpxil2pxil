"""Gene, transcript and protein sequences from a gpxil file: the reference, or one of a sample's haplotypes.

DNA is kept as gpxil reference codes, which are also the indices of pxil's
standard DNA alphabet (``-ACMGRSVTWYHKDBN-acmgrsvtwyhkdbn``): soft-masked
bases (lower case) keep their mark, and ambiguity codes and N survive.

A haplotype's sequence applies the sample's variants as gpxiltools does
(:func:`gpxiltools.export.haplotype_sequence`): substitutions replace the
reference base, deletions drop it, insertions follow their anchor base, and
bases whose state is unknown become N. A transcript of a haplotype is the
haplotype's sequence over each exon; its protein is translated from the
annotated start codon until the first stop, which a frameshift can move
into the 3' UTR.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from gpxiltools import codes as C
from gpxiltools.export import haplotype_sequence

from .annotation import Gene, Transcript

DNA = "-ACMGRSVTWYHKDBN-acmgrsvtwyhkdbn"  # pxil's standard DNA alphabet = gpxil reference codes
_INDEX = np.full(256, 255, dtype=np.uint8)
for _i, _ch in enumerate(DNA):
    if _INDEX[ord(_ch)] == 255:
        _INDEX[ord(_ch)] = _i
_INDEX[ord("U")], _INDEX[ord("u")] = _INDEX[ord("T")], _INDEX[ord("t")]


def _complement_table() -> np.ndarray:
    out = np.zeros(256, dtype=np.uint8)
    for code in range(32):
        low, mask = code & 0x0F, code & 0x10
        swapped = ((low & 1) << 3) | ((low & 2) << 1) | ((low & 4) >> 1) | ((low & 8) >> 3)  # A<->T, C<->G
        out[code] = mask | swapped
    return out


_COMPLEMENT = _complement_table()


def letters_to_codes(seq: bytes | str) -> np.ndarray:
    raw = np.frombuffer(seq.encode() if isinstance(seq, str) else seq, dtype=np.uint8)
    codes = _INDEX[raw]
    if (codes == 255).any():
        bad = chr(int(raw[np.flatnonzero(codes == 255)[0]]))
        raise ValueError(f"{bad!r} is not a DNA symbol")
    return codes


def codes_to_letters(codes: np.ndarray) -> str:
    return "".join(DNA[int(c)] for c in codes) if codes.size < 64 else C.ref_to_ascii(codes).decode("ascii")


def reverse_complement(codes: np.ndarray) -> np.ndarray:
    return _COMPLEMENT[codes[::-1]]


# ----------------------------------------------------------------- genetic codes

_BASES = "TCAG"
_STANDARD = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_TABLES = {
    "standard": _STANDARD,
    # NCBI table 2: AGA and AGG are stops, ATA is Met, TGA is Trp
    "vertebrate mitochondrial": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
}
_CODONS = {
    name: {a + b + c: t[16 * i + 4 * j + k] for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}
    for name, t in _TABLES.items()
}


def genetic_code_for(contig: str) -> str:
    return "vertebrate mitochondrial" if contig.removeprefix("chr").upper() in ("M", "MT") else "standard"


def translate(dna: str, code: str = "standard", *, sec_codons: set[int] = frozenset(), until_stop: bool = True) -> str:
    """Translate DNA (letters, any case) from its first base. Codons with anything but A, C, G or T become X.

    ``sec_codons`` are codon numbers (0-based) read as selenocysteine (U) when they are TGA.
    """
    table = _CODONS[code]
    dna = dna.upper().replace("U", "T")
    out = []
    for k in range(len(dna) // 3):
        codon = dna[3 * k : 3 * k + 3]
        aa = "U" if (k in sec_codons and codon == "TGA") else table.get(codon, "X")
        out.append(aa)
        if aa == "*" and until_stop:
            break
    return "".join(out)


# ---------------------------------------------------------------- the sequences


@dataclass
class Spliced:
    """A transcript's spliced sequence (5' to 3', as gpxil codes) and where its coding sequence starts in it."""

    codes: np.ndarray
    cds_start: int | None  # offset of the first base of the CDS (after skipping ``frame``), or None
    exon_lengths: list[int] = field(default_factory=list)  # in transcript order

    def letters(self) -> str:
        return codes_to_letters(self.codes)


class Source:
    """Where sequences come from: the reference of a gpxil contig, or one haplotype of its sample (1 or 2)."""

    def __init__(self, gfile, haplotype: int | None = None):
        if haplotype not in (None, 1, 2):
            raise ValueError("haplotype must be 1, 2 or None (the reference)")
        self.g = gfile
        self.haplotype = haplotype
        self.sample = gfile.sample_id if haplotype else None

    @property
    def label(self) -> str:
        return f"{self.sample or 'sample'}-hap{self.haplotype}" if self.haplotype else "reference"

    def contig(self, name: str):
        if name not in self.g:
            raise KeyError(f"{name} is not in {self.g.path}")
        c = self.g[name]
        if self.haplotype and not c.has("gt"):
            raise KeyError(f"{c.name} has no genotypes in {self.g.path}, so it has no haplotype {self.haplotype}")
        return c

    def dna(self, contig: str, start: int, end: int) -> np.ndarray:
        """[start, end) on the forward strand, as codes."""
        c = self.contig(contig)
        if self.haplotype is None:
            return c.read("ref", start, end)
        return letters_to_codes(haplotype_sequence(c, self.haplotype, start, end, keep_case=True))

    def gene(self, gene: Gene) -> np.ndarray:
        """The gene's genomic sequence, 5' to 3' (introns and all)."""
        seq = self.dna(gene.contig, gene.start, gene.end)
        return reverse_complement(seq) if gene.strand == -1 else seq

    def spliced(self, t: Transcript) -> Spliced:
        exons = sorted(t.exons)
        pieces = [self.dna(t.contig, s, e) for s, e in exons]
        if t.strand == -1:
            pieces = [reverse_complement(p) for p in reversed(pieces)]
        codes = np.concatenate(pieces) if pieces else np.zeros(0, dtype=np.uint8)
        cds_start = None
        if t.cds:
            lo, hi = t.cds
            if t.strand == 1:
                before = [(s, min(e, lo)) for s, e in exons if s < lo]
            else:
                before = [(max(s, hi), e) for s, e in exons if e > hi]
            cds_start = sum(len(self.dna(t.contig, s, e)) for s, e in before if e > s) + t.frame
        return Spliced(codes, cds_start, [len(p) for p in pieces])

    def protein(self, t: Transcript, spliced: Spliced | None = None) -> str:
        """The protein from the annotated start until the first stop codon (with it, as ``*``), or to the transcript's end."""
        if not t.cds:
            raise ValueError(f"{t.name} ({t.id}) is not protein coding")
        sp = spliced or self.spliced(t)
        cds = sp.letters()[sp.cds_start :]
        sec = self._sec_codons(t) if t.sec else frozenset()
        return translate(cds, genetic_code_for(t.contig), sec_codons=sec)

    def _sec_codons(self, t: Transcript) -> set[int]:
        """Codon numbers of the annotated selenocysteines, counted from the CDS start in this sequence."""
        exons = sorted(t.exons)
        lo, hi = t.cds
        out = set()
        for p in t.sec:
            if t.strand == 1:
                parts = [(max(s, lo), min(e, p)) for s, e in exons]
            else:
                parts = [(max(s, p + 1), min(e, hi)) for s, e in exons]
            n = sum(len(self.dna(t.contig, s, e)) for s, e in parts if e > s) - t.frame
            if n >= 0 and n % 3 == 0:
                out.add(n // 3)
        return out


def protein_differences(ref: str, alt: str) -> list[str]:
    """Differences between two proteins: substitutions (P72R) where they line up, then a note on any change of length."""
    out = []
    n = min(len(ref), len(alt))
    for i in range(n):
        if ref[i] != alt[i]:
            if ref[i] == "*" or alt[i] == "*":
                out.append(f"{ref[i]}{i + 1}{alt[i]}")
                break
            out.append(f"{ref[i]}{i + 1}{alt[i]}")
            if len(out) > 50:
                out.append("…")
                break
    if len(ref) != len(alt):
        out.append(f"length {len(ref)} → {len(alt)}")
    return out


def haplotype_changes(gfile, contig: str, h: int, spans: list[tuple[int, int]]) -> int:
    """Positions where haplotype ``h`` differs from the reference within the spans (substitutions, deleted bases, insertions)."""
    c = gfile[contig]
    if not c.has("gt"):
        return 0
    n = 0
    for s, e in spans:
        gt = c.read("gt", s, e)
        st = C.normalize_states(C.hap1(gt) if h == 1 else C.hap2(gt))
        n += int(np.isin(st, C.HAP_SUB_CODES + (C.HAP_DEL,)).sum())
        if c.has("flags"):
            fl = c.read("flags", s, e)
            n += int(((fl & (C.FLAG_INS1 if h == 1 else C.FLAG_INS2)) != 0).sum())
    return n
