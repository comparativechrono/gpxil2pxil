"""Writing pxil files of genes, transcripts and proteins.

from gpxiltools import GpxilFile
from gpxil2pxil import BrowserData, make

with GpxilFile("HG002.gpxil.tif") as g:
    gene = BrowserData().find("TP53")
    for f in make(g, gene, haplotypes=(None, 1, 2), out_dir="TP53"):
        print(f.path, f.notes)
"""

from __future__ import annotations

import math
import os
import re
from collections.abc import Iterable
from dataclasses import dataclass, field

import pxil

from . import __version__
from .annotation import Gene, Transcript
from .sequence import DNA, Source, codes_to_letters, genetic_code_for, haplotype_changes, protein_differences

KINDS = ("gene", "transcript", "protein")


@dataclass
class Written:
    """One pxil file written."""

    path: str
    kind: str
    name: str  # the gene's, transcript's or protein's label
    id: str
    length: int
    haplotype: int | None = None
    notes: list[str] = field(default_factory=list)


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(text))


def default_width(n: int, kind: str) -> int:
    """About square, in steps of 10 residues (proteins), 30 bases (transcripts: whole codons) or 100 bases (genes)."""
    step = {"gene": 100, "transcript": 30}.get(kind, 10)
    return max(step, math.ceil(math.sqrt(max(n, 1)) / step) * step)


def file_name(kind: str, gene: Gene, t: Transcript | None, source: Source) -> str:
    who = f"_{safe_name(source.label)}" if source.haplotype else ""
    if kind == "gene":
        return f"{safe_name(gene.name)}_{gene.id}{who}.gene.pxil.png"
    if kind == "transcript":
        return f"{safe_name(t.name)}_{t.id}{who}.transcript.pxil.png"
    return f"{safe_name(t.name)}_{t.protein or t.id + '_protein'}{who}.protein.pxil.png"


def make(
    gfile,
    gene: Gene,
    *,
    transcripts: str = "canonical",
    kinds: Iterable[str] = KINDS,
    haplotypes: Iterable[int | None] = (None,),
    out_dir: str = ".",
    width: int | None = None,
) -> list[Written]:
    """Write pxil files for a gene: its genomic sequence, and its transcripts' spliced sequences and proteins.

    ``haplotypes`` lists where the sequences come from: None for the reference, 1 and 2 for the sample's
    haplotypes. ``transcripts`` is ``canonical``, ``mane``, ``coding``, ``all`` or a comma-separated list
    of transcript IDs or names. Files are named with the Ensembl IDs and, for a haplotype, the sample.
    """
    asked = list(kinds)
    unknown = [k for k in asked if k not in KINDS]
    if unknown:
        raise ValueError(f"unknown kinds: {', '.join(unknown)} (choose from {', '.join(KINDS)})")
    kinds = [k for k in KINDS if k in asked]
    os.makedirs(out_dir, exist_ok=True)
    chosen = gene.select(transcripts) if ("transcript" in kinds or "protein" in kinds) else []
    assembly = (gfile.metadata.get("reference") or {}).get("name")
    common = {
        "source": gene.source or None,
        "assembly": assembly,
        "contig": gfile[gene.contig].name if gene.contig in gfile else gene.contig,  # the gpxil file's name for it
        "gene": gene.id,
        "generator": f"gpxil2pxil {__version__}",
    }
    written: list[Written] = []
    for h in haplotypes:
        src = Source(gfile, h)
        ref = Source(gfile, None) if h else None
        who = {"sample": src.sample, "haplotype": h} if h else {}
        if "gene" in kinds:
            codes = src.gene(gene)
            seq = codes_to_letters(codes)
            notes = []
            if h:
                n = haplotype_changes(gfile, gene.contig, h, [(gene.start, gene.end)])
                notes.append(f"{n} change{'s' if n != 1 else ''} from the reference")
            path = os.path.join(out_dir, file_name("gene", gene, None, src))
            pxil.write(
                path, seq, kind="gene", alphabet=DNA, width=width or default_width(len(seq), "gene"),
                id=gene.id, name=gene.name, start=gene.start, end=gene.end, strand=gene.strand,
                orientation="5' to 3' of the gene", description=gene.description or None, **common, **who,
                **({"notes": notes} if notes else {}),
            )  # fmt: skip
            written.append(Written(path, "gene", gene.name, gene.id, len(seq), h, notes))
        for t in chosen:
            sp = src.spliced(t)
            protein = src.protein(t, sp) if t.coding else None
            ref_protein = ref.protein(t) if (ref and t.coding) else None
            cds = None
            if protein is not None and sp.cds_start is not None:
                cds = [sp.cds_start, min(len(sp.codes), sp.cds_start + 3 * len(protein))]
            tnotes = []
            if h:
                n = haplotype_changes(gfile, t.contig, h, t.exons)
                tnotes.append(f"{n} change{'s' if n != 1 else ''} from the reference in its exons")
            if "transcript" in kinds:
                seq = sp.letters()
                path = os.path.join(out_dir, file_name("transcript", gene, t, src))
                pxil.write(
                    path, seq, kind="transcript", alphabet=DNA, width=width or default_width(len(seq), "transcript"),
                    id=t.id, name=t.name, transcript=t.id, protein=t.protein, strand=t.strand,
                    start=t.start, end=t.end, exons=[list(e) for e in t.exons], cds=cds, biotype=t.biotype,
                    orientation="5' to 3' of the transcript", tags=t.tags() or None, **common, **who,
                    **({"notes": tnotes} if tnotes else {}),
                )  # fmt: skip
                written.append(Written(path, "transcript", t.name, t.id, len(seq), h, tnotes))
            if "protein" in kinds and protein is not None:
                pnotes = []
                if not protein.endswith("*"):
                    pnotes.append("no stop codon before the end of the transcript")
                if ref_protein is not None:
                    diffs = protein_differences(ref_protein, protein)
                    pnotes.append("the same as the reference" if not diffs else "differs from the reference: " + ", ".join(diffs))
                path = os.path.join(out_dir, file_name("protein", gene, t, src))
                pxil.write(
                    path, protein, kind="protein", width=width or default_width(len(protein), "protein"),
                    id=t.protein or f"{t.id}_protein", name=t.name, transcript=t.id, protein=t.protein,
                    genetic_code=genetic_code_for(t.contig), **common, **who,
                    **({"notes": pnotes} if pnotes else {}),
                )  # fmt: skip
                written.append(Written(path, "protein", t.name, t.protein or t.id, len(protein), h, pnotes))
    return written
