"""gpxil2pxil on a made-up genome whose genes, transcripts, proteins and variants are known exactly."""

from __future__ import annotations

import contextlib
import gzip
import io
import json
import os
import random

import pxil
import pytest
from gpxiltools import GpxilFile
from gpxiltools.encode import encode

from gpxil2pxil import GTF, BrowserData, Source, make, translate
from gpxil2pxil.cli import main

COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")
STOPS = {"TAA", "TAG", "TGA"}


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def coding(rng: random.Random, n_codons: int, sec_at: int | None = None) -> str:
    """ATG, then n_codons - 2 codons that are not stops (TGA at sec_at), then TAA."""
    codons = ["ATG"]
    while len(codons) < n_codons - 1:
        c = "".join(rng.choice("ACGT") for _ in range(3))
        if c not in STOPS:
            codons.append(c)
    if sec_at is not None:
        codons[sec_at] = "TGA"
    return "".join(codons) + "TAA"


def to_genomic(exons, strand, offset):
    """The genomic position of a transcript offset (0-based)."""
    k = offset
    for a, b in sorted(exons) if strand == 1 else sorted(exons, reverse=True):
        if k < b - a:
            return a + k if strand == 1 else b - 1 - k
        k -= b - a
    raise IndexError(offset)


def changed_base(codon: str, rev: bool = False) -> str:
    """A base for the codon's first position that changes its amino acid (and is not a stop); forward-strand if rev."""
    aa = translate(codon)
    for nb in "ACGT":
        alt = (nb.translate(COMP) if rev else nb) + codon[1:]
        if alt[0] != codon[0] and translate(alt) not in (aa, "*"):
            return nb
    raise ValueError(codon)


class Genome:
    """chr7: GENEA on the plus strand (with a selenocysteine), GENEB on the minus strand."""

    def __init__(self, seed=1):
        rng = random.Random(seed)
        self.seq = [rng.choice("ACGT") for _ in range(6000)]
        for i in range(1700, 1900):  # a soft-masked stretch inside GENEA's last intron... and exon
            self.seq[i] = self.seq[i].lower()
        # GENEA, + strand: exons [1000,1100) [1300,1420) [1600,1750); CDS from 1050 to 1660 (stop codon included)
        self.a_exons = [(1000, 1100), (1300, 1420), (1600, 1750)]
        cds = coding(rng, 76, sec_at=20)  # 228 bases
        self.a_cds_seq = cds
        spliced = list("".join(rng.choice("ACGT") for _ in range(370)))
        spliced[50 : 50 + len(cds)] = cds  # the CDS starts 50 bases into the transcript
        self.a_spliced = "".join(spliced)
        self._place(self.a_exons, self.a_spliced, 1)
        self.a_cds = (1050, to_genomic(self.a_exons, 1, 50 + len(cds) - 1) + 1)
        self.a_spliced = self._spliced(self.a_exons, 1)  # with the soft-mask case of the genome
        # GENEB, - strand: exons [3000,3150) [3400,3500); transcript 5' end at 3500
        self.b_exons = [(3000, 3150), (3400, 3500)]
        cds = coding(rng, 51)  # 153 bases: 60 in the first exon (3400..3460 reversed) ... the rest in the second
        spliced = list("".join(rng.choice("ACGT") for _ in range(250)))
        spliced[40 : 40 + len(cds)] = cds  # starts 40 bases into the transcript: genomic 3460
        self.b_cds_seq = cds
        self._place(self.b_exons, "".join(spliced), -1)
        self.b_spliced = self._spliced(self.b_exons, -1)
        self.b_cds = (to_genomic(self.b_exons, -1, 40 + len(cds) - 1), to_genomic(self.b_exons, -1, 40) + 1)
        self.text = "".join(self.seq)

    def _place(self, exons, spliced, strand):
        s = spliced if strand == 1 else rc(spliced)
        k = 0
        for a, b in sorted(exons):
            for i in range(a, b):
                self.seq[i] = s[k].lower() if self.seq[i].islower() else s[k]
                k += 1

    def _spliced(self, exons, strand):
        s = "".join("".join(self.seq[a:b]) for a, b in sorted(exons))
        return s if strand == 1 else rc(s)


def gtf_lines(g: Genome) -> str:
    def line(kind, a, b, strand, attrs, frame="."):
        return f"7\ttest\t{kind}\t{a + 1}\t{b}\t.\t{strand}\t{frame}\t{attrs}"

    out = ["#!genome-build test"]
    ga = 'gene_id "ENSG00000000001"; gene_version "1"; gene_name "GENEA"; gene_biotype "protein_coding";'
    ta = (
        ga
        + ' transcript_id "ENST00000000001"; transcript_name "GENEA-201"; transcript_biotype "protein_coding"; tag "basic"; tag "Ensembl_canonical"; tag "MANE_Select";'
    )
    out.append(line("gene", 1000, 1750, "+", ga))
    out.append(line("transcript", 1000, 1750, "+", ta))
    for a, b in g.a_exons:
        out.append(line("exon", a, b, "+", ta))
    lo, hi = g.a_cds
    for a, b in g.a_exons:
        s, e = max(a, lo), min(b, hi - 3)
        if s < e:
            out.append(line("CDS", s, e, "+", ta + ' protein_id "ENSP00000000001";', "0"))
    sec = to_genomic(g.a_exons, 1, 50 + 60)  # codon 20: 60 bases into the CDS
    out.append(line("Selenocysteine", sec, sec + 3, "+", ta))
    out.append(line("stop_codon", hi - 3, hi, "+", ta))
    gb = 'gene_id "ENSG00000000002"; gene_name "GENEB"; gene_biotype "protein_coding";'
    tb = gb + ' transcript_id "ENST00000000002"; transcript_name "GENEB-201"; transcript_biotype "protein_coding"; tag "Ensembl_canonical";'
    tb2 = gb + ' transcript_id "ENST00000000003"; transcript_name "GENEB-202"; transcript_biotype "retained_intron";'
    out.append(line("gene", 3000, 3500, "-", gb))
    out.append(line("transcript", 3000, 3500, "-", tb))
    for a, b in reversed(g.b_exons):
        out.append(line("exon", a, b, "-", tb))
    lo, hi = g.b_cds
    for a, b in reversed(g.b_exons):
        s, e = max(a, lo + 3), min(b, hi)
        if s < e:
            out.append(line("CDS", s, e, "-", tb + ' protein_id "ENSP00000000002";', "0"))
    out.append(line("stop_codon", lo, lo + 3, "-", tb))
    out.append(line("transcript", 3000, 3500, "-", tb2))
    out.append(line("exon", 3000, 3500, "-", tb2))
    return "\n".join(out) + "\n"


def browser_files(folder, gtf_path):
    """The same genes in the Gpxil browser's format (what its pipeline writes)."""
    genes = GTF(gtf_path, grouped=False).find_all(["GENEA", "GENEB"])
    os.makedirs(os.path.join(folder, "genes"))
    os.makedirs(os.path.join(folder, "transcripts"))
    rows, index, blob = [], [], io.BytesIO()
    for g in genes.values():
        g0 = g.start
        trs = []
        for t in g.transcripts:
            flat = [x - g0 for e in t.exons for x in e]
            cds = [t.cds[0] - g0, t.cds[1] - g0, t.frame, t.protein or ""] if t.cds else None
            r = [t.id, t.name, t.biotype, t.flags, t.start - g0, t.end - g0, flat, cds, ",".join(t.incomplete)]
            if t.sec:
                r.append([p - g0 for p in t.sec])
            trs.append(r)
        member = gzip.compress(json.dumps(trs).encode())
        where = [blob.tell(), len(member), len(trs)]
        blob.write(member)
        rows.append([g.id, g.name, g.start, g.end, g.strand, g.biotype, "", None, None, where])
        index.append([g.name, g.id, "chr7", g.start, g.end, g.strand, g.biotype])
    with gzip.open(os.path.join(folder, "genes", "chr7.json.gz"), "wt") as fh:
        json.dump({"contig": "chr7", "genes": rows}, fh)
    with open(os.path.join(folder, "transcripts", "chr7.bin"), "wb") as fh:
        fh.write(blob.getvalue())
    with gzip.open(os.path.join(folder, "genes", "index.json.gz"), "wt") as fh:
        json.dump(index, fh)
    with gzip.open(os.path.join(folder, "genes", "aliases.json.gz"), "wt") as fh:
        json.dump([["GA", ["GENEA"]]], fh)
    with open(os.path.join(folder, "manifest.json"), "w") as fh:
        json.dump(
            {
                "assembly": "TEST",
                "sources": {"genes": "Ensembl release 999: test.gtf"},
                "genes": {"index": "genes/index.json.gz", "aliases": "genes/aliases.json.gz"},
            },
            fh,
        )


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("world")
    g = Genome()
    fa = tmp / "ref.fa"
    fa.write_text(">7\n" + g.text + "\n")
    gtf = tmp / "genes.gtf"
    gtf.write_text(gtf_lines(g))
    # variants: GENEA codon 5 changes on haplotype 1; a 1-base deletion in codon 30 on haplotype 2 (a frameshift);
    # GENEB codon 3 changes on both; an insertion in GENEA's first intron on haplotype 1
    a_codon5 = to_genomic(g.a_exons, 1, 50 + 4 * 3)  # codon 5, first base
    old = g.text[a_codon5].upper()
    new = changed_base(g.a_cds_seq[12:15])
    a_del = to_genomic(g.a_exons, 1, 50 + 90)  # a base of codon 31, in exon 2: deleted on haplotype 2
    b_codon3 = to_genomic(g.b_exons, -1, 40 + 6)  # codon 3 of GENEB, first base (minus strand)
    bold = g.text[b_codon3].upper()
    bnew = changed_base(g.b_cds_seq[6:9], rev=True)
    vcf = tmp / "s.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n##contig=<ID=7,length=6000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        f"7\t{a_codon5 + 1}\t.\t{old}\t{new}\t50\tPASS\t.\tGT\t1|0\n"
        f"7\t{1150 + 1}\t.\t{g.text[1150].upper()}\t{g.text[1150].upper()}TTT\t50\tPASS\t.\tGT\t1|0\n"
        f"7\t{a_del}\t.\t{g.text[a_del - 1 : a_del + 1].upper()}\t{g.text[a_del - 1].upper()}\t50\tPASS\t.\tGT\t0|1\n"
        f"7\t{b_codon3 + 1}\t.\t{bold}\t{bnew}\t50\tPASS\t.\tGT\t1|1\n"
    )
    gpx = tmp / "s.gpxil.tif"
    encode(fa, gpx, vcf=vcf, ucsc_names=True, log=lambda m: None)
    ref = tmp / "ref.gpxil.tif"
    encode(fa, ref, ucsc_names=True, log=lambda m: None)
    browser = tmp / "browser"
    browser_files(str(browser), gtf)
    return {"g": g, "gtf": gtf, "gpx": gpx, "ref": ref, "browser": browser, "tmp": tmp, "a_codon5": a_codon5, "a_del": a_del, "new": new}


# ------------------------------------------------------------------ annotation


def test_gtf_and_browser_files_give_the_same_genes(world):
    a = GTF(world["gtf"]).find_all(["GENEA", "ENSG00000000002"])
    b = BrowserData(str(world["browser"]))
    for q, gene in a.items():
        other = b.find(q)
        assert other.id == gene.id and other.strand == gene.strand
        assert [(t.id, t.exons, t.cds, t.frame, t.protein, t.flags, t.sec) for t in other.transcripts] == [
            (t.id, t.exons, t.cds, t.frame, t.protein, t.flags, t.sec) for t in gene.transcripts
        ]
    assert b.find("GA").name == "GENEA"  # an alias
    assert a["GENEA"].canonical().name == "GENEA-201" and a["GENEA"].transcripts[0].sec == (to_genomic(world["g"].a_exons, 1, 110),)


def test_selecting_transcripts(world):
    gene = GTF(world["gtf"]).find("GENEB")
    assert [t.name for t in gene.select("canonical")] == ["GENEB-201"]
    assert [t.name for t in gene.select("all")] == ["GENEB-201", "GENEB-202"]
    assert [t.name for t in gene.select("coding")] == ["GENEB-201"]
    assert [t.id for t in gene.select("GENEB-202,ENST00000000002")] == ["ENST00000000003", "ENST00000000002"]


# -------------------------------------------------------------------- sequences


def test_reference_sequences_are_the_designed_ones(world):
    g = world["g"]
    genes = GTF(world["gtf"]).find_all(["GENEA", "GENEB"])
    with GpxilFile(world["ref"]) as gx:
        src = Source(gx)
        a, b = genes["GENEA"], genes["GENEB"]
        from gpxil2pxil.sequence import codes_to_letters

        assert codes_to_letters(src.gene(a)) == g.text[1000:1750]
        assert codes_to_letters(src.gene(b)) == rc(g.text[3000:3500])
        sa = src.spliced(a.transcripts[0])
        assert sa.letters() == g.a_spliced and sa.cds_start == 50
        sb = src.spliced(b.transcripts[0])
        assert sb.letters() == g.b_spliced and sb.cds_start == 40
        pa = src.protein(a.transcripts[0])
        assert pa == translate(g.a_cds_seq, sec_codons={20}) and pa[20] == "U" and pa.endswith("*")
        assert len(pa) == len(g.a_cds_seq) // 3
        assert src.protein(b.transcripts[0]) == translate(g.b_cds_seq)


def test_haplotypes(world):
    g = world["g"]
    genes = GTF(world["gtf"]).find_all(["GENEA", "GENEB"])
    a, b = genes["GENEA"].transcripts[0], genes["GENEB"].transcripts[0]
    with GpxilFile(world["gpx"]) as gx:
        ref = Source(gx).protein(a)
        h1 = Source(gx, 1).protein(a)
        h2 = Source(gx, 2).protein(a)
        assert h1 != ref and sum(x != y for x, y in zip(h1, ref, strict=False)) == 1 and h1[4] != ref[4] and len(h1) == len(ref)
        # haplotype 2 lost one base in codon 30: the same up to it, then a new frame until a new stop
        assert h2[:29] == ref[:29] and h2 != ref
        t2 = Source(gx, 2).spliced(a).letters()
        k = 50 + 90  # the deleted base, in the transcript
        assert t2 == g.a_spliced[:k] + g.a_spliced[k + 1 :]
        assert h2 == translate(t2[50:], sec_codons={20})
        # the insertion in the intron changes the gene, not the transcript
        assert len(Source(gx, 1).gene(genes["GENEA"])) == 750 + 3
        assert Source(gx, 1).spliced(a).letters() != g.a_spliced  # (the codon 5 change)
        assert len(Source(gx, 1).spliced(a).codes) == len(g.a_spliced)
        # GENEB: the same change on both haplotypes, on the minus strand
        rb = Source(gx).protein(b)
        for h in (1, 2):
            pb = Source(gx, h).protein(b)
            assert len(pb) == len(rb) and [i for i, (x, y) in enumerate(zip(pb, rb, strict=False)) if x != y] == [2]


# ---------------------------------------------------------------------- files


def test_files_named_with_ids_and_read_back(world, tmp_path):
    genes = GTF(world["gtf"]).find_all(["GENEA"])
    with GpxilFile(world["gpx"]) as gx:
        written = make(gx, genes["GENEA"], haplotypes=(None, 1, 2), out_dir=str(tmp_path))
    names = sorted(os.path.basename(w.path) for w in written)
    assert names == sorted(
        [
            "GENEA_ENSG00000000001.gene.pxil.png",
            "GENEA_ENSG00000000001_S1-hap1.gene.pxil.png",
            "GENEA_ENSG00000000001_S1-hap2.gene.pxil.png",
            "GENEA-201_ENST00000000001.transcript.pxil.png",
            "GENEA-201_ENST00000000001_S1-hap1.transcript.pxil.png",
            "GENEA-201_ENST00000000001_S1-hap2.transcript.pxil.png",
            "GENEA-201_ENSP00000000001.protein.pxil.png",
            "GENEA-201_ENSP00000000001_S1-hap1.protein.pxil.png",
            "GENEA-201_ENSP00000000001_S1-hap2.protein.pxil.png",
        ]
    )
    for w in written:
        seq, meta = pxil.read(w.path)
        assert len(seq) == w.length == meta["length"] and meta["kind"] == w.kind
        assert meta["gene"] == "ENSG00000000001" and meta["contig"] == "chr7"
        if w.haplotype:
            assert meta["sample"] == "S1" and meta["haplotype"] == w.haplotype
    p1 = next(w for w in written if w.kind == "protein" and w.haplotype == 1)
    assert any(n.startswith("differs from the reference: ") for n in p1.notes)
    seq, meta = pxil.read(next(w.path for w in written if w.kind == "transcript" and w.haplotype is None))
    assert seq == world["g"].a_spliced and meta["cds"] == [50, 278] and meta["exons"][0] == [1000, 1100]


def test_command_line(world, tmp_path):
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        assert main([str(world["gpx"]), "GENEB", "-a", str(world["gtf"]), "-o", str(tmp_path), "--kinds", "protein"]) == 0
    text = out.getvalue()
    assert "GENEB-201_ENSP00000000002_S1-hap1.protein.pxil.png" in text and "2 pxil files" in text
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        assert main([str(world["ref"]), "GENEA", "GENEB", "-a", str(world["browser"]), "-o", str(tmp_path / "ref"), "-t", "all"]) == 0
    assert out.getvalue().strip().endswith("7 pxil files")  # 2 genes, 3 transcripts, 2 proteins
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        assert main([str(world["gpx"]), "GENEB", "-a", str(world["gtf"]), "--list"]) == 0
    assert "* GENEB-201" in out.getvalue() and "GENEB-202" in out.getvalue()
    err = io.StringIO()
    with contextlib.redirect_stderr(err):
        assert main([str(world["gpx"]), "NOPE", "-a", str(world["gtf"])]) == 1
        assert main([str(world["ref"]), "GENEA", "-a", str(world["gtf"]), "--haplotypes", "1"]) == 1
    assert "NOPE" in err.getvalue() and "no sample" in err.getvalue()


@pytest.mark.skipif(not os.environ.get("GPXIL2PXIL_NETWORK"), reason="set GPXIL2PXIL_NETWORK=1 to use the Gpxil browser's files")
def test_tp53_from_the_gpxil_browser(tmp_path):
    with GpxilFile("https://comparativechrono.github.io/Gpxil-browser/data/grch38/gpxil/chr17.gpxil.tif") as gx:
        gene = BrowserData().find("TP53")
        written = make(gx, gene, kinds=("protein",), out_dir=str(tmp_path))
    seq, meta = pxil.read(written[0].path)
    assert seq.startswith("MEEPQSDPSV") and len(seq) == 394 and meta["protein"] == "ENSP00000269305"
