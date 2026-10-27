"""Genes and their transcripts, from the Gpxil browser's gene files or from a GTF.

Coordinates are 0-based and half-open, on the contig's forward strand. A
transcript's coding region, ``cds``, spans its CDS *and* its stop codon,
as the Gpxil browser stores it; ``frame`` is the number of bases to skip
before the first whole codon (non-zero only when the start is incomplete).
"""

from __future__ import annotations

import gzip
import io
import json
import os
import re
import urllib.request
from dataclasses import dataclass, field

CANONICAL, MANE_SELECT, MANE_CLINICAL, BASIC = 1, 2, 4, 8
INCOMPLETE = {"cds_start_NF", "cds_end_NF", "mRNA_start_NF", "mRNA_end_NF"}

#: the Gpxil browser's gene files: Ensembl release 116 on GRCh38, every transcript
BROWSER_DATA = "https://comparativechrono.github.io/Gpxil-browser/data/grch38"
#: an environment variable naming another copy of them (a URL or a folder), used in place of BROWSER_DATA
BROWSER_DATA_ENV = "GPXIL_BROWSER_DATA"


def browser_data_base() -> str:
    """Where the Gpxil browser's files are read from: $GPXIL_BROWSER_DATA, else the browser's own site."""
    return os.environ.get(BROWSER_DATA_ENV) or BROWSER_DATA


class AnnotationError(LookupError):
    """A gene that cannot be found, or found more than once."""


@dataclass
class Transcript:
    id: str
    name: str
    biotype: str
    contig: str
    strand: int  # 1 or -1
    exons: list[tuple[int, int]]  # ascending
    cds: tuple[int, int] | None = None  # CDS and stop codon
    frame: int = 0
    protein: str | None = None  # Ensembl protein ID
    flags: int = 0
    incomplete: tuple[str, ...] = ()
    sec: tuple[int, ...] = ()  # selenocysteine codons: their first base, in transcript orientation

    @property
    def start(self) -> int:
        return self.exons[0][0]

    @property
    def end(self) -> int:
        return self.exons[-1][1]

    @property
    def canonical(self) -> bool:
        return bool(self.flags & CANONICAL)

    @property
    def mane(self) -> bool:
        return bool(self.flags & (MANE_SELECT | MANE_CLINICAL))

    @property
    def coding(self) -> bool:
        return self.cds is not None

    @property
    def length(self) -> int:
        """Length of the spliced transcript."""
        return sum(e - s for s, e in self.exons)

    def tags(self) -> list[str]:
        out = []
        if self.flags & CANONICAL:
            out.append("Ensembl canonical")
        if self.flags & MANE_SELECT:
            out.append("MANE Select")
        if self.flags & MANE_CLINICAL:
            out.append("MANE Plus Clinical")
        return out


@dataclass
class Gene:
    id: str
    name: str
    contig: str
    start: int
    end: int
    strand: int
    biotype: str
    description: str = ""
    transcripts: list[Transcript] = field(default_factory=list)
    source: str = ""

    def transcript(self, which: str) -> Transcript:
        """A transcript by Ensembl ID (version optional) or name (TP53-201)."""
        w = which.strip()
        for t in self.transcripts:
            if _same_id(t.id, w) or t.name.lower() == w.lower():
                return t
        raise AnnotationError(f"{self.name} has no transcript {which!r} (it has {', '.join(t.name for t in self.transcripts[:12])}…)")

    def canonical(self) -> Transcript | None:
        """The Ensembl canonical transcript, else the MANE Select one, else the longest coding one."""
        for test in (lambda t: t.canonical, lambda t: t.flags & MANE_SELECT):
            for t in self.transcripts:
                if test(t):
                    return t
        if not self.transcripts:
            return None
        return max(self.transcripts, key=lambda t: (t.coding, (t.cds[1] - t.cds[0]) if t.cds else 0, t.length))

    def select(self, which: str = "canonical") -> list[Transcript]:
        """``canonical``, ``mane`` (MANE Select and Plus Clinical), ``coding``, ``all``, or IDs and names separated by commas."""
        if which == "canonical":
            t = self.canonical()
            return [t] if t else []
        if which == "mane":
            return [t for t in self.transcripts if t.mane] or self.select("canonical")
        if which == "coding":
            return [t for t in self.transcripts if t.coding]
        if which == "all":
            return list(self.transcripts)
        return [self.transcript(w) for w in which.split(",") if w.strip()]


def _same_id(a: str, b: str) -> bool:
    return a.split(".")[0].upper() == b.split(".")[0].upper()


def _is_url(path: str) -> bool:
    return path.lower().startswith(("http://", "https://"))


# ----------------------------------------------------------- the Gpxil browser


class BrowserData:
    """The Gpxil browser's gene files: a URL (by default, the browser's own) or a folder holding a copy.

    One request finds a gene in the index (78,941 genes, about 1 MB); another
    fetches its contig's gene list, and a range request just its transcripts.
    """

    def __init__(self, base: str | None = None, *, timeout: float = 60):
        self.base = str(base or browser_data_base()).rstrip("/")
        self.timeout = timeout
        self._index: list | None = None
        self._aliases: dict | None = None
        self._rows: dict[str, list] = {}
        self._manifest: dict | None = None

    def __repr__(self) -> str:
        return f"<BrowserData {self.base}>"

    # reading --------------------------------------------------------------
    def _read(self, rel: str, offset: int | None = None, length: int | None = None) -> bytes:
        if _is_url(self.base):
            headers = {"User-Agent": "gpxil2pxil"}
            if offset is not None:
                headers["Range"] = f"bytes={offset}-{offset + length - 1}"
            req = urllib.request.Request(f"{self.base}/{rel}", headers=headers)
            with urllib.request.urlopen(req, timeout=self.timeout) as r:
                data = r.read()
                if offset is not None and r.status == 200:  # the server ignored the range
                    data = data[offset : offset + length]
            return data
        path = os.path.join(self.base, *rel.split("/"))
        with open(path, "rb") as fh:
            if offset is None:
                return fh.read()
            fh.seek(offset)
            return fh.read(length)

    def _json(self, rel: str):
        data = self._read(rel)
        if data[:2] == b"\x1f\x8b":
            data = gzip.decompress(data)
        return json.loads(data)

    @property
    def manifest(self) -> dict:
        if self._manifest is None:
            self._manifest = self._json("manifest.json")
        return self._manifest

    @property
    def source(self) -> str:
        m = self.manifest
        genes = (m.get("sources") or {}).get("genes", "")
        rel = re.search(r"release (\d+)", genes)
        return f"{m.get('assembly', '')}, Ensembl release {rel.group(1)}" if rel else f"{m.get('assembly', '')}, {genes}"

    def _load_index(self) -> None:
        if self._index is None:
            g = self.manifest.get("genes") or {}
            self._index = self._json(g.get("index", "genes/index.json.gz"))
            self._aliases = {a: syms for a, syms in self._json(g.get("aliases", "genes/aliases.json.gz"))}

    # finding --------------------------------------------------------------
    def find(self, query: str) -> Gene:
        """A gene by Ensembl ID (ENSG…, version optional), symbol (TP53) or alias."""
        self._load_index()
        q = query.strip()
        hits = [r for r in self._index if _same_id(r[1], q)] if q.upper().startswith("ENSG") else []
        if not hits:
            hits = [r for r in self._index if r[0].upper() == q.upper()]
        if not hits and self._aliases and q in self._aliases:
            syms = self._aliases[q]
            if len(syms) > 1:
                raise AnnotationError(f"{q} is an alias of {', '.join(syms)}: give one of them")
            hits = [r for r in self._index if r[0] == syms[0]]
        if not hits:
            raise AnnotationError(f"no gene {q!r} in the Gpxil browser's genes ({self.base})")
        if len(hits) > 1:
            primary = [r for r in hits if r[2].startswith("chr") and "_" not in r[2]]
            hits = primary if len(primary) == 1 else hits
        if len(hits) > 1:
            raise AnnotationError(f"{q} names {len(hits)} genes: " + ", ".join(f"{r[1]} ({r[2]}:{r[3] + 1})" for r in hits[:6]))
        return self._gene(hits[0][2], hits[0][1])

    def genes_in(self, contig: str, start: int, end: int) -> list[Gene]:
        """The genes overlapping a region."""
        rows = self._contig_rows(contig)
        return [self._gene_from_row(contig, r) for r in rows if r[2] < end and r[3] > start]

    def _contig_rows(self, contig: str) -> list:
        if contig not in self._rows:
            self._rows[contig] = self._json(f"genes/{contig}.json.gz")["genes"]
        return self._rows[contig]

    def _gene(self, contig: str, gene_id: str) -> Gene:
        for r in self._contig_rows(contig):
            if r[0] == gene_id:
                return self._gene_from_row(contig, r)
        raise AnnotationError(f"{gene_id} is in the index but not in genes/{contig}.json.gz")

    def _gene_from_row(self, contig: str, row: list) -> Gene:
        gid, name, start, end, strand, biotype, description = row[:7]
        where = row[9] if len(row) > 9 else None
        g = Gene(gid, name, contig, start, end, strand, biotype, description or "", source=self.source)
        if where:
            offset, length, _count = where
            rows = json.loads(gzip.decompress(self._read(f"transcripts/{contig}.bin", offset, length)))
            g.transcripts = [_browser_transcript(r, contig, start, strand) for r in rows]
        elif row[7]:  # older files: only the canonical transcript, in absolute coordinates
            g.transcripts = [_browser_canonical(row[7], contig, strand)]
        return g


def _browser_transcript(r: list, contig: str, g0: int, strand: int) -> Transcript:
    tid, name, biotype, flags, _s, _e, flat, cds, incomplete = r[:9]
    sec = r[9] if len(r) > 9 else []
    exons = [(g0 + flat[i], g0 + flat[i + 1]) for i in range(0, len(flat), 2)]
    return Transcript(
        id=tid,
        name=name,
        biotype=biotype,
        contig=contig,
        strand=strand,
        exons=sorted(exons),
        cds=(g0 + cds[0], g0 + cds[1]) if cds else None,
        frame=int(cds[2]) if cds else 0,
        protein=(cds[3] or None) if cds and len(cds) > 3 else None,
        flags=int(flags),
        incomplete=tuple(x for x in (incomplete or "").split(",") if x),
        sec=tuple(g0 + p for p in sec),
    )


def _browser_canonical(r: list, contig: str, strand: int) -> Transcript:
    tid, name, biotype, mane, ts, _te, flat, cds = r[:8]
    exons = [(ts + flat[i], ts + flat[i + 1]) for i in range(0, len(flat), 2)]
    return Transcript(tid, name, biotype, contig, strand, sorted(exons), (cds[0], cds[1]) if cds else None, int(cds[2]) if cds else 0,
                      (cds[3] or None) if cds and len(cds) > 3 else None, CANONICAL | (MANE_SELECT if mane else 0))  # fmt: skip


# --------------------------------------------------------------------------- GTF

_ATTR = re.compile(r'(\S+)\s+"([^"]*)"')


def _attrs(text: str) -> tuple[dict[str, str], set[str]]:
    a: dict[str, str] = {}
    tags: set[str] = set()
    for k, v in _ATTR.findall(text):
        if k == "tag":
            tags.add(v)
        else:
            a.setdefault(k, v)
    return a, tags


def _open_text(path) -> io.TextIOBase:
    fh = open(path, "rb")
    magic = fh.read(2)
    fh.seek(0)
    if magic == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.GzipFile(fileobj=fh), encoding="utf-8")
    return io.TextIOWrapper(fh, encoding="utf-8")


class GTF:
    """Genes from a GTF file (Ensembl or GENCODE), plain or gzipped.

    The file is read once for each call of :meth:`find_all`; Ensembl and
    GENCODE keep each gene's lines together, so reading stops after the last
    gene asked for (``grouped=False`` reads to the end).
    """

    def __init__(self, path, *, grouped: bool = True, source: str | None = None):
        self.path = os.fspath(path)
        self.grouped = grouped
        self.source = source or os.path.basename(self.path)

    def __repr__(self) -> str:
        return f"<GTF {self.path}>"

    def find(self, query: str) -> Gene:
        return self.find_all([query])[query]

    def find_all(self, queries: list[str]) -> dict[str, Gene]:
        want = {q.strip(): q for q in queries}
        genes: dict[str, Gene] = {}  # gene_id -> gene
        txs: dict[str, dict] = {}
        matched: dict[str, str] = {}  # query -> gene_id
        current = None
        with _open_text(self.path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                kind = f[2]
                if kind == "gene":
                    if self.grouped and current is not None and len(matched) == len(want):
                        break
                    a, _ = _attrs(f[8])
                    gid, name = a.get("gene_id", ""), a.get("gene_name", a.get("gene_id", ""))
                    hit = [q for q in want if _same_id(gid, q) or name.upper() == q.upper()]
                    current = None
                    if hit:
                        for q in hit:
                            if q in matched and matched[q] != gid:
                                raise AnnotationError(f"{q} names more than one gene in {self.path}: {matched[q]} and {gid}")
                            matched[q] = gid
                        current = gid
                        genes[gid] = Gene(gid, name, f[0], int(f[3]) - 1, int(f[4]), 1 if f[6] == "+" else -1,
                                          a.get("gene_biotype", a.get("gene_type", "")), source=self.source)  # fmt: skip
                    continue
                if current is None:
                    if self.grouped:
                        continue
                    a, _ = _attrs(f[8])
                    gid = a.get("gene_id")
                    if gid not in genes:
                        continue
                    current_line_gene = gid
                else:
                    current_line_gene = current
                a, tags = _attrs(f[8])
                if a.get("gene_id") != current_line_gene:
                    continue
                start0, end = int(f[3]) - 1, int(f[4])
                if kind == "transcript":
                    t = Transcript(
                        id=a["transcript_id"],
                        name=a.get("transcript_name", a["transcript_id"]),
                        biotype=a.get("transcript_biotype", a.get("transcript_type", "")),
                        contig=f[0],
                        strand=1 if f[6] == "+" else -1,
                        exons=[],
                        flags=(CANONICAL if "Ensembl_canonical" in tags else 0)
                        | (MANE_SELECT if "MANE_Select" in tags else 0)
                        | (MANE_CLINICAL if "MANE_Plus_Clinical" in tags else 0)
                        | (BASIC if "basic" in tags else 0),
                        incomplete=tuple(sorted(x for x in tags if x in INCOMPLETE)),
                    )
                    txs[t.id] = {"t": t, "frame": None, "sec": []}
                    genes[current_line_gene].transcripts.append(t)
                    continue
                rec = txs.get(a.get("transcript_id", ""))
                if rec is None:
                    continue
                t = rec["t"]
                if kind == "exon":
                    t.exons.append((start0, end))
                elif kind == "Selenocysteine":
                    rec["sec"].append(start0 if t.strand == 1 else end - 1)
                elif kind in ("CDS", "stop_codon"):
                    lo, hi = t.cds or (start0, end)
                    t.cds = (min(lo, start0), max(hi, end))
                    if kind == "CDS":
                        t.protein = t.protein or a.get("protein_id")
                        five = start0 if t.strand == 1 else -end
                        if rec["frame"] is None or five < rec["frame"][0]:
                            rec["frame"] = (five, int(f[7]) if f[7] not in (".", "") else 0)
        for rec in txs.values():
            t = rec["t"]
            t.exons.sort()
            t.frame = rec["frame"][1] if rec["frame"] else 0
            t.sec = tuple(sorted(rec["sec"]))
        missing = [q for q in want if q not in matched]
        if missing:
            raise AnnotationError(f"not in {self.path}: {', '.join(missing)}")
        return {q: genes[matched[q]] for q in want}


def annotation(source: str | None):
    """The annotation for a source: a GTF file, a folder or URL of the Gpxil browser's files, or None for the browser's own
    (or the copy that $GPXIL_BROWSER_DATA names)."""
    if source is None:
        return BrowserData()
    s = str(source)
    if _is_url(s) or os.path.isdir(s):
        return BrowserData(s)
    return GTF(s)
