"""Command line: ``gpxil2pxil FILE GENE [GENE ...]``."""

from __future__ import annotations

import argparse
import os
import sys

from . import __version__
from .annotation import BROWSER_DATA, BROWSER_DATA_ENV, AnnotationError, annotation
from .build import KINDS, make


def _haplotypes(text: str, has_sample: bool) -> list[int | None]:
    if text == "auto":
        return [1, 2] if has_sample else [None]
    out: list[int | None] = []
    for part in text.split(","):
        part = part.strip().lower()
        if part in ("ref", "reference"):
            out.append(None)
        elif part in ("1", "2"):
            out.append(int(part))
        elif part == "both":
            out += [1, 2]
        elif part == "all":
            out += [None, 1, 2]
        else:
            raise SystemExit(f"gpxil2pxil: error: --haplotypes takes ref, 1, 2, both, all or auto, not {part!r}")
    return list(dict.fromkeys(out))


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="gpxil2pxil",
        description="Pictures (pxil files) of genes, transcripts and proteins from a gpxil genome: its reference, or a sample's haplotypes.",
        epilog=(
            "examples:\n"
            "  gpxil2pxil https://comparativechrono.github.io/Gpxil-browser/data/grch38/gpxil/chr17.gpxil.tif TP53\n"
            "  gpxil2pxil HG002.gpxil.tif TP53 BRCA1 --haplotypes all -o pictures\n"
            "  gpxil2pxil HG002.gpxil.tif TP53 --transcripts all --kinds protein\n"
            "  gpxil2pxil genome.gpxil.tif MYGENE --annotation genes.gtf.gz\n"
            "  gpxil2pxil HG002.gpxil.tif TP53 --list"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--version", action="version", version=f"gpxil2pxil {__version__}")
    p.add_argument("file", help="a gpxil file (path or URL)")
    p.add_argument("genes", nargs="+", help="genes: symbols (TP53) or Ensembl IDs (ENSG00000141510)")
    p.add_argument(
        "-t",
        "--transcripts",
        default="canonical",
        help="canonical (default), mane, coding, all, or transcript IDs or names separated by commas",
    )
    p.add_argument("-k", "--kinds", default=",".join(KINDS), help=f"what to write, separated by commas (default {','.join(KINDS)})")
    p.add_argument(
        "--haplotypes",
        default="auto",
        help="ref, 1, 2, both or all, separated by commas; auto (default): both haplotypes for a file with a sample, else the reference",
    )
    p.add_argument(
        "-a",
        "--annotation",
        help=f"a GTF file, or a folder or URL of the Gpxil browser's gene files (default: ${BROWSER_DATA_ENV} if set, else {BROWSER_DATA})",
    )
    p.add_argument("-o", "--output", default=".", help="folder for the files (default: here)")
    p.add_argument("--width", type=int, help="residues per row (default: about square)")
    p.add_argument("--reference", help="the reference (FASTA or gpxil, path or URL) for a file written without one")
    p.add_argument("--list", action="store_true", help="list each gene's transcripts and stop")
    return p


def main(argv=None) -> int:
    try:
        code = _main(argv)
        sys.stdout.flush()  # (here, so that a reader that has gone is noticed here, not as Python exits)
        return code
    except BrokenPipeError:
        # the reader stopped early (gpxil2pxil … | head): not an error of ours, so no traceback. Python's own
        # advice: point stdout at /dev/null, so that the flush at exit doesn't fail again
        try:
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, sys.stdout.fileno())
        except (OSError, ValueError, AttributeError):
            pass
        return 1


def _contig_names(path, reference):
    """the file's own name for each contig name it knows (17 → chr17), to list genes as the file names them"""
    from gpxiltools import GpxilFile

    try:
        with GpxilFile(path, reference=reference) as g:
            return lambda name: g[name].name if name in g else name
    except Exception:  # (a file that can't be opened just now: the annotation's names)
        return lambda name: name


def _main(argv=None) -> int:
    from gpxiltools import GpxilFile

    a = build_parser().parse_args(argv)
    ann = annotation(a.annotation)
    try:
        genes = ann.find_all(a.genes) if hasattr(ann, "find_all") else {q: ann.find(q) for q in a.genes}
    except AnnotationError as exc:
        print(f"gpxil2pxil: error: {exc}", file=sys.stderr)
        return 1
    if a.list:
        named = _contig_names(a.file, a.reference)
        for _q, g in genes.items():
            print(
                f"{g.name} {g.id} {named(g.contig)}:{g.start + 1:,}-{g.end:,} ({'+' if g.strand == 1 else '-'}), {g.biotype}, {len(g.transcripts)} transcripts"
            )
            canon = g.canonical()
            for t in g.transcripts:
                mark = "*" if t is canon else " "
                print(f" {mark} {t.name:<14} {t.id:<18} {t.biotype:<24} {t.protein or '':<18} {', '.join(t.tags())}")
        return 0
    try:
        with GpxilFile(a.file, reference=a.reference) as g:
            haps = _haplotypes(a.haplotypes, bool(g.sample_id))
            if any(haps) and not g.sample_id:
                print(f"gpxil2pxil: error: {a.file} has no sample, so it has no haplotypes; use --haplotypes ref", file=sys.stderr)
                return 1
            total = 0
            for _q, gene in genes.items():
                files = make(g, gene, transcripts=a.transcripts, kinds=a.kinds.split(","), haplotypes=haps, out_dir=a.output, width=a.width)
                for f in files:
                    print(
                        f"wrote {f.path} ({f.length:,} {'amino acids' if f.kind == 'protein' else 'bases'})"
                        + (f": {'; '.join(f.notes)}" if f.notes else "")
                    )
                total += len(files)
            print(f"{total} pxil file{'s' if total != 1 else ''}")
    except (AnnotationError, KeyError, ValueError, OSError) as exc:
        msg = exc.args[0] if isinstance(exc, KeyError) and exc.args else exc
        print(f"gpxil2pxil: error: {msg}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
