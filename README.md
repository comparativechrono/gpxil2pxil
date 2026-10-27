# gpxil2pxil

Pictures of genes, transcripts and proteins, as
[pxil](https://github.com/comparativechrono/pxil) files, from a
[gpxil](https://github.com/comparativechrono/Gpxil) genome: from its
reference, or from either of a person's haplotypes, so that someone's own
version of a protein can be seen and compared with the reference.

Genes and transcripts come from the
[Gpxil browser](https://comparativechrono.github.io/Gpxil-browser/)'s gene
files (every Ensembl release 116 transcript on GRCh38, read from the web a gene
at a time), or from any GTF file (Ensembl or GENCODE).

```sh
pip install gpxil2pxil

# TP53 from the GRCh38 reference, read from the Gpxil browser (nothing to download first)
gpxil2pxil https://comparativechrono.github.io/Gpxil-browser/data/grch38/gpxil/chr17.gpxil.tif TP53

# a person's TP53 and BRCA1, both haplotypes and the reference, every transcript
gpxil2pxil HG002.gpxil.tif TP53 BRCA1 --haplotypes all --transcripts all -o pictures

# the transcripts of a gene, to choose from
gpxil2pxil HG002.gpxil.tif TP53 --list

# another assembly or annotation: a GTF
gpxil2pxil genome.gpxil.tif MYGENE --annotation Homo_sapiens.GRCh38.116.gtf.gz
```

To read the Gpxil browser's gene files from a copy of them (a folder, or
another web address), set `GPXIL_BROWSER_DATA` to it, or give it as
`--annotation`.

Files are named with the Ensembl IDs, and with the sample and haplotype when
they come from one:

| File | What it holds |
|---|---|
| `TP53_ENSG00000141510.gene.pxil.png` | the gene's genomic sequence, 5′ to 3′, introns and soft-masked repeats included |
| `TP53-201_ENST00000269305.transcript.pxil.png` | the spliced transcript, 5′ to 3′ |
| `TP53-201_ENSP00000269305.protein.pxil.png` | the protein, from the start codon to the stop codon (`*`) |
| `TP53-201_ENSP00000269305_HG002-hap1.protein.pxil.png` | the same protein, made from HG002's first haplotype |

Each file's metadata (`python -m pxil --info FILE`) gives its coordinates,
exons, coding region, genetic code and, for a haplotype, the sample and how it
differs from the reference (`"notes": ["differs from the reference: P72R"]`).

## How a haplotype's sequence is made

The sample's genotypes in the gpxil file are applied to the reference, as
`gpxil fasta --which 1` does: substitutions replace bases, deletions remove
them, insertions follow their anchor base, and bases whose state is not known
become N. A transcript of a haplotype is its sequence over each exon, joined.
Its protein is translated from the annotated start codon until the first stop
codon, so a frameshift is followed into the 3′ UTR, and a protein whose stop
codon is lost ends at the transcript's end (and the file says so).
Selenocysteines (TGA read as U) are kept where Ensembl annotates them, and
mitochondrial genes use the vertebrate mitochondrial code.

## Python

```python
from gpxiltools import GpxilFile
from gpxil2pxil import BrowserData, Source, make

gene = BrowserData().find("TP53")               # or GTF("genes.gtf.gz").find("TP53")
with GpxilFile("HG002.gpxil.tif") as g:
    for f in make(g, gene, haplotypes=(None, 1, 2), out_dir="TP53"):
        print(f.path, f.notes)
    hap1 = Source(g, 1).protein(gene.canonical())  # the sequence itself, as a string
```

## Tests

```sh
pip install -e ".[test]"
pytest                                  # a made-up genome whose genes and variants are known exactly
GPXIL2PXIL_NETWORK=1 pytest             # also TP53 from the Gpxil browser
```

## Licence

Copyright © 2026 Chronomic Medicine Group, Department of Genomic Medicine,
University of Cambridge. Released under the [MIT licence](LICENSE).
