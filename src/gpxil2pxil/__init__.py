"""gpxil2pxil: pictures of genes, transcripts and proteins (pxil files) from a gpxil genome.

The sequences come from a gpxil file's reference, or from one of its sample's
haplotypes, so a person's own version of a protein can be drawn and compared
with the reference. Genes and transcripts come from the Gpxil browser's gene
files (Ensembl, GRCh38) or from any GTF.

Copyright (c) 2026 Chronomic Medicine Group, Department of Genomic Medicine,
University of Cambridge. MIT licence.
"""

__version__ = "0.2.0"

from .annotation import BROWSER_DATA, GTF, AnnotationError, BrowserData, Gene, Transcript, annotation  # noqa: E402
from .build import Written, make  # noqa: E402
from .sequence import Source, translate  # noqa: E402

__all__ = [
    "BROWSER_DATA",
    "AnnotationError",
    "BrowserData",
    "GTF",
    "Gene",
    "Source",
    "Transcript",
    "Written",
    "__version__",
    "annotation",
    "make",
    "translate",
]
