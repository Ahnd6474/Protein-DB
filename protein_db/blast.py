"""BLAST utilities for protein sequences."""

from typing import List, Tuple

from Bio.Blast import NCBIWWW

from .database import ProteinDB
from .query import ProteinQuery


def protein_blast(sequence: str, program: str = "blastp", database: str = "nr") -> str:
    """Run an NCBI BLAST query for a single protein sequence.

    Parameters
    ----------
    sequence:
        Amino acid sequence to query.
    program:
        BLAST program to use (default: ``blastp``).
    database:
        Target database (default: ``nr``).

    Returns
    -------
    str
        Raw BLAST result as a string.
    """
    handle = NCBIWWW.qblast(program=program, database=database, sequence=sequence)
    return handle.read()


def deep_blast(
    sequence: str,
    db: ProteinDB,
    threshold: float,
    metric: str = "euclidean",
    program: str = "blastp",
    database: str = "nr",
) -> List[Tuple[str, str]]:
    """Filter proteins by embedding similarity and BLAST the matches.

    Parameters
    ----------
    sequence:
        Query amino acid sequence.
    db:
        ``ProteinDB`` instance to search.
    threshold:
        Distance or similarity cutoff for latent vectors.
    metric:
        ``euclidean`` for L2 distance or ``cosine`` for cosine similarity.
    program:
        BLAST program to use.
    database:
        BLAST target database.

    Returns
    -------
    list of tuple
        Each tuple contains the accession and the raw BLAST result for a
        matching sequence.
    """
    query = ProteinQuery(db)
    good = query.by_embedding(sequence, threshold, metric)

    results: List[Tuple[str, str]] = []
    for protein in good:
        results.append(
            (protein.accession, protein_blast(protein.sequence, program, database))
        )
    return results
