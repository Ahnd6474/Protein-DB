"""Query utilities for the SQLAlchemy protein database."""
"""Query utilities for the SQLAlchemy protein database."""

from typing import List

import numpy as np

from .database import ProteinDB, ProteinModel
from .protein import Protein


class ProteinQuery:
    """Convenience wrapper providing common protein queries."""

    def __init__(self, db: ProteinDB):
        self.Session = db.Session

    def by_accession(self, accession: str) -> List[ProteinModel]:
        """Return proteins matching a specific accession."""

        with self.Session() as session:
            return (
                session.query(ProteinModel)
                .filter(ProteinModel.accession == accession)
                .all()
            )

    def by_organism(self, organism: str) -> List[ProteinModel]:
        """Return proteins whose organism contains the given text."""

        pattern = f"%{organism}%"
        with self.Session() as session:
            return (
                session.query(ProteinModel)
                .filter(ProteinModel.organism.ilike(pattern))
                .all()
            )

    def description_contains(self, text: str) -> List[ProteinModel]:
        """Return proteins whose description contains the given text."""

        pattern = f"%{text}%"
        with self.Session() as session:
            return (
                session.query(ProteinModel)
                .filter(ProteinModel.description.ilike(pattern))
                .all()
            )

    def similar_sequence(self, sequence: str, top: int = 5) -> List[ProteinModel]:
        """Return proteins with embeddings closest to the given sequence."""

        query_protein = Protein("query", "", "", "", sequence)
        target = query_protein.Z
        with self.Session() as session:
            proteins = session.query(ProteinModel).all()
        proteins.sort(key=lambda p: np.linalg.norm(p.embedding - target))
        return proteins[:top]

    def by_embedding(
        self, sequence: str, threshold: float, metric: str = "euclidean"
    ) -> List[ProteinModel]:
        """Return proteins whose embeddings satisfy a similarity threshold.

        Parameters
        ----------
        sequence:
            Amino acid sequence to compare against the database.
        threshold:
            For ``metric='euclidean'`` this represents the maximum L2 distance.
            For ``metric='cosine'`` it is the minimum cosine similarity.
        metric:
            Either ``'euclidean'`` or ``'cosine'``.
        """

        query_protein = Protein("query", "", "", "", sequence)
        target = query_protein.Z
        with self.Session() as session:
            proteins = session.query(ProteinModel).all()

        good = []
        for p in proteins:
            emb = p.embedding
            if metric == "cosine":
                sim = float(
                    np.dot(emb, target)
                    / (np.linalg.norm(emb) * np.linalg.norm(target))
                )
                if sim >= threshold:
                    good.append(p)
            elif metric == "euclidean":
                dist = float(np.linalg.norm(emb - target))
                if dist <= threshold:
                    good.append(p)
            else:
                raise ValueError("metric must be 'euclidean' or 'cosine'")
        return good
