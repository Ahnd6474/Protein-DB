"""Query utilities for the SQLAlchemy protein database."""

from typing import List

import numpy as np

from Database import ProteinDB, ProteinModel
from Protein import Protein


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
