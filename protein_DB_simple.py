"""SQLAlchemy-backed protein database storing sequences and embeddings.

This module mirrors :class:`ProteinDB` from :mod:`Database` but keeps only the
protein sequence and its vector embedding.  A FAISS index is used to perform
similarity search directly on the stored embeddings.
"""

from __future__ import annotations

from typing import Iterable, List, Sequence, Tuple

import faiss
import numpy as np
from Bio import SeqIO
from sqlalchemy import Column, Integer, String, create_engine
from sqlalchemy.orm import declarative_base, sessionmaker
from sqlalchemy.types import PickleType


Base = declarative_base()


class ProteinSimpleModel(Base):
    """SQLAlchemy model representing a protein with only sequence and embedding."""

    __tablename__ = "proteins_simple"

    id = Column(Integer, primary_key=True)
    sequence = Column(String, unique=True)
    embedding = Column(PickleType)


def read_fasta(filename: str) -> List[str]:
    """Return all sequences contained in ``filename``.

    Parameters
    ----------
    filename:
        Path to a FASTA file.
    """

    return [str(record.seq) for record in SeqIO.parse(filename, "fasta")]


class ProteinDBSimple:
    """Persistent protein database backed by SQLAlchemy."""

    def __init__(self, url: str = "sqlite:///proteins_simple.db") -> None:
        self.engine = create_engine(url)
        Base.metadata.create_all(self.engine)
        self.Session = sessionmaker(bind=self.engine)

        self.sequences: List[str] = []
        self._embeddings: np.ndarray | None = None
        self._index: faiss.Index | None = None
        self._load_cache()

    # ------------------------------------------------------------------
    # Loading helpers
    # ------------------------------------------------------------------
    def _load_cache(self) -> None:
        """Load all sequences and embeddings from the database."""

        with self.Session() as session:
            rows = session.query(ProteinSimpleModel).all()
            self.sequences = [r.sequence for r in rows]
            if rows:
                self._embeddings = np.vstack(
                    [np.asarray(r.embedding, dtype=np.float32) for r in rows]
                )
            else:
                self._embeddings = None
        self._index = None  # invalidate FAISS index

    # ------------------------------------------------------------------
    # Insertion helpers
    # ------------------------------------------------------------------
    def add(self, sequence: str, embedding: Sequence[float]) -> None:
        """Insert a single sequence and embedding."""

        vec = np.asarray(embedding, dtype=np.float32)
        if vec.ndim != 1:
            raise ValueError("Embedding must be a 1D vector")

        with self.Session() as session:
            session.add(ProteinSimpleModel(sequence=sequence, embedding=vec.tolist()))
            session.commit()

        self._load_cache()

    def add_many(self, items: Iterable[Tuple[str, Sequence[float]]]) -> None:
        """Insert many sequences and embeddings at once."""

        models = []
        for seq, emb in items:
            vec = np.asarray(emb, dtype=np.float32)
            if vec.ndim != 1:
                raise ValueError("Embedding must be a 1D vector")
            models.append(ProteinSimpleModel(sequence=seq, embedding=vec.tolist()))

        with self.Session() as session:
            session.add_all(models)
            session.commit()

        self._load_cache()

    def add_from_fasta(
        self, fasta_path: str, embeddings: Iterable[Sequence[float]]
    ) -> None:
        """Read sequences from a FASTA file and insert them with embeddings.

        Parameters
        ----------
        fasta_path:
            Path to the FASTA file.
        embeddings:
            Iterable of embedding vectors matching the sequences in the FASTA
            file.
        """

        sequences = read_fasta(fasta_path)
        embedding_list = list(embeddings)
        if len(sequences) != len(embedding_list):
            raise ValueError("Number of embeddings must match number of sequences")
        self.add_many(zip(sequences, embedding_list))

    # ------------------------------------------------------------------
    # Search helpers
    # ------------------------------------------------------------------
    def _ensure_index(self) -> None:
        if self._embeddings is None:
            raise ValueError("Database is empty")
        if self._index is None:
            dim = self._embeddings.shape[1]
            self._index = faiss.IndexFlatL2(dim)
            self._index.add(self._embeddings)

    def search(self, vector: Sequence[float], k: int = 1000) -> List[str]:
        """Return the ``k`` closest sequences to ``vector``."""

        self._ensure_index()
        query = np.asarray(vector, dtype=np.float32)[None, :]
        k = min(k, len(self.sequences))
        _dists, idx = self._index.search(query, k)
        return [self.sequences[i] for i in idx[0]]


__all__ = ["ProteinDBSimple", "read_fasta"]
