"""Simple in-memory protein database storing sequences and embeddings.

This module provides a lightweight alternative to the SQLAlchemy-backed
:class:`ProteinDB` in :mod:`Database`. It keeps only the protein sequence and its
vector embedding and exposes a FAISS based similarity search.
"""

from __future__ import annotations

import faiss
import numpy as np
from typing import Iterable, List, Sequence, Tuple


class ProteinDBSimple:
    """In-memory collection of protein sequences and their embeddings."""

    def __init__(self) -> None:
        self.sequences: List[str] = []
        self._embeddings: np.ndarray | None = None
        self._index: faiss.Index | None = None

    def add(self, sequence: str, embedding: Sequence[float]) -> None:
        """Add a single protein to the database.

        Parameters
        ----------
        sequence:
            Amino acid sequence.
        embedding:
            Vector embedding for the sequence. Must be 1‑dimensional.
        """

        vec = np.asarray(embedding, dtype=np.float32)
        if vec.ndim != 1:
            raise ValueError("Embedding must be a 1D vector")

        self.sequences.append(sequence)
        if self._embeddings is None:
            self._embeddings = vec[None, :]
        else:
            self._embeddings = np.vstack([self._embeddings, vec])

        self._index = None  # invalidate FAISS index

    def add_many(self, items: Iterable[Tuple[str, Sequence[float]]]) -> None:
        """Add many proteins at once."""

        for seq, emb in items:
            self.add(seq, emb)

    def _ensure_index(self) -> None:
        if self._embeddings is None:
            raise ValueError("Database is empty")
        if self._index is None:
            dim = self._embeddings.shape[1]
            self._index = faiss.IndexFlatL2(dim)
            self._index.add(self._embeddings)

    def search(self, vector: Sequence[float], k: int = 1000) -> List[str]:
        """Return the ``k`` closest sequences to ``vector``.

        The search is performed in the embedding space using a FAISS index.  If
        the index has not been built yet it will be constructed on first use.

        Parameters
        ----------
        vector:
            Query embedding vector.
        k:
            Number of results to return (default 1000).
        """

        self._ensure_index()
        query = np.asarray(vector, dtype=np.float32)[None, :]
        k = min(k, len(self.sequences))
        _dists, idx = self._index.search(query, k)
        return [self.sequences[i] for i in idx[0]]


__all__ = ["ProteinDBSimple"]
