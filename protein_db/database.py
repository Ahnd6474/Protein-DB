"""SQLAlchemy-backed protein database utilities with optional simple mode."""

from __future__ import annotations

from typing import Iterable, List, Sequence, Tuple

import faiss
import numpy as np
from Bio import SeqIO
from sqlalchemy import Column, Integer, String, create_engine
from sqlalchemy.orm import declarative_base, sessionmaker
from sqlalchemy.types import PickleType

from .protein import Protein, describe

Base = declarative_base()


class ProteinModel(Base):
    """SQLAlchemy model representing a full protein record."""

    __tablename__ = "proteins"

    id = Column(Integer, primary_key=True)
    accession = Column(String, unique=True, index=True)
    description = Column(String)
    locus = Column(String)
    organism = Column(String, index=True)
    sequence = Column(String)
    embedding = Column(PickleType)


class ProteinSimpleModel(Base):
    """SQLAlchemy model storing only sequence and embedding."""

    __tablename__ = "proteins_simple"

    id = Column(Integer, primary_key=True)
    sequence = Column(String, unique=True)
    embedding = Column(PickleType)


def read_fasta(filename: str, simple: bool = False):
    """Read a FASTA file.

    Parameters
    ----------
    filename:
        Path to the FASTA file.
    simple:
        When ``True`` return only sequences, otherwise return ``Protein``
        objects.
    """

    if simple:
        return [str(record.seq) for record in SeqIO.parse(filename, "fasta")]

    data = []
    for record in SeqIO.parse(filename, "fasta"):
        accession, description, locus, organism = describe(record.description)
        data.append(Protein(accession, description, locus, organism, str(record.seq)))
    return data


class ProteinDB:
    """High level helper around a SQLAlchemy session.

    Parameters
    ----------
    url:
        Database URL. If ``None`` it defaults to ``sqlite:///proteins.db`` or
        ``sqlite:///proteins_simple.db`` depending on ``simple``.
    simple:
        Store only sequence and embedding and enable FAISS based search.
    """

    def __init__(self, url: str | None = None, simple: bool = False) -> None:
        self.simple = simple
        if url is None:
            url = "sqlite:///proteins_simple.db" if simple else "sqlite:///proteins.db"
        self.engine = create_engine(url)
        Base.metadata.create_all(self.engine)
        self.Session = sessionmaker(bind=self.engine)

        if self.simple:
            self.sequences: List[str] = []
            self._embeddings: np.ndarray | None = None
            self._index: faiss.Index | None = None
            self._load_cache()

    # ------------------------------------------------------------------
    # Loading helpers
    # ------------------------------------------------------------------
    def _load_cache(self) -> None:
        """Load all sequences and embeddings from the database."""

        if not self.simple:
            return

        with self.Session() as session:
            rows = session.query(ProteinSimpleModel).all()
            self.sequences = [r.sequence for r in rows]
            if rows:
                self._embeddings = np.vstack(
                    [np.asarray(r.embedding, dtype=np.float32) for r in rows]
                )
            else:
                self._embeddings = None
        self._index = None

    # ------------------------------------------------------------------
    # Insertion helpers
    # ------------------------------------------------------------------
    def add(self, item, embedding: Sequence[float] | None = None) -> None:
        """Insert a single protein or sequence and embedding."""

        if self.simple:
            if embedding is None:
                raise ValueError("embedding must be provided in simple mode")
            vec = np.asarray(embedding, dtype=np.float32)
            if vec.ndim != 1:
                raise ValueError("Embedding must be a 1D vector")
            with self.Session() as session:
                session.add(ProteinSimpleModel(sequence=item, embedding=vec.tolist()))
                session.commit()
            self._load_cache()
        else:
            protein: Protein = item
            with self.Session() as session:
                session.add(
                    ProteinModel(
                        accession=protein.accession,
                        description=protein.description,
                        locus=protein.locus,
                        organism=protein.organism,
                        sequence=str(protein.sequence),
                        embedding=protein.Z,
                    )
                )
                session.commit()

    def add_many(self, items: Iterable[Tuple[str, Sequence[float]]]) -> None:
        """Insert many sequences and embeddings at once (simple mode only)."""

        if not self.simple:
            raise RuntimeError("add_many is only available in simple mode")

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
        self, fasta_path: str, embeddings: Iterable[Sequence[float]] | None = None
    ) -> None:
        """Read from a FASTA file and insert records into the database."""

        if self.simple:
            if embeddings is None:
                raise ValueError("embeddings must be provided in simple mode")
            sequences = read_fasta(fasta_path, simple=True)
            embedding_list = list(embeddings)
            if len(sequences) != len(embedding_list):
                raise ValueError(
                    "Number of embeddings must match number of sequences"
                )
            self.add_many(zip(sequences, embedding_list))
        else:
            proteins = read_fasta(fasta_path)
            with self.Session() as session:
                session.add_all(
                    [
                        ProteinModel(
                            accession=p.accession,
                            description=p.description,
                            locus=p.locus,
                            organism=p.organism,
                            sequence=str(p.sequence),
                            embedding=p.Z,
                        )
                        for p in proteins
                    ]
                )
                session.commit()

    def session(self):
        """Return a new session object (non-simple mode only)."""

        if self.simple:
            raise RuntimeError("session is only available in non-simple mode")
        return self.Session()

    # ------------------------------------------------------------------
    # Search helpers
    # ------------------------------------------------------------------
    def _ensure_index(self) -> None:
        if not self.simple:
            raise RuntimeError("search is only available in simple mode")
        if self._embeddings is None:
            raise ValueError("Database is empty")
        if self._index is None:
            dim = self._embeddings.shape[1]
            self._index = faiss.IndexFlatL2(dim)
            self._index.add(self._embeddings)

    def search(self, vector: Sequence[float], k: int = 1000) -> List[str]:
        """Return the ``k`` closest sequences to ``vector`` (simple mode)."""

        self._ensure_index()
        query = np.asarray(vector, dtype=np.float32)[None, :]
        k = min(k, len(self.sequences))
        _dists, idx = self._index.search(query, k)
        return [self.sequences[i] for i in idx[0]]


__all__ = [
    "ProteinDB",
    "ProteinModel",
    "ProteinSimpleModel",
    "read_fasta",
]
