"""SQLAlchemy-backed protein database utilities."""

from Bio import SeqIO
from sqlalchemy import Column, Integer, String, create_engine
from sqlalchemy.orm import declarative_base, sessionmaker
from sqlalchemy.types import PickleType

from Protein import Protein, describe

Base = declarative_base()


class ProteinModel(Base):
    """SQLAlchemy model representing a protein."""

    __tablename__ = "proteins"

    id = Column(Integer, primary_key=True)
    accession = Column(String, unique=True, index=True)
    description = Column(String)
    locus = Column(String)
    organism = Column(String, index=True)
    sequence = Column(String)
    embedding = Column(PickleType)


def read_fasta(filename):
    """Read a FASTA file and return a list of Protein objects."""

    data = []
    for record in SeqIO.parse(filename, "fasta"):
        accession, description, locus, organism = describe(record.description)
        data.append(Protein(accession, description, locus, organism, str(record.seq)))
    return data


class ProteinDB:
    """High level helper around a SQLAlchemy session."""

    def __init__(self, url: str = "sqlite:///proteins.db"):
        self.engine = create_engine(url)
        Base.metadata.create_all(self.engine)
        self.Session = sessionmaker(bind=self.engine)

    def add_protein(self, protein: Protein) -> None:
        """Insert a single Protein into the database."""

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

    def add_from_fasta(self, fasta_path: str) -> None:
        """Read proteins from a FASTA file and insert them into the DB."""

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
        """Return a new session object."""

        return self.Session()

