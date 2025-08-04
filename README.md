# Protein-DB

A lightweight toolkit for building a protein sequence database and representing proteins with embeddings produced by a variational autoencoder (VAE).

## Installation

Clone the repository and install the dependencies:

```bash
pip install -r requirements.txt
```

## Usage

### Working with proteins

```python
from Protein import Protein
protein = Protein(accession="ABC123", description="Sample protein", locus="locus1", organism="E. coli", seq="MSEQN...")
print(len(protein))          # Sequence length
print(protein.Z.shape)       # Latent representation from the VAE
```

### Loading a FASTA file into a database

```python
from Database import ProteinDB

db = ProteinDB(fasta="proteins.fasta")
print(len(db))
```

The repository includes helper utilities for parsing FASTA files and computing sequence embeddings. More advanced querying functionality can be added via the `Query` module.

