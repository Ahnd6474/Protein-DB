# Protein-DB

A lightweight toolkit for building a protein sequence database and representing proteins with embeddings produced by a variational autoencoder (VAE).

## Installation

Clone the repository and install the dependencies:

```bash
pip install -r requirements.txt
```

## Usage


## Dataset & Embeddings

### Source dataset (UniRef50, 1M subset)
The sample data in this repository is derived from **UniRef50**, using an approximately **1 million sequence subset**. This subset is intentionally lightweight to make vector indexing, nearest‑neighbor search, and quick demos feasible on a single machine.

### Embedding generation (LATTE encoder)
All sequences were passed through the **LATTE** paper’s **encoder** to produce **latent embeddings** (e.g., 256‑dim). The encoder was trained with structural/perceptual objectives so that the latent space reflects coarse structural similarity in addition to sequence patterns. These embeddings can be stored in a vector database (e.g., FAISS) for fast similarity search and downstream analysis.

### Intended use: “Deep BLAST” (latent prefilter → alignment)
As outlined in the LATTE paper, a practical extension is a two‑stage search workflow:
1) **Prefilter in latent space** (e.g., cosine similarity) to quickly shortlist candidates, then  
2) **Run alignment (e.g., BLAST)** only on the shortlisted sequences to keep interpretability and precision.
This can reduce fan‑out and wall‑clock time substantially on large corpora while preserving alignment‑based interpretability.

### Reproducibility notes
- Re‑generate embeddings by running the **LATTE encoder** on the same sequences with a fixed model/version.  
- For large‑scale search, store only the **(sequence, embedding)** pairs and index the embeddings with your preferred ANN library.

### Working with proteins

```python
from protein_db import Protein
protein = Protein(accession="ABC123", description="Sample protein", locus="locus1", organism="E. coli", seq="MSEQN...")
print(len(protein))          # Sequence length
print(protein.Z.shape)       # Latent representation from the VAE
```

### Loading a FASTA file into a database

```python
from protein_db import ProteinDB

db = ProteinDB()
db.add_from_fasta("proteins.fasta")
```

Set ``simple=True`` to store only sequences and embeddings for fast vector search.
The repository includes helper utilities for parsing FASTA files and computing sequence embeddings. More advanced querying functionality can be added via the `ProteinQuery` class.

### BLAST searches

Two helper functions simplify running BLAST queries.

```python
from protein_db import protein_blast, deep_blast, ProteinDB

# Run a plain BLASTP search
result = protein_blast("MSEQN...")

# Filter by latent vector similarity then BLAST the matches
db = ProteinDB("sqlite:///proteins.db")
hits = deep_blast("MSEQN...", db, threshold=0.5, metric="cosine")
```

The `threshold` and `metric` arguments control how candidate sequences are
selected. Use `metric="euclidean"` for an L2 distance cut-off or
`metric="cosine"` for a cosine similarity threshold.

### Sequence generation

Generate candidate sequences that approach a target embedding using a simple
genetic algorithm.

```python
from protein_db import Protein, ProteinDB, generate_sequences

db = ProteinDB(url="sqlite:///proteins_simple.db", simple=True)
target_vec = Protein("tmp", "", "", "", "MSEQN...").Z
seqs = generate_sequences(db, target_vec, cos_threshold=0.2, rmse_threshold=0.2)
print(seqs)  # Five optimized sequences
```

The `cos_threshold` and `rmse_threshold` arguments control acceptance of
mutations during evolution.

### Visualization and UI

A simple Streamlit application provides a graphical interface for searching
the database and visualizing embedding results.

```bash
streamlit run app.py
```

Enter a query sequence, choose a distance metric and threshold, and the app
will display a PCA scatter plot of the matching embeddings. Optionally, the
interface can run BLAST on the hits and show the raw results.

The plotting utility can also be used programmatically:

```python
from protein_db import plot_embeddings

fig = plot_embeddings(hits, "MSEQN...")
fig.show()
```

