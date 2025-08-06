from protein_db import ProteinDB, Protein, ProteinQuery, generate_sequences


if __name__ == "__main__":
    # Vector-based search in simple mode
    db = ProteinDB(url="sqlite:///proteins_simple.db", simple=True)
    try:
        query = Protein("query", "", "", "", "MSEQN...")
        hits = db.search(query.Z, k=5)
        print("Closest sequences:", hits)
    except Exception as exc:
        # Database may be empty if no sequences were loaded
        print(f"Vector search failed: {exc}")

    # Sequence generation using the genetic algorithm
    try:
        seqs = generate_sequences(db, query.Z, cos_threshold=0.2, rmse_threshold=0.2)
        print("Generated sequences:", seqs)
    except Exception as exc:
        print(f"Generation failed: {exc}")

    # Text-based queries on the full database
    db_full = ProteinDB(url="sqlite:///proteins.db")
    pq = ProteinQuery(db_full)
    for protein in pq.by_accession("ABC123"):
        print(protein.accession, protein.organism)
